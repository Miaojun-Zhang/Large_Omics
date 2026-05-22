
import sys
import csv
import asyncio
import time
import vcfpy
from fake_enformer import async_predict

"""
Concurrency Approach: asyncio + async_predict + Semaphore

Why this approach:
- fake_enformer provides async_predict natively — no need for HTTP or multiprocessing
- The task is purely I/O-bound (5-20s sleep per call), NOT CPU-bound
  Multiprocessing would waste process-spawn overhead
  Asyncio has near-zero overhead and handles thousands of concurrent waits easily
- A Semaphore(50) caps simultaneous in-flight requests to be a polite neighbour
  on the shared server, as advised in the docs

Expected speedup:
- Sequential: 3400 SNPs × avg 12.5s ≈ 42500s (~12 hours)
- Async(50):  3400 SNPs / 50 concurrent × avg 12.5s ≈ 850s (~14 min)
- Speedup factor: ~50×
"""

# Tune this: 10-50 is reasonable for the shared server
CONCURRENCY = 50
   
def read_snps(vcf_path):
    """Read and deduplicate SNP coordinates from a VCF file."""
    reader = vcfpy.Reader.from_path(vcf_path)
    seen = set()
    snps = []
    for record in reader:
        #record = each variant
        #read each line of the VCF 
        #get the first alternate 
        alt = record.ALT[0]
        # Format expected by fake_enformer:
        # hg<build>:chr<N|X|Y>:<position>:<REF>:<ALT>
        if len(record.REF) != 1 or len(alt.value) != 1:
            continue
        snp = f"hg38:{record.CHROM}:{record.POS}:{record.REF}:{alt.value}"
        if snp not in seen:
            seen.add(snp)
            snps.append(snp)
  
    return snps

async def fetch_one(semaphore, snp, index, total):
    """Fetch score for one SNP, bounded by semaphore."""
    async with semaphore:
        print(f"[{index}/{total}] Fetching: {snp}")
        score = await async_predict(snp)
        print(f"[{index}/{total}] Done:     {snp} -> {score:.4f}")
        return (snp, score)


async def run(vcf_path, concurrency):
    snps = read_snps(vcf_path)
    total = len(snps)
    print(f"Found {total} unique SNPs.")
    print(f"Concurrency = {concurrency}")

    semaphore = asyncio.Semaphore(concurrency)

    start = time.time()

    tasks = [
        fetch_one(semaphore, snp, i + 1, total)
        for i, snp in enumerate(snps)
    ]
    results = await asyncio.gather(*tasks)

    elapsed = time.time() - start
    avg = elapsed / total if total else 0
    print(f"\nFinished {total} SNPs in {elapsed:.1f}s (avg {avg:.2f}s/SNP)")
    print(f"Theoretical sequential time: {total * 12.5:.0f}s")
    print(f"Speedup estimate: ~{total * 12.5 / elapsed:.1f}×")

    #Write to TSV file
    output_file = f"enformer_results_c{concurrency}.tsv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        # write the header
        writer.writerow(['coordinate', 'score'])
        #write the result returned from parallelism
        writer.writerows(results)
    
    print(f'Done! Results saved to {output_file}')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python run_enformer.py <path_to_vcf>")
        sys.exit(1)
    concurrency = int(sys.argv[2]) if len(sys.argv) > 2 else CONCURRENCY
    asyncio.run(run(sys.argv[1], concurrency))
 
