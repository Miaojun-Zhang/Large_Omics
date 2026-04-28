import os
import pandas as pd
import sqlite3
import vcfpy

os.makedirs("150.db", exist_ok=True)
VCF = snakemake.input[0]
MAX_ALLELE_LENGTH = 20

# create database and make a connection to the database 
conn = sqlite3.Connection('150.db/snps.sqlite')

# use vcfpy to open the vcf file, but didn't read the file, 
#reader needs to work with for loop to get the data
reader = vcfpy.Reader.from_path(VCF)
# a list to 
snps = []
effects = []
calls = []

for record in reader:
#for each loop, read one line(record) of the file
#each record.calls gives per-sample genotype
    
    assert len(record.REF) <= MAX_ALLELE_LENGTH
    assert len(record.ALT) == 1
    
    #acquire info related to SNPs table
    alt = record.ALT[0].value
    snp_id = "_".join([record.CHROM, str(record.POS), record.REF, alt])
    snps.append([
        snp_id, record.CHROM, record.POS, record.REF, alt])
    
    #acquire info related to Effect table
    #INFO in record is a dictionary, and get the value by the key-"ANN"
    #If ANN can't be found, return empty list
    ann_list = record.INFO.get('ANN',[])
    for entry in ann_list:
        parts = entry.split('|')
        effect = parts[1]
        impact = parts[2]
        gene = parts[3]
        effects.append([snp_id, gene, effect, impact])
    
     #acquire info related to call table
    for call in record.calls:
        sample_name = call.sample
        sample_name = sample_name.split('/')[-1].split('.')[0]
        bases = call.gt_bases
        genotype = '/'.join(str(b) for b in bases) if bases else None
        parts = sample_name.split('_')
        phenotype = parts[1]
        calls.append([snp_id,sample_name,genotype,phenotype])

reader.close()
#First table: SNP.     1    snp_id, chr, pos, ref, alt  ..... (genome build)
df_snps = pd.DataFrame(snps, columns='snp_id chrom pos ref alt'.split())
df_snps

#Second table:effect.  +    snp_id, gene, effect, impact, ...
df_effect = pd.DataFrame(effects, columns= 'snp_id gene effect impact'.split())
df_effect

#Third table: call.    +    snp_id, sample_name, phenotype, patient, ...
df_call = pd.DataFrame(calls, columns= 'snp_id sample_name genotype phenotype'.split())
df_call

#write dataframe into the database 
df_snps.to_sql('snps', conn, if_exists='replace', index=False)
df_effect.to_sql('effect', conn, if_exists='replace', index=False)
df_call.to_sql('call', conn, if_exists='replace', index=False)

conn.close()

# SQL query
conn = sqlite3.Connection('150.db/snps.sqlite')

#1.count the number of SNPs
number_snps = conn.execute('select Count(*) from snps').fetchone()[0]
print(f"Number of SNPs: {number_snps}")

#2.check that all calls have a unknown genotype
check_calls = conn.execute('''select * from call
                           where genotype IS NULL 
                           OR genotype = './.'
                           OR genotype = 'None/None'
                           ''').fetchall()
print(f"Calls with unknown genotype: {check_calls}")

#2.count the number of calls whose phenotype is 'T'
t_calls = conn.execute('''select Count(*)from call
                           where phenotype = 'T' 
                           ''').fetchone()[0]
print(f"Number of T calls: {t_calls}")

#To verify that every effect links to a known SNP
verify_links = conn.execute('''select * from effect
                           where snp_id  NOT IN (select snp_id from snps)
                           ''').fetchall()
print(f"Effects with unknown SNP: {verify_links}")

conn.close()