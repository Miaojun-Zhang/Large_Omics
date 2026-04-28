import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

output_dir = "160.figures"
os.makedirs(output_dir, exist_ok=True)

db_path = snakemake.input[0]
conn = sqlite3.connect(db_path)

# Figure 1: SNP impact severity per sample
query1 = '''
    SELECT c.sample_name, e.impact, COUNT(DISTINCT s.snp_id) as snp_count
    FROM snps s
    JOIN effect e ON s.snp_id = e.snp_id
    JOIN call c ON s.snp_id = c.snp_id
    WHERE c.genotype != (s.ref || '/' || s.ref)
    GROUP BY c.sample_name, e.impact
'''

#pd.read_sql_query：Read SQL query into a DataFrame
df1 = pd.read_sql_query(query1, conn)

#categorize the types of impact
impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
df1['impact'] = pd.Categorical(df1['impact'], categories=impact_order, ordered=True)
df1['sample_name'] = df1['sample_name'].map({'TLE66_N': 'Normal', 'TLE66_T': 'Tumour'})

#create a canvas with 10X6 inches2
plt.figure(figsize=(10, 6))
#set background as whitegrid
sns.set_theme(style="whitegrid")
#plot the bar chart
sns.barplot(data=df1, 
            x='impact',
            y='snp_count', 
            hue='sample_name',  #distinguish samples by colors
            palette='colorblind')

plt.title('Figure 1: SNP Impact Severity per Sample', fontsize=14, fontweight='bold')
plt.xlabel('Impact Category', fontsize=12)
plt.ylabel('Number of Distinct SNPs', fontsize=12)
plt.ylim(0, None)  # 坐标轴不截断

plt.tight_layout()
plt.savefig(f'{output_dir}/fig1_impact_severity.svg')
plt.close()

# Figure 2: Top 10 genes affected by HIGH or MODERATE variants
# Justification: Identifying most affected genes helps pinpoint potential functional 
# consequences and key drivers in the tumour sample.
query2 = '''
    SELECT e.gene, COUNT(DISTINCT s.snp_id) as snp_count
    FROM snps s
    JOIN effect e ON s.snp_id = e.snp_id
    WHERE e.impact IN ('HIGH', 'MODERATE')
    GROUP BY e.gene
    ORDER BY snp_count DESC
    LIMIT 10
'''
df2 = pd.read_sql_query(query2, conn)

plt.figure(figsize=(10, 6))
sns.barplot(data=df2, x='snp_count', y='gene', palette='magma')

plt.title('Figure 2: Top 10 Genes with High/Moderate Impact Variants', fontsize=14, fontweight='bold')
plt.xlabel('Number of Unique SNPs', fontsize=12)
plt.ylabel('Gene Symbol', fontsize=12)

plt.tight_layout()
plt.savefig(f'{output_dir}/fig2_top_affected_genes.svg')
plt.close()

#Figure 3: Distribution of SNP quality scores (QUAL) across impact categories
# Justification: This violin plot shows the confidence level of variant calls 
# across different functional impact groups, helping to assess if high-impact 
# variants meet the necessary quality thresholds.
query3 = '''
    SELECT s.qual, e.impact
    FROM snps s
    JOIN effect e ON s.snp_id = e.snp_id
'''
#create df for qual but drop NA
df3 = pd.read_sql_query(query3, conn).dropna(subset=['qual'])
df3['impact'] = pd.Categorical(df3['impact'], categories=impact_order, ordered=True)

plt.figure(figsize=(10, 6))
sns.set_theme(style="whitegrid")
#plot kdeplot-Kernel Density Estimate Plot
# common_norm=False: Calculate the density of each curve independently 
# (to prevent the HIGH curve, which has a small amount of data, from being flattened and invisible).
# clip=(0, 500): set cutoff
# hue='impact': plot 4 diferent-color lines 
sns.kdeplot(data=df3, x='qual', hue='impact', palette='colorblind', 
            fill=True, common_norm=False, alpha=0.3, linewidth=2, clip=(0, 500))

plt.title('Figure 3: Density of SNP Quality Scores by Impact Category', fontsize=14, fontweight='bold')
plt.xlabel('Quality Score (QUAL)', fontsize=12) 
plt.ylabel('Density (Frequency)', fontsize=12)
plt.xlim(0, 520)

plt.tight_layout()
plt.savefig(f'{output_dir}/fig3_qual_density.svg')
plt.close()

conn.close()