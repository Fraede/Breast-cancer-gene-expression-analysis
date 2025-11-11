# -*- coding: utf-8 -*-
"""
Created on Sat Nov  8 01:52:07 2025

@author: reedt
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss 
from statsmodels.stats.multitest import multipletests
#read dataset without metadata
file_path = r"C:\Users\reedt\Downloads\Breat Cancer Analysis\GSE16873_series_matrix.txt"
df = pd.read_csv(file_path, sep='\t', comment='!', index_col=0)

print("Data shape:", df.shape)
print(df.head())
 
#extract important metadata
sample_states =[]
with open(file_path) as file:
    for line in file:
        if line.startswith('!Sample_characteristics_ch1	"disease state:'):
            states = line.split('"disease state:')
            for item in states[1:]:
                sample_states.append(item.strip().replace('"', ''))

print("Unique states found:", sample_states[:4])

#Match unique groups to each column
sample_names = list(df.columns)
normal =[]
simple_dh = []
a_dh = []
dcis = []
for i,state in enumerate(sample_states):
    if state == "histologically normal":
        normal.append(sample_names[i])
    elif state == "simple ductal hyperplasia":
        simple_dh.append(sample_names[i])
    elif state == "atypical ductal hyperplasia":
        a_dh.append(sample_names[i])
    elif state == "ductal carcinoma in situ":
        dcis.append(sample_names[i])
print(f"Normal: {len(normal)}, Simple DH: {len(simple_dh)}, ADH: {len(a_dh)}, DCIS: {len(dcis)}")
        
#Compute means for each group

mean_normal = df.loc[:, normal].mean(axis=1)
mean_simple_dh = df.loc[:, simple_dh].mean(axis=1)
mean_a_dh = df.loc[:, a_dh].mean(axis=1)
mean_dcis = df.loc[:, dcis].mean(axis=1)

#Calculate log2 flod change
eps = 1e-6  # to avoid division by zero
fc_dcis_normal = np.log2((mean_dcis + eps) / (mean_normal + eps))

#perform t-test for each gene
t_stat, p_values = ss.ttest_ind(df.loc[:, dcis], df.loc[:, normal], axis=1, equal_var=False)

#conduct fdr correction
_, p_values_fdr, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

#Tabulate results into a new data frame
results = pd.DataFrame({"log2FoldChange":fc_dcis_normal, "p values":p_values, "p_value_FDR" : p_values_fdr}, index = df.index)


#define thresholds for significant genes
fc_dcis_normal_thold = 1.00             #>= 2-fold change
pval_thold = 0.05

# create a column marking significance
sig = (abs(results["log2FoldChange"]) > fc_dcis_normal_thold) & (results["p values"] < pval_thold)
results['significant'] = sig

#check how many genes passed
sig_count = results['significant'].sum()
print(f"Significant genes: {sig_count} / {len(results)}")


# Scatter plot: all genes
plt.figure(figsize=(7,5)) 
plt.scatter(results['log2FoldChange'], -np.log10(results['p values']), s=8, color='grey', alpha=0.5)

# Highlight significant genes
sig_genes = results[results['significant']]
plt.scatter(sig_genes['log2FoldChange'], -np.log10(sig_genes['p values']), s=8, color='red')

# Add threshold lines
plt.axvline(x=fc_dcis_normal_thold, color='blue', linestyle='--')
plt.axvline(x=-fc_dcis_normal_thold, color='blue', linestyle='--')
plt.axhline(y=-np.log10(pval_thold), color='green', linestyle='--')

plt.xlabel('log2 Fold Change (DCIS / Normal)')
plt.ylabel('-log10(p values)')
plt.title('Volcano Plot — DCIS vs Normal')
plt.tight_layout()
plt.show()


#Read file containing gene IDs
ID_path = r"C:\Users\reedt\Downloads\Breat Cancer Analysis\GPL570-55999.txt"
gpl = pd.read_csv(ID_path, sep='\t', comment='#', low_memory=False)
gpl = gpl[['ID', 'Gene Symbol', 'Gene Title']]
gpl.head()

#Merge with results
annotated_results = results.merge(gpl, left_index=True, right_on='ID', how='left')
print(annotated_results.sort_values('p values')[['log2FoldChange', 'p values', 'Gene Symbol']])

#Show progression for top genes.
top10_symbols = annotated_results.sort_values('p values').head(10)['Gene Symbol'].tolist()
top10_genes = results.sort_values('p values').head(10).index.tolist()

gene_lis = [] 
for i in range(len(top10_genes)):
    gene_lis.append((top10_genes[i], top10_symbols[i]))
gene_map = dict(gene_lis)


stages = ['Histologically Normal', 'Simple DH', 'ADH', 'DCIS']

plt.figure(figsize=(9,6))
for gene in top10_genes:
    means = [mean_normal[gene],mean_simple_dh[gene], mean_a_dh[gene], mean_dcis[gene]]
    label = f"{gene_map[gene]} ({gene})"
    plt.plot(stages, means, marker='o', linewidth=2, alpha=0.8, label=label)

plt.title('Expression Progression of Top 10 Significant Genes')
plt.xlabel('Breast Tissue Stage')
plt.ylabel('Mean Expression Level')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
plt.tight_layout()
plt.show()

    






























