# -- coding:utf-8 --
# tspex: https://apcamargo.github.io/tspex/python_api/

import tspex
import pandas as pd
import seaborn as sns

# Data poccessing...
gtex_link = 'https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz'
expression_data = pd.read_csv(gtex_link, sep='\t', index_col='gene_id', skiprows=2)
expression_data = expression_data.loc[:, ['Bladder', 'Liver', 'Lung', 'Pancreas', 'Stomach']]
expression_data = expression_data.loc[(expression_data > 0).any(axis=1)]
expression_data

# General scoring metrics : 
##  tau...
tso_tau = tspex.TissueSpecificity(expression_data, 'tau', log=True)
tso_tau.expression_data
tso_tau.tissue_specificity

tso_tau_plot = tso_tau.plot_histogram()
tso_tau_heatmap = tso_tau.plot_heatmap(threshold=0.8, sort_genes=True, use_zscore=True, gene_names=False)

##  roku...
tso_roku = tspex.TissueSpecificity(expression_data, 'roku_specificity', log=True)
tso_roku.expression_data
tso_roku.tissue_specificity

tso_roku.plot_histogram()
tso_roku.plot_heatmap(threshold=0.8, sort_genes=True, use_zscore=True, gene_names=False)

tso_roku = tspex.TissueSpecificity(expression_data, 'roku_specificity', log=True, transform=False)
tso_roku.plot_histogram()


# Individualized scoring metric : 
## spm...
tso_spm = tspex.TissueSpecificity(expression_data, 'spm', log=True)
tso_spm.expression_data
tso_spm.tissue_specificity

sns.violinplot(x="variable", y="value", inner=None, scale="count", color='C0', data=tso_spm.tissue_specificity.melt())

###selecting genes that are preferentially expressed in the bladder and lung
selected_genes_specificity = tso_spm.tissue_specificity.query('Bladder >= 0.7 & Lung >= 0.7')
selected_genes_expression = tso_spm.expression_data.reindex(selected_genes_specificity.index)
selected_genes_specificity
selected_genes_expression
sns.clustermap(selected_genes_expression, figsize=(6,6), z_score=0);


# The threshold parameter of the counts metric : dafault 0
tso_counts = tspex.TissueSpecificity(expression_data, 'counts', threshold=50)
sum(tso_counts.tissue_specificity >= 0.8)
tso_counts.plot_heatmap(0.8, sort_genes=True, use_zscore=True, gene_names=False)




