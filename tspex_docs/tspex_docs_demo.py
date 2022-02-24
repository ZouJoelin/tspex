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



# General scoring metrics : roku...
tso_roku = tspex.TissueSpecificity(expression_data, 'roku_specificity', log=True)
tso_roku.expression_data
tso_roku.tissue_specificity

tso_roku.plot_histogram()
tso_roku.plot_heatmap(threshold=0.8, sort_genes=True, use_zscore=True, gene_names=False)

tso_roku = tspex.TissueSpecificity(expression_data, 'roku_specificity', log=True, transform=False)
tso_roku.plot_histogram()



# Individualized scoring metric : spm...
tso_spm = tspex.TissueSpecificity(expression_data, 'spm', log=True)
tso_spm.expression_data
tso_spm.tissue_specificity

sns.violinplot(x="variable", y="value", inner=None, scale="count", color='C0', data=tso_spm.tissue_specificity.melt())
sns.clustermap(tso_spm.expression_data, figsize=(6,6), z_score=0)


