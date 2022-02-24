# tspex
---
## 概述

[tspex](https://apcamargo.github.io/tspex/)是一个计算基因在组织表达特异性的工具。为使用者提供了Python API，命令行接口(CLI)，以及[网页版](https://tspex.lge.ibi.unicamp.br)三种方式来对基因表达数据进行 tissue-specificity 指标计算。


!!! warning ""
      tspex是基于GPL-3.0的开源软件，源代码可从[GitHub repository](https://github.com/apcamargo/tspex/)获取，以及[网页版源代码](https://github.com/apcamargo/tspex-webapp/)。



## Tissue-specificity 指标

tspex提供了12个不同的Tissue-specificity 指标，使用者可根据需求进行选择。这12个指标可大致被分为两类：

   * **-General scoring metrics：** 总体上衡量一个基因在所有组织中的表达特异性水平（`Counts`,`Tau`,`Gini coefficient`,`Simpson index`,`Shannon entropy specificity`,`ROKU specificity`,`Specificity measure dispersion`,`Jensen-Shannon specificity dispersion`）
   * **-Individualized scoring metrics：** 对基因在每个组织中的表达特异性水平进行量化（`Tissue-specificity index`,`Z-score`,`Specificity measure`,`Jensen-Shannon specificity`）


### -General scoring metrics

#### Counts

$$
\begin{align}
    & t(x_i) =
    \begin{cases}
        1, & \text{if } x_i > threshold \\
        0, & \text{if } x_i \le threshold
    \end{cases} \\[9pt]
    & \operatorname{Counts} = \frac{n-\sum_{i=1}^n(t(x_i))}{n-1}
\end{align}
$$

#### Tau

$$
\begin{align}
    & \widehat{x_i} = \frac{x_i}{\max\limits_{0\le i \le n}(x_i)} \\[9pt]
    & \operatorname{Tau} = \frac{\sum_{i=1}^n(1-\widehat{x_i})}{(n-1)}
\end{align}
$$

#### Gini coefficient

$$
\begin{align}
    & X = (x_1,x_2,\ldots,x_i,\ldots,x_{n-1},x_{n}) \\[9pt]
    & \qquad \;\; x_1 \le x_2 \le \ldots x_i \le \ldots \le x_{n-1} \le x_{n} \\[9pt]
    & \operatorname{Gini} = \frac{\sum_{i=1}^n(2i-n-1)x_i}{n\sum_{i=1}^n(x_i)} \\[9pt]
    & \operatorname{Gini'} = \operatorname{Gini} \frac{n}{n-1}
\end{align}
$$

#### Simpson index

$$
\begin{align}
    & p_i = \frac{x_i}{\sum_{i=1}^n(x_i)} \\[9pt]
    & \operatorname{Simpson} = \sum_{i=1}^n(p_i^2) \\[9pt]
    & \operatorname{Simpson'} = \frac{\operatorname{Simpson} - \frac{1}{n}}{1 - \frac{1}{n}}
\end{align}
$$

#### Shannon entropy specificity (HS)

$$
\begin{align}
    & p_i = \frac{x_i}{\sum_{i=1}^n(x_i)} \\[9pt]
    & H = -\sum_{i=1}^n(p_i \log_2(p_i)) \\[9pt]
    & \operatorname{HS} = log_2(n) - H \\[9pt]
    & \operatorname{HS'} = \frac{\operatorname{HS}}{log_2(n)}
\end{align}
$$

#### ROKU specificity

$$
\begin{align}
    & M = \operatorname{median}_{0\le i \le n}(x_i) \\[9pt]
    & S = \operatorname{median}_{0\le i \le n}(|x_i - M|) \\[9pt]
    & u_i = \frac{x_i - M}{5S + 10^{-4}} \\[9pt]
    & w(u_i) =
    \begin{cases}
        (1-u^2)^2, & \text{if } |u_i| \le 1 \\
        0, & \text{if } |u_i| > 1
    \end{cases} \\[9pt]
    & t = \frac{\sum_{i=1}^n(x_i w(u_i))}{\sum_{i=1}^n(w(u_i))} \\[9pt]
    & X' = (|x_1-t|,|x_2-t|,\ldots,|x_i-t|,\ldots,|x_{n-1}-t|,|x_{n}-t|) \\[9pt]
    & p_i = \frac{x_i'}{\sum_{i=1}^n(x_i')} \\[9pt]
    & H = -\sum_{i=1}^n(p_i \log_2(p_i)) \\[9pt]
    & \operatorname{ROKU} = log_2(n) - H \\[9pt]
    & \operatorname{ROKU'} = \frac{\operatorname{ROKU}}{log_2(n)}
\end{align}
$$

#### Specificity measure dispersion (SPM DPM)

$$
\begin{align}
    & \operatorname{\overline{SPM}} = \frac{\sum_{i=1}^n(\operatorname{SPM_i})}{n} \\[9pt]
    & \sigma = \sqrt{\frac{\sum_{i=1}^n(\operatorname{SPM_i} - \operatorname{\overline{SPM}})^2}{n-1}} \\[9pt]
    & \operatorname{SPM\ DPM} = \sigma \sqrt{n}
\end{align}
$$

#### Jensen-Shannon specificity dispersion (JSS DPM)

$$
\begin{align}
    & \operatorname{\overline{JSS}} = \frac{\sum_{i=1}^n(\operatorname{JSS_i})}{n} \\[9pt]
    & \sigma = \sqrt{\frac{\sum_{i=1}^n(\operatorname{JSS_i} - \operatorname{\overline{JSS}})^2}{n-1}} \\[9pt]
    & \operatorname{JSS\ DPM} = \sigma \sqrt{n}
\end{align}
$$


<br>

###  -Individualized scoring metrics

#### Tissue-specificity index (TSI)

$$
\begin{align}
    & \operatorname{TSI_i} = \frac{x_i}{\sum_{i=1}^n(x_i)}
\end{align}
$$

#### Z-score

$$
\begin{align}
    & \overline{x} = \frac{\sum_{i=1}^n(x_i)}{n} \\[9pt]
    & \sigma = \sqrt{\frac{\sum_{i=1}^n(x_i - \overline{x})^2}{n-1}} \\[9pt]
    & \operatorname{Z-Score_i} = \frac{x_i - \overline{x}}{\sigma} \\[9pt]
    & \operatorname{Z-Score_i'} = \frac{\operatorname{Z-Score_i} + \frac{(n-1)}{\sqrt{n}}}{2 \frac{(n-1)}{\sqrt{n}}}
\end{align}
$$

#### Specificity measure (SPM)

$$
\begin{align}
    & X = (x_1,x_2,\ldots,x_i,\ldots,x_{n-1},x_{n}) \\[9pt]
    & \operatorname{SPM_i} = \frac{x_i^2}{x_i*||X||_2}
\end{align}
$$

#### Jensen-Shannon specificity (JSS)

$$
\begin{align}
    & X = (x_1,x_2,\ldots,x_i,\ldots,x_{n-1},x_{n}) \\[9pt]
    & X' = (0,0,\ldots,x_i,\ldots,0,0) \\[9pt]
    & p_i = \frac{x_i}{\sum_{i=1}^n(x_i)} \\[9pt]
    & q_i = \frac{x'_i}{\sum_{i=1}^n(x'_i)} \\[9pt]
    & \operatorname{JSS_i} = 1 - \sqrt{\frac{\sum_{i=1}^n(p_i \log_2(p_i))}{2} - \sum_{i=1}^n\left(\frac{p_i+q_i}{2} \log_2\left(\frac{p_i+q_i}{2}\right)\right)}
\end{align}
$$

---

## Python API 使用教程

### 安装API

有以下两种方式安装tspex包：
1. Using pip
```
pip install tspex
```
2. Using conda
```
conda install -c conda-forge -c bioconda tspex
```
还将会用到pandas, seaborn库。在这里推荐直接使用Anaconda平台，已经内置了许多常用的第三方库，其中就已涵括上述两个库。当然也可以自行选择用pip逐个安装，只是可能出现兼容性等各种问题。

!!! warning ""
      **注意：** Anaconda已自带了Python。写下本文时Anaconda自带的Python是3.9版本，而Python官网最新版是3.10。若设备同时安装了Python3.10和Anaconda，在命令行终端Anaconda很可能被自动覆盖，导致不可用。此时可以选择卸载Python3.10或者在Anaconda安装目录下启动其自带Python。

### 进入Python

```
import tspex
import pandas as pd
import seaborn as sns
```

本教程将使用[Genotype-Tissue Expression (GTEx) project](https://gtexportal.org/home/index.html)中的数据进行示范。简便起见，对数据进行处理。

```
gtex_link = 'https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz'
expression_data = pd.read_csv(gtex_link, sep='\t', index_col='gene_id', skiprows=2)
expression_data = expression_data.loc[:, ['Bladder', 'Liver', 'Lung', 'Pancreas', 'Stomach']]
expression_data = expression_data.loc[(expression_data > 0).any(axis=1)]
```

```
expression_data.head()
```

| 	 gene_id 		   |Bladder 	|Liver 	   |Lung 	   |Pancreas 	|Stomach    |
|---                 |---        |---        |---	      |---        |---        |
|ENSG00000223972.4 	|0.05878 	|0.06259 	|0.06655 	|0.027255 	|0.063895   |
|ENSG00000227232.4 	|14.24000 	|5.40600 	|13.68000 	|5.553000 	|9.342500   |
|ENSG00000243485.2 	|0.06097 	|0.08316 	|0.06216 	|0.034055 	|0.078575   |
|ENSG00000237613.2 	|0.04113 	|0.03354 	|0.03790 	|0.022915 	|0.043800   |
|ENSG00000268020.2 	|0.00000 	|0.02959 	|0.00000 	|0.000000 	|0.000000   |

### API： `TissueSpecificity`

可使用如下代码查看`TissueSpecificity()`方法详情
```
print(tspex.TissueSpecificity.__doc__)
```

|参数             |数据类型             |值                         |
|---              |---             |---                        |     
|expression_data  |pandas.core.frame.DataFrame|要进行计算的数据，列：基因；行：组织|
|method           |str              |将计算何种指标：'counts', 'tau', 'gini', 'simpson', 'shannon_specificity', 'roku_specificity', 'tsi', 'zscore', 'spm', 'spm_dpm', 'js_specificity', 'js_specificity_dpm'|
|log              |bool             |是否对数据进行取对数处理，默认false|
|transform        |bool             |是否对结果值进行转换至[0,1]范围，默认true（此参数会影响如下指标：`gini`, `simpson`, `shannon_specificity`, `roku_specificity` and `zscore`）|
|threshold        |int or float     |当表达量超过何值时视其表达，默认0(此参数只会影响`counts`指标)|

|属性             |数据类型              |值|
|---              |---              |---|
|expression_data  |pandas.DataFrame |用于计算的数据。若`log=true`,则已被取对数|
|tissue_specificity|pandas.Series / pandas.DataFrame|根据输入得到的计算结果。由于指标类型的不同(General, Individualized)，数据可能为一维或二维|

### 示范

* #### General scoring metric

我们将以ROKU specificity指标为例，示范如何计算基因表达的组织特异性：
```
tso_roku = tspex.TissueSpecificity(expression_data, 'roku_specificity', log=True)
```
以下两个对象分别储存了取对数处理之后的表达数据（如果设置了参数log=True的话）和计算得出的ROKU specificity指标（由于是General scoring metric，所以会是一维的）：
```
tso_roku.expression_data
tso_roku.tissue_specificity
```
还可用以下两个方法绘制直方图和热图：
```
tso_roku.plot_histogram()
```
![histgram_1](/histogram_1.png)

```
tso_roku.plot_heatmap(threshold=0.8, sort_genes=True, use_zscore=True, gene_names=False)
```
![heatmap](/heatmap.png)

对于ROKU specificity，`transform`参数的改变还会影响计算所得结果的形式：
```
tso_roku = tspex.TissueSpecificity(expression_data, 'roku_specificity', log=True, transform=False)
tso_roku.plot_histogram()
```
![histogram_2](/histogram_2.png)


* #### Individualized scoring metric

接下来以SPM指标为例演示Individualized scoring metric的计算：
与上同理，只不过对于Individualized scoring metric，计算结果将会是二维数组。
```
tso_spm = tspex.TissueSpecificity(expression_data, 'spm', log=True)
tso_spm.expression_data
tso_spm.tissue_specificity
```

同样可用一些可视化方法（如seaborn库）对数据进行呈现：
```
sns.violinplot(x="variable", y="value", inner=None, scale="count", color='C0', data=tso_spm.tissue_specificity.melt())
```
![violinplot](/violinplot.png)

```
sns.clustermap(tso_spm.expression_data, figsize=(6,6), z_score=0)
```
![clustermap](/clustermap.png)



