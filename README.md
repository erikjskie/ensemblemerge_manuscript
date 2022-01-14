### The code and datasets for the manuscript:

# EnsembleMerge: Integration of single cell RNA-seq datasets using an ensemble approach

The R package EnsembleMerge can be found at https://github.com/erikjskie/ensemblemerge

## Processed datasets

| | Dataset | Format | Script | 
| --- | --- | --- | --- | 
| Fig 1| [link](https://s3.msi.umn.edu/skiex003/datasets/EnsembleMerge/fig1.csv) | csv | [R](https://github.com/erikjskie/ensemblemerge_manuscript/blob/main/Figure1.R) |
| latent ensemblemerge | [link](https://s3.msi.umn.edu/skiex003/datasets/EnsembleMerge/latent_ensemblemerge.csv) | csv | [R](https://github.com/erikjskie/ensemblemerge_manuscript/blob/main/Figure6.R) |
| scVI ensemblemerge | [link](https://s3.msi.umn.edu/skiex003/datasets/EnsembleMerge/SCVI_ensemblemerge.csv) | csv | [R](https://github.com/erikjskie/ensemblemerge_manuscript/blob/main/Figure7.R) |
| CLR ensemblemerge | [link](https://s3.msi.umn.edu/skiex003/datasets/EnsembleMerge/CLR_normalization.csv) | csv | [R](https://github.com/erikjskie/ensemblemerge_manuscript/blob/main/Figure8.R) |
| fastMNN ensemblemerge | [link](https://s3.msi.umn.edu/skiex003/datasets/EnsembleMerge/Methods_EnsembleMerge.csv) | csv | [R](https://github.com/erikjskie/ensemblemerge_manuscript/blob/main/Figure10.R) |
| Fig 9 | [link](https://s3.msi.umn.edu/skiex003/datasets/EnsembleMerge/fig9.csv) | csv | [R](https://github.com/erikjskie/ensemblemerge_manuscript/blob/main/Figure9.R) |
| fastMNN | [link](https://s3.msi.umn.edu/skiex003/datasets/EnsembleMerge/fastMNN.csv) | csv | [R](https://github.com/erikjskie/ensemblemerge_manuscript/blob/main/fastMNN.R) |
| scVI | [link](https://s3.msi.umn.edu/skiex003/datasets/EnsembleMerge/scVI.csv) | csv | [R](https://github.com/erikjskie/ensemblemerge_manuscript/blob/main/scVI.R) ||

## Data Sources
Datasets sourced from public data
| | Dataset | Format | Script | 
| --- | --- | --- | --- | 
| Villani 2017 | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_1.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| Han 2018 | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_2.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| Hemberg Panc | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_4.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| PBMC | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_5.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| 293t_Jurkat | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_6.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| Hemberg Retina | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_7.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| Saunders 2018 | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_8.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| HCA Blood | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_9.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| Paul 2015 | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_10.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| Zillionis 2019 | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_12.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| Shekhar 2016 | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_13.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| Nestorowa 2016 | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_14.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| Polanski 2019 | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_15.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| Zheng 2017 | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_16.rds) | SummarizedExperiment | [Source](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking) |
| HBDC | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_17.rds) | SummarizedExperiment | [Source](https://github.com/satijalab/seurat-data) |
| Panc8 | [link](https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_18.rds) | SummarizedExperiment | [Source](https://github.com/satijalab/seurat-data) ||

# Notebooks
| Description | link |  
| --- | --- |
| Figure 1 | [[preview](Fig1.ipynb)] [[colab](https://colab.research.google.com/github/erikjskie/ensemblemerge_manuscript/blob/main/Fig1.ipynb)]|
| Supplementary Figure 3 | [[preview](Supplementary_Fig3.ipynb)] [[colab](https://colab.research.google.com/github/erikjskie/ensemblemerge_manuscript/blob/main/Supplementary_Fig3.ipynb)]|
| Supplementary Figure 4 | [[preview](Supplementarty_Fig4.ipynb)] [[colab](https://colab.research.google.com/github/erikjskie/ensemblemerge_manuscript/blob/main/Supplementarty_Fig4.ipynb)]|
| Supplementary Figure 5 | [[preview](supplementary_Fig5.ipynb)] [[colab](https://colab.research.google.com/github/erikjskie/ensemblemerge_manuscript/blob/main/supplementary_Fig5.ipynb)]|
| Supplementary Figure 6 | [[preview](Supplementary_Fig6.ipynb)] [[colab](https://colab.research.google.com/github/erikjskie/ensemblemerge_manuscript/blob/main/Supplementary_Fig6.ipynb)]|
| Supplementary Figure 7 | [[preview](Supplementary_Fig7.ipynb)] [[colab](https://colab.research.google.com/github/erikjskie/ensemblemerge_manuscript/blob/main/Supplementary_Fig7.ipynb)]|
| Supplementary Figure 8 | [[preview](Supplementary_Fig8.ipynb)] [[colab](https://colab.research.google.com/github/erikjskie/ensemblemerge_manuscript/blob/main/Supplementary_Fig8.ipynb)]|
| Supplementary Figure 9 | [[preview](Supplementary_Fig9.ipynb)] [[colab](https://colab.research.google.com/github/erikjskie/ensemblemerge_manuscript/blob/main/Supplementary_Fig9.ipynb)]|
| Supplementary Figure 10 | [[preview](Supplementary_Fig10.ipynb)] [[colab](https://colab.research.google.com/github/erikjskie/ensemblemerge_manuscript/blob/main/Supplementary_Fig10.ipynb)]|
| Supplementary Figure 11 | [[preview](Supplementary_Fig11.ipynb)] [[colab](https://colab.research.google.com/github/erikjskie/ensemblemerge_manuscript/blob/main/Supplementary_Fig11.ipynb)]|

## Instructions for generating the intermediate results

* Intermediate results are generating using R version 4.0.4 with the following packages
  * [EnsembleMerge](https://github.com/erikjskie/ensemblemerge)
  * [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
  * [Seurat](https://cran.r-project.org/web/packages/Seurat/index.html)
  * [Seurat Wrappers](https://github.com/satijalab/seurat-wrappers)
  * [Liger](https://github.com/welch-lab/liger)
  * [Harmony](https://github.com/immunogenomics/harmony)
  * [bbknn](https://github.com/Teichlab/bbknn)
  * [scVI](https://scvi-tools.org/)
  * [Scanorama](https://github.com/brianhie/scanorama)
  * [fastMNN](https://github.com/LTLA/batchelor)
  * [magrittr](https://cran.r-project.org/web/packages/magrittr/index.html)
  * [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)

* The scripts used to generate the data for each figure are named as ```FigureX.R```.
* The scripts will pull all the needed data from AWS storage if they are not already in current directory.
* By default the script runs for ensemblemerge and produces the results in .csv format with the same name as the script in the current directory.
* scripts can be run by:
  *  ```Rscript Figure1.R```
  *  ```Rscript Figure6.R```
  *  ```Rscript Figure7.R```
  *  ```Rscript Figure8.R```
  *  ```Rscript Figure9.R```
  *  ```Rscript Figure10.R```


