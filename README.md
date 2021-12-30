---
title: "README: scRNAseq_Seurat_Integration_Tutorial v1.0"
author: "Andrew Davis"
date: "2021-12-30"
output:
  github_document:
    toc: TRUE
---

## Overview

This directory contains a tutorial for Seurat's single cell RNA-seq analysis methods, including anchor-based integration. To illustrate these methods, this tutorial includes a comparative analysis of human immune cells (PBMC) in either a resting or interferon-stimulated state 

The following tutorial is designed to give you an overview of the kinds of comparative analyses on complex cell types that are possible using the Seurat integration procedure. Here, we address a few key goals:

  + Create an ‘integrated’ data assay for downstream analysis
  + Identify cell types that are present in both datasets
  + Obtain cell type markers that are conserved in both control and stimulated cells
  + Perform differential expression analysis on cell type clusters to find cell-type specific responses to stimulation


## Lead contact(s):

* [Andrew J. Davis, PhD: 88adavis@gmail.com](mailto:88adavis@gmail.com)

## Input Data

The input data was acquired from ______ on  0/0/0000. Note that the contents of this folder will not be tracked by the remote repository by default (these files tend to be large to store on GitHub). To track the contents of this folder, one must remove the line "/Input/" from the [.gitignore](.gitignore) file.

## Analysis pipeline:

The following pipeline was adapted from the guided tutorials and vignettes listed below:

  + https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html
  + https://github.com/hbctraining/scRNA-seq/tree/master/lessons

Single-Cell RNA-seq Analysis workflow

1. 
2.
...


## Mechanics:

To run this analysis, first create/clean the results output folder by running the codes in [95_Make_Clean.Rmd](). Subsequently, run the code chunks in [99_Run_All.Rmd](R_scripts/99_Run_All.Rmd). This will run the Rmarkdown (.Rmd) files containing the actual code for this analysis in numerical order (i.e., 1.03_Gather-and-Tidy-Data.Rmd, followed by 2.03_Downstream_Analysis.Rmd). Note that one must modify the variable "files_in_r_to_run" in [99_Run_All.Rmd](R_scripts/99_Run_All.Rmd) if one edits or add/deletes filenames of .Rmd scripts associated with this analysis. 

Running [99_Run_All.Rmd](R_scripts/99_Run_All.Rmd) will also render html files of each .Rmd file, which will be saved to the results folder, making useful reports of this analysis. Finally, this README.Rmd files will also be knitted to an html file, as well as a markdown (.md) file, in the working directory of this repository. This markdown file makes for easy viewing on GitHub, and acts as the "home page" for this repo.


## Output:

The resulting output files were saved to the [Results](Results) folder in this repository. Note that the contents of this folder will not be tracked by the remote repository by default (these files tend to be large to store on GitHub). To track the contents of this folder, one must remove the line "/Results/" from the [.gitignore](.gitignore) file.
 

```
## [1] "1.03_Gather-and-Tidy-Data.html"          "2.03_QC_Filter_Data.html"                "3.03_Normalize_Data.html"               
## [4] "4.04_Clustering_finding_Biomarkers.html"
```


## Summary: 

Describe any and all major insights generated from this analysis here....

Presentations and reports shared with other members of our team are stored in [Presentations_Reports](Presentations_Reports). Note that the contents of this folder will not be tracked by the remote repository by default (as these files tend to be large to store on GitHub). To track the contents of this folder, one must remove the line "/Presentations_Reports/" from the [.gitignore](.gitignore) file.


## To do list:

1. 
...

## Template used:
This repository was generated from [10adavis/Rmarkdown_Template_2021](https://github.com/10adavis/Rmarkdown_Template).  

### Session information


```r
sessionInfo()
```

```
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.3 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
##  [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices datasets  utils     methods   base     
## 
## other attached packages:
##  [1] metap_1.7             rmarkdown_2.11        here_1.0.1            cowplot_1.1.1         ggplot2_3.3.5         dplyr_1.0.7          
##  [7] ifnb.SeuratData_3.1.0 patchwork_1.1.1       SeuratData_0.2.1      SeuratObject_4.0.4    Seurat_4.0.6          BiocManager_1.30.16  
## 
## loaded via a namespace (and not attached):
##   [1] sn_2.0.1              plyr_1.8.6            igraph_1.2.10         lazyeval_0.2.2        splines_4.1.2         listenv_0.8.0        
##   [7] scattermore_0.7       TH.data_1.1-0         usethis_2.1.5         digest_0.6.29         htmltools_0.5.2       fansi_0.5.0          
##  [13] magrittr_2.0.1        memoise_2.0.1         tensor_1.5            cluster_2.1.2         ROCR_1.0-11           remotes_2.4.2        
##  [19] globals_0.14.0        matrixStats_0.61.0    sandwich_3.0-1        spatstat.sparse_2.1-0 prettyunits_1.1.1     colorspace_2.0-2     
##  [25] rappdirs_0.3.3        ggrepel_0.9.1         rbibutils_2.2.7       xfun_0.29             callr_3.7.0           crayon_1.4.2         
##  [31] jsonlite_1.7.2        spatstat.data_2.1-2   survival_3.2-13       zoo_1.8-9             glue_1.6.0            polyclip_1.10-0      
##  [37] gtable_0.3.0          leiden_0.3.9          pkgbuild_1.3.1        future.apply_1.8.1    BiocGenerics_0.40.0   abind_1.4-5          
##  [43] scales_1.1.1          mvtnorm_1.1-3         miniUI_0.1.1.1        Rcpp_1.0.7            plotrix_3.8-2         viridisLite_0.4.0    
##  [49] xtable_1.8-4          tmvnsim_1.0-2         reticulate_1.22       spatstat.core_2.3-2   stats4_4.1.2          htmlwidgets_1.5.4    
##  [55] httr_1.4.2            RColorBrewer_1.1-2    TFisher_0.2.0         ellipsis_0.3.2        ica_1.0-2             pkgconfig_2.0.3      
##  [61] farver_2.1.0          sass_0.4.0            uwot_0.1.11           deldir_1.0-6          utf8_1.2.2            tidyselect_1.1.1     
##  [67] labeling_0.4.2        rlang_0.4.12          reshape2_1.4.4        later_1.3.0           munsell_0.5.0         tools_4.1.2          
##  [73] cachem_1.0.6          cli_3.1.0             generics_0.1.1        mathjaxr_1.4-0        devtools_2.4.3        ggridges_0.5.3       
##  [79] evaluate_0.14         stringr_1.4.0         fastmap_1.1.0         yaml_2.2.1            goftest_1.2-3         processx_3.5.2       
##  [85] knitr_1.37            fs_1.5.2              fitdistrplus_1.1-6    purrr_0.3.4           RANN_2.6.1            pbapply_1.5-0        
##  [91] future_1.23.0         nlme_3.1-153          mime_0.12             compiler_4.1.2        rstudioapi_0.13       plotly_4.10.0        
##  [97] curl_4.3.2            png_0.1-7             testthat_3.1.1        spatstat.utils_2.3-0  tibble_3.1.6          bslib_0.3.1          
## [103] stringi_1.7.6         highr_0.9             ps_1.6.0              desc_1.4.0            RSpectra_0.16-0       lattice_0.20-45      
## [109] Matrix_1.4-0          multtest_2.50.0       vctrs_0.3.8           mutoss_0.1-12         pillar_1.6.4          lifecycle_1.0.1      
## [115] Rdpack_2.1.3          jquerylib_0.1.4       spatstat.geom_2.3-1   lmtest_0.9-39         RcppAnnoy_0.0.19      data.table_1.14.2    
## [121] irlba_2.3.5           httpuv_1.6.4          R6_2.5.1              promises_1.2.0.1      renv_0.14.0           KernSmooth_2.23-20   
## [127] gridExtra_2.3         parallelly_1.30.0     sessioninfo_1.2.2     codetools_0.2-18      MASS_7.3-54           pkgload_1.2.4        
## [133] rprojroot_2.0.2       withr_2.4.3           mnormt_2.0.2          sctransform_0.3.2     multcomp_1.4-17       mgcv_1.8-38          
## [139] parallel_4.1.2        grid_4.1.2            rpart_4.1-15          tidyr_1.1.4           Rtsne_0.15            Biobase_2.54.0       
## [145] numDeriv_2016.8-1.1   shiny_1.7.1
```

This document was processed on: 2021-12-30.




