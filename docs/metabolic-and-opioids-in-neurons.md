---
title: "Expression analysis of PVN neurons with focus on Crh"
author: "Evgenii O. Tretiakov"
date: "2024-11-26"
format:
  html:
    toc: true
    df-print: paged
    code-fold: true
    fig-width: 9
    fig-height: 12
    fig-format: retina
    fig-responsive: true
    fig-dpi: 300
  pdf:
    colorlinks: true
    fontsize: 12pt
execute:
  keep-md: true
  echo: true
  error: false
  message: false
  warning: false
  debug: false
knitr:
  opts_chunk:
    autodep: true
    fig.align: center
    fig.retina: 2
    fig.width: 14
    fig.height: 12
---





## Setup parameters


::: {.cell layout-align="center"}

```{.r .cell-code}
# Load tidyverse infrastructure packages
suppressPackageStartupMessages({
  library(future)
  library(here)
  library(tidyverse)
  library(magrittr)
  library(stringr)
  library(skimr)
  library(RColorBrewer)
  library(viridis)
})


# Load packages for scRNA-seq analysis and visualisation
suppressPackageStartupMessages({
  library(UpSetR)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(ggstatsplot)
  library(anndata)
  library(sceasy)
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(scCustomize)
})

sc <- import("scanpy", convert = FALSE)
```
:::


### Set paths


::: {.cell layout-align="center"}

```{.r .cell-code}
src_dir <- here("code")
data_dir <- here("data")
output_dir <- here("output")
plots_dir <- here(output_dir, "figures/")
tables_dir <- here(output_dir, "tables/")
```
:::


### Load helper functions and gene-sets


::: {.cell layout-align="center"}

```{.r .cell-code}
source(here(src_dir, "genes.R"))
source(here(src_dir, "functions.R"))
```
:::


### Set fixed variables


::: {.cell layout-align="center"}

```{.r .cell-code}
# set seed
reseed <- 42
set.seed(seed = reseed)

# Parameters for parallel execution
n_cores <- 8
plan("multisession", workers = n_cores)
options(
  future.globals.maxSize = 100000 * 1024^2,
  future.rng.onMisuse = "ignore"
)
plan()
```

::: {.cell-output .cell-output-stdout}
```
multisession:
- args: function (..., workers = 8, envir = parent.frame())
- tweaked: TRUE
- call: plan("multisession", workers = n_cores)
```
:::

```{.r .cell-code}
# ggplot2 theme
theme_set(ggmin::theme_powerpoint())
```
:::


## Load selected astrocytes data from Lopez JP et al (2021)


::: {.cell layout-align="center"}

```{.r .cell-code}
anndata <- sc$read(here(
  data_dir,
  "class_cello/PRJNA679294-whole_dataset-0.001-cello_annotation.h5ad"
))
```
:::


### Convert adata object to R AnnDataR6 object.

::: {.cell layout-align="center"}

```{.r .cell-code}
adata <- py_to_r(anndata)
class(adata)
```

::: {.cell-output .cell-output-stdout}
```
[1] "AnnDataR6" "R6"       
```
:::

```{.r .cell-code}
class(adata$X)
```

::: {.cell-output .cell-output-stdout}
```
[1] "dgCMatrix"
attr(,"package")
[1] "Matrix"
```
:::

```{.r .cell-code}
adata
```

::: {.cell-output .cell-output-stdout}
```
AnnData object with n_obs × n_vars = 9572 × 22835
    obs: 'nCount_RAW', 'nFeature_RAW', 'nCount_RNA', 'nFeature_RNA', 'orig.ident', 'nFeature_Diff', 'nCount_Diff', 'percent_mito', 'percent_ribo', 'percent_mito_ribo', 'percent_hb', 'log10GenesPerUMI', 'cell_name', 'barcode', 'latent_RT_efficiency', 'latent_cell_probability', 'latent_scale', 'doublet_score', 'predicted_doublets', 'QC', 'var_regex', 'RNA_snn_res.0.5', 'RNA_snn_res.0.741576329532297', 'RNA_snn_res.1.24913691074005', 'RNA_snn_res.3.000001', 'seurat_clusters', 'k_tree', 'comb_clstr1', 'S.Score', 'G2M.Score', 'Phase', 'nCount_SCT', 'nFeature_SCT', 'SCT_snn_res.1', 'SCT_snn_res.1.1014548330962', 'SCT_snn_res.1.22223389211258', 'SCT_snn_res.1.36843938820807', 'SCT_snn_res.1.54904506877152', 'SCT_snn_res.1.7778105049838', 'SCT_snn_res.2.07696233957953', 'SCT_snn_res.2.48489124123347', 'SCT_snn_res.3.07411084005006', 'SCT_snn_res.4.00000100000002', 'bioproject', 'project', 'model', 'tech', 'region', 'sex', 'stage', 'libname', 'expbtch', 'condit', 'ora_celltype'
    var: 'vst.mean', 'vst.variance', 'vst.variance.expected', 'vst.variance.standardized', 'vst.variable'
    uns: 'k_tree_colors', 'name', 'ora_celltype_colors'
    obsm: 'X_pacmap', 'X_pca', 'X_umap', 'ora_estimate', 'ora_pvals'
```
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
expr_mtx <- t(as.matrix(adata$raw$X))
colnames(expr_mtx) <- rownames(adata$X)
rownames(expr_mtx) <- adata$var_names
srt <- CreateSeuratObject(
  expr_mtx,
  assay = "RNA",
  project = "individual_hypothalamic_nuclei_astrocytes_evaluation_dataset",
  meta.data = as.data.frame(adata$obs)
)

Idents(srt) <- "ora_celltype"
table(Idents(srt))
```

::: {.cell-output .cell-output-stdout}
```

                      Astrocytes                Endothelial cells 
                            2838                              493 
                 Ependymal cells                      Macrophages 
                             472                               30 
                 Mesangial cells                        Microglia 
                              52                              463 
                         Neurons Oligodendrocyte progenitor cells 
                            3309                              383 
                Oligodendrocytes                        Pericytes 
                            1275                              257 
```
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt <- subset(srt, idents = c("Neurons"))

Idents(srt) <- "libname"

print(srt)
```

::: {.cell-output .cell-output-stdout}
```
An object of class Seurat 
22835 features across 3309 samples within 1 assay 
Active assay: RNA (22835 features, 0 variable features)
 1 layer present: counts
```
:::

```{.r .cell-code}
rm(adata, anndata, expr_mtx)
invisible(gc())
table(Idents(srt))
```

::: {.cell-output .cell-output-stdout}
```

control    csds 
   1677    1632 
```
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
opioid_system_genes <- c(
    # Classical opioid receptors
    "Oprd1",  # Delta opioid receptor
    "Oprk1",  # Kappa opioid receptor 
    "Oprl1",  # Nociceptin/orphanin FQ receptor
    "Oprm1",  # Mu opioid receptor
    
    # Processing enzymes
    "Pcsk1",  # Proprotein convertase 1
    "Pcsk2",   # Proprotein convertase 2
    
    # Endogenous opioid precursors
    "Pdyn",   # Prodynorphin
    "Penk",   # Proenkephalin
    #"Pomc",   # Proopiomelanocortin
    "Pnoc"   # Prepronociceptin
)
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
metabolic_signaling_genes <- c(
    # Receptor tyrosine kinases and ligands
    "Alk",      # Anaplastic lymphoma kinase - neural development, metabolism
    "Fam150a",  # ALK ligand 1/Augmentor-α - ALK receptor activator
    "Fam150b",  # ALK ligand 2/Augmentor-β - ALK receptor activator
    
    # Melanocortin system
    "Mc3r",     # Melanocortin 3 receptor - energy homeostasis, inflammation
    "Mc4r",     # Melanocortin 4 receptor - appetite control, energy balance
    
    # Metabolic hormone receptors
    "Lepr",     # Leptin receptor - energy balance, satiety
    "Insr",     # Insulin receptor - glucose homeostasis
    #"Igf1r",    # Insulin-like growth factor 1 receptor - growth, development
    
    # Signaling adaptors/regulators
    "Lmo4",     # LIM domain only 4 - transcriptional regulation, metabolism
    "Irs1",     # Insulin receptor substrate 1 - insulin signaling
    "Irs4"      # Insulin receptor substrate 4 - insulin/leptin signaling
)
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
all.genes <- rownames(srt)
gene.scale <- c(
  cnbn,
  opioid_system_genes,
  metabolic_signaling_genes,
  np,
  npr,
  nmr,
  neurotrans
  ) |> 
  unique() %>%
  .[. %in% all.genes]
srt <- NormalizeData(srt)
srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 3000)
srt <- ScaleData(srt, features = c(VariableFeatures(srt), gene.scale))
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.bck <- srt
srt <- subset(srt, subset = Oxt >= 1 | Trh > 0 | Crh > 0 | Sst >= 1 | Avp >= 1)
Idents(srt) <- "ora_celltype"
table(Idents(srt))
```

::: {.cell-output .cell-output-stdout}
```

Neurons 
   1034 
```
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt <- RunPCA(srt, npcs = 20, verbose = FALSE)


invisible(gc())
set.seed(reseed)
if (!file.exists(here(data_dir, glue::glue("{DOCNAME}-init/{DOCNAME}-init-umap-search-ref.Rds")))) {
  
  source(here(src_dir, "scDEED.R"))
  library(furrr)

  permuted.srt <- Permuted(srt, K = 20)

  invisible(gc())
  set.seed(reseed)

  umap_example <- scDEED(
    input_data = srt,
    K = 20,
    n_neighbors = seq(from = 15, to = 35, by = 10),
    min.dist = c(0.05, 0.1, 0.25, 0.5, 0.8),
    reduction.method = "umap",
    rerun = FALSE,
    permuted = permuted.srt,
    default_assay = "RNA"
  )

  dir.create(here(data_dir, sprintf("%s-init", DOCNAME)))
  readr::write_rds(
    x = umap_example,
    file = here(data_dir, glue::glue("{DOCNAME}-init/{DOCNAME}-init-umap-search-ref.Rds"))
  )
} else {
  umap_example <-
    read_rds(here(data_dir, glue::glue("{DOCNAME}-init/{DOCNAME}-init-umap-search-ref.Rds")))
}
```

::: {.cell-output .cell-output-stdout}
```

  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |=======                                                               |  10%
  |                                                                            
  |============================                                          |  40%
  |                                                                            
  |============================                                          |  41%
  |                                                                            
  |=============================                                         |  41%
  |                                                                            
  |=============================                                         |  42%
  |                                                                            
  |==============================                                        |  42%
  |                                                                            
  |==============================                                        |  43%
  |                                                                            
  |==============================                                        |  44%
  |                                                                            
  |===============================                                       |  44%
  |                                                                            
  |===============================                                       |  45%
  |                                                                            
  |================================                                      |  45%
  |                                                                            
  |================================                                      |  46%
  |                                                                            
  |=================================                                     |  46%
  |                                                                            
  |=================================                                     |  47%
  |                                                                            
  |=================================                                     |  48%
  |                                                                            
  |==================================                                    |  48%
  |                                                                            
  |==================================                                    |  49%
  |                                                                            
  |===================================                                   |  49%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |===================================                                   |  51%
  |                                                                            
  |====================================                                  |  51%
  |                                                                            
  |====================================                                  |  52%
  |                                                                            
  |=====================================                                 |  52%
  |                                                                            
  |=====================================                                 |  53%
  |                                                                            
  |=====================================                                 |  54%
  |                                                                            
  |======================================                                |  54%
  |                                                                            
  |======================================                                |  55%
  |                                                                            
  |=======================================                               |  55%
  |                                                                            
  |=======================================                               |  56%
  |                                                                            
  |========================================                              |  56%
  |                                                                            
  |========================================                              |  57%
  |                                                                            
  |========================================                              |  58%
  |                                                                            
  |=========================================                             |  58%
  |                                                                            
  |=========================================                             |  59%
  |                                                                            
  |==========================================                            |  59%
  |                                                                            
  |==========================================                            |  60%
  |                                                                            
  |==========================================                            |  61%
  |                                                                            
  |===========================================                           |  61%
  |                                                                            
  |===========================================                           |  62%
  |                                                                            
  |============================================                          |  62%
  |                                                                            
  |============================================                          |  63%
  |                                                                            
  |============================================                          |  64%
  |                                                                            
  |=============================================                         |  64%
  |                                                                            
  |=============================================                         |  65%
  |                                                                            
  |==============================================                        |  65%
  |                                                                            
  |==============================================                        |  66%
  |                                                                            
  |===============================================                       |  66%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |===============================================                       |  68%
  |                                                                            
  |================================================                      |  68%
  |                                                                            
  |================================================                      |  69%
  |                                                                            
  |=================================================                     |  69%
  |                                                                            
  |=================================================                     |  70%
  |                                                                            
  |=================================================                     |  71%
  |                                                                            
  |==================================================                    |  71%
  |                                                                            
  |==================================================                    |  72%
  |                                                                            
  |===================================================                   |  72%
  |                                                                            
  |===================================================                   |  73%
  |                                                                            
  |===================================================                   |  74%
  |                                                                            
  |====================================================                  |  74%
  |                                                                            
  |====================================================                  |  75%
  |                                                                            
  |=====================================================                 |  75%
  |                                                                            
  |=====================================================                 |  76%
  |                                                                            
  |======================================================                |  76%
  |                                                                            
  |======================================================                |  77%
  |                                                                            
  |======================================================                |  78%
  |                                                                            
  |=======================================================               |  78%
  |                                                                            
  |=======================================================               |  79%
  |                                                                            
  |========================================================              |  79%
  |                                                                            
  |========================================================              |  80%
```
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt <-
  srt |>
  FindNeighbors(
    dims = 1:20,
    k.param = umap_example$num_dubious |> 
      dplyr::slice_min(order_by = number_dubious_cells, n = 1) |> 
      pull(n_neighbors),
    annoy.metric = "cosine",
    n.trees = 100,
    verbose = FALSE
  ) |>
  RunUMAP(
    dims = 1:20,
    reduction.name = "umap",
    reduction.key = "UMAP_",
    return.model = FALSE,
    umap.method = "uwot",
    n.epochs = 1000L,
    n.neighbors = umap_example$num_dubious |> 
      dplyr::slice_min(order_by = number_dubious_cells, n = 1) |> 
      pull(n_neighbors),
    min.dist = umap_example$num_dubious |> 
      dplyr::slice_min(order_by = number_dubious_cells, n = 1) |> 
      pull(min.dist),
    seed.use = reseed,
    verbose = FALSE
  )
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
Idents(srt) <- "libname"
DimPlot(srt)
```

::: {.cell-output-display}
![](metabolic-and-opioids-in-neurons_files/figure-html/plot-lopez2021-pvn-condit-astro-1.png){fig-align='center' width=1800}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
FeaturePlot_scCustom(
  srt,
  reduction = "umap",
  features = c(
    metabolic_signaling_genes, 
    opioid_system_genes,
    "Slc17a6", "Gad1", "Gad2", "Crh", "Trh", "Oxt", "Sst"
    ),
  layer = "data",
  label = F,
  num_columns = 4
) #* NoLegend()
```

::: {.cell-output-display}
![](metabolic-and-opioids-in-neurons_files/figure-html/plot-lopez2021-pvn-feature-cb-1.png){fig-align='center' width=5400}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sbs_mtx <-
  srt@assays$RNA@layers$data %>%
  as.data.frame() %>%
  t()

rownames(sbs_mtx) <- colnames(srt)
colnames(sbs_mtx) <- rownames(srt)

# Filter features
filt_low_genes <-
  colSums(sbs_mtx) %>%
  .[. > quantile(., 0.4)] %>%
  names()
sbs_mtx %<>% .[, filt_low_genes]

min_filt_vector2 <-
  sbs_mtx %>%
  as_tibble() %>%
  select(all_of(filt_low_genes)) %>%
  summarise(across(.fns = ~ quantile(.x, .005))) %>%
  as.list() %>%
  map(as.double) %>%
  simplify() %>%
  .[filt_low_genes]

# Prepare table of intersection sets analysis
content_sbs_mtx <-
  (sbs_mtx > min_filt_vector2) %>%
  as_tibble() %>%
  mutate_all(as.numeric) %>%
  bind_cols(
    srt@meta.data |> select(condit)
  )
```
:::

::: {.cell layout-align="center" fig.asp='1.214'}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
      filter(condit == 0) |>
      select(
        c(metabolic_signaling_genes, "Crh", "Trh", "Oxt") %>% .[. %in% colnames(content_sbs_mtx)]
      )
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 30,
  sets = c(metabolic_signaling_genes, "Crh", "Trh", "Oxt") %>%
    .[. %in% colnames(content_sbs_mtx)],
  empty.intersections = NULL
)
```

::: {.cell-output-display}
![](metabolic-and-opioids-in-neurons_files/figure-html/upset-group-e-cb-pvn-Adult-f2-astro-norm-1.png){fig-align='center' width=4200}
:::

```{.r .cell-code}
skim(as.data.frame(
    content_sbs_mtx |>
      filter(condit == 0) |>
      select(
        c(metabolic_signaling_genes, "Crh", "Trh", "Oxt") %>% .[. %in% colnames(content_sbs_mtx)]
      )
  ))
```

::: {.cell-output-display}
Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |559                |
|Number of columns        |9                  |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |9                  |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Alk           |         0|             1| 0.35| 0.48|  0|   0|   0|   1|    1|▇▁▁▁▅ |
|Lepr          |         0|             1| 0.04| 0.19|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Insr          |         0|             1| 0.40| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▆ |
|Lmo4          |         0|             1| 0.27| 0.44|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Irs1          |         0|             1| 0.10| 0.30|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Irs4          |         0|             1| 0.23| 0.42|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Crh           |         0|             1| 0.12| 0.32|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Trh           |         0|             1| 0.33| 0.47|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Oxt           |         0|             1| 0.14| 0.34|  0|   0|   0|   0|    1|▇▁▁▁▁ |
:::
:::

::: {.cell layout-align="center" fig.asp='1.214'}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
      filter(condit == 0) |>
      select(
        c(opioid_system_genes, "Crh", "Trh", "Oxt") %>% .[. %in% colnames(content_sbs_mtx)]
      )
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 30,
  sets = c(opioid_system_genes, "Crh", "Trh", "Oxt") %>%
    .[. %in% colnames(content_sbs_mtx)],
  empty.intersections = NULL
)
```

::: {.cell-output-display}
![](metabolic-and-opioids-in-neurons_files/figure-html/upset-group-e-cb-pvn-Adult-f3-astro-norm-1.png){fig-align='center' width=4200}
:::

```{.r .cell-code}
skim(as.data.frame(
    content_sbs_mtx |>
      filter(condit == 0) |>
      select(
        c(opioid_system_genes, "Crh", "Trh", "Oxt") %>% .[. %in% colnames(content_sbs_mtx)]
      )
  ))
```

::: {.cell-output-display}
Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |559                |
|Number of columns        |12                 |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |12                 |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Oprd1         |         0|             1| 0.06| 0.23|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprk1         |         0|             1| 0.16| 0.36|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Oprl1         |         0|             1| 0.42| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▆ |
|Oprm1         |         0|             1| 0.27| 0.44|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Pcsk1         |         0|             1| 0.33| 0.47|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Pcsk2         |         0|             1| 0.68| 0.47|  0|   0|   1|   1|    1|▃▁▁▁▇ |
|Pdyn          |         0|             1| 0.19| 0.39|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Penk          |         0|             1| 0.18| 0.39|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Pnoc          |         0|             1| 0.15| 0.36|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Crh           |         0|             1| 0.12| 0.32|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Trh           |         0|             1| 0.33| 0.47|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Oxt           |         0|             1| 0.14| 0.34|  0|   0|   0|   0|    1|▇▁▁▁▁ |
:::
:::

::: {.cell layout-align="center" fig.asp='1.214'}

```{.r .cell-code}
upset(
  as.data.frame(
    content_sbs_mtx |>
      filter(condit == 0) |>
      select(
        c(metabolic_signaling_genes, opioid_system_genes, "Crh", "Trh", "Oxt") %>% .[. %in% colnames(content_sbs_mtx)]
      )
  ),
  order.by = "freq",
  cutoff = 3,
  sets.x.label = "Number of cells",
  number.angles = 0,
  point.size = 3.5, line.size = 2,
  text.scale = c(2, 1.6, 2, 1.3, 2, 1.1),
  nsets = 30,
  nintersects = 30,
  sets = c(metabolic_signaling_genes, opioid_system_genes, "Crh", "Trh", "Oxt") %>%
    .[. %in% colnames(content_sbs_mtx)],
  empty.intersections = NULL
)
```

::: {.cell-output-display}
![](metabolic-and-opioids-in-neurons_files/figure-html/upset-group-e-cb-pvn-Adult-f4-astro-norm-1.png){fig-align='center' width=4200}
:::

```{.r .cell-code}
skim(as.data.frame(
    content_sbs_mtx |>
      filter(condit == 0) |>
      select(
        c(metabolic_signaling_genes, opioid_system_genes, "Crh", "Trh", "Oxt") %>% .[. %in% colnames(content_sbs_mtx)]
      )
  ))
```

::: {.cell-output-display}
Table: Data summary

|                         |                   |
|:------------------------|:------------------|
|Name                     |as.data.frame(...) |
|Number of rows           |559                |
|Number of columns        |18                 |
|_______________________  |                   |
|Column type frequency:   |                   |
|numeric                  |18                 |
|________________________ |                   |
|Group variables          |None               |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Alk           |         0|             1| 0.35| 0.48|  0|   0|   0|   1|    1|▇▁▁▁▅ |
|Lepr          |         0|             1| 0.04| 0.19|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Insr          |         0|             1| 0.40| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▆ |
|Lmo4          |         0|             1| 0.27| 0.44|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Irs1          |         0|             1| 0.10| 0.30|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Irs4          |         0|             1| 0.23| 0.42|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Oprd1         |         0|             1| 0.06| 0.23|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Oprk1         |         0|             1| 0.16| 0.36|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Oprl1         |         0|             1| 0.42| 0.49|  0|   0|   0|   1|    1|▇▁▁▁▆ |
|Oprm1         |         0|             1| 0.27| 0.44|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Pcsk1         |         0|             1| 0.33| 0.47|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Pcsk2         |         0|             1| 0.68| 0.47|  0|   0|   1|   1|    1|▃▁▁▁▇ |
|Pdyn          |         0|             1| 0.19| 0.39|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Penk          |         0|             1| 0.18| 0.39|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Pnoc          |         0|             1| 0.15| 0.36|  0|   0|   0|   0|    1|▇▁▁▁▂ |
|Crh           |         0|             1| 0.12| 0.32|  0|   0|   0|   0|    1|▇▁▁▁▁ |
|Trh           |         0|             1| 0.33| 0.47|  0|   0|   0|   1|    1|▇▁▁▁▃ |
|Oxt           |         0|             1| 0.14| 0.34|  0|   0|   0|   0|    1|▇▁▁▁▁ |
:::
:::
