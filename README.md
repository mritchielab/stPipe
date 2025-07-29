# stPipe

<img  src="vignettes/stPipe_logo.png">

The stPipe package provides a comprehensive pipeline for preprocessing sequencing-based spatial transcriptomics data, including 10X Visium, BGI Stereo-seq, Slide-seq, and Curio-seeker.

## Installation

### From Bioconductor (recommended approach)

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("stPipe")
```

### From GitHub (Developmental version)

```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("mritchielab/stPipe")

or

if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("mritchielab/stPipe")
```

The sample data is downsampled from a 10X Visium mouse spleen sample (Sample 709). It is a probe-based, FFPE tissue dataset generated from paired-end FASTQ files. The data can be downloaded from [Zenodo](https://zenodo.org/records/14920583), and a detailed sample description is available in the [SpatialBenchVisium paper by Du *et al.* (2025) Genome Biol 26:77](https://doi.org/10.1186/s13059-025-03543-4).

Sample HTML report for demo data can be found [here](https://github.com/YangXuuu/demo_data_stPipe/blob/main/report.html).
