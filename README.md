# stPipe

<img  src="vignettes/stPipe_logo.png">

The stPipe package provides a comprehensive pipeline for preprocessing sequencing-based spatial transcriptomics data, including 10X Visium, BGI Stereo-seq, Slide-seq, and Curio-seeker.

The latest development version can be installed with:
```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("mritchielab/stPipe")

or

if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("mritchielab/stPipe")
```

The sample data is downsampled from a 10X Visium mouse spleen sample (Sample 709). It is a probe-based, FFPE tissue dataset generated from paired-end FASTQ files. The data can be downloaded from [Zenodo](https://zenodo.org/records/14920583), and a detailed sample description is available in the [SpatialBench paper](https://www.biorxiv.org/content/10.1101/2024.03.13.584910v1.abstract).

Sample HTML report for demo data can be found [here](https://github.com/YangXuuu/demo_data_stPipe/blob/main/report.html).
