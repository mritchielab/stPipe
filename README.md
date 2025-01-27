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

The stPipe vignette can be found using:
```
# stPipe must be installed
browseVignettes("stPipe")
```

Sample 10X Visium data used in the vignette can be found [here](https://github.com/YangXuuu/demo_data_stPipe).
