# define Basilisk environment
stPipe_env <- basilisk::BasiliskEnvironment(
  envname = "stPipe_env",
  pkgname = "stPipe",
  packages = c("python=3.9"),
  channels = c("conda-forge"),
  pip = c(
    "numpy==1.22.1",
    "pandas==1.4.3",
    "pysam==0.22.1",
    "opencv-python==4.5.5.64"
  )
)
