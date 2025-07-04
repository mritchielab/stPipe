# Helper: write a minimal config yaml
write_config <- function(cfg_list) {
  cfg_file <- tempfile(fileext = ".yml")
  yaml::write_yaml(cfg_list, cfg_file)
  cfg_file
}

# Test that mismatched data and output directories throw an error
test_that("Run_ST stops when data and output dirs count mismatch", {
  cfg_list <- list(
    data_directory    = "dir1,dir2",
    output_directory  = "out1",
    technology_version= "Visium_probe_v1",
    species           = "mouse",
    analysis_source   = "",
    index_fa          = "",
    index_gff3        = "",
    scpipe_nthreads   = 1,
    max_reads         = 10,
    min_count         = 1,
    number_of_locations = 5,
    bs1 = -1, bl1 = 0, bs2 = 0, bl2 = 16, us = 16, ul = 12
  )
  cfg_file <- write_config(cfg_list)
  expect_error(
    Run_ST(config = cfg_file, show.config = FALSE),
    "The number of data directories and output directories must match!"
  )
})

# Visium inconsistency: more/less R1 or R2 files
test_that("Run_ST Visium stops when multiple or missing FASTQ files", {
  tmp_data <- tempfile("data")
  tmp_out  <- tempfile("out")
  dir.create(tmp_data); dir.create(tmp_out)
  # create two R1 and one R2
  file.create(file.path(tmp_data, "sample_R1_1.fastq"))
  file.create(file.path(tmp_data, "sample_R1_2.fastq"))
  file.create(file.path(tmp_data, "sample_R2.fastq"))
  cfg <- list(
    data_directory    = tmp_data,
    output_directory  = tmp_out,
    technology_version= "Visium_v1",
    species           = "mouse",
    analysis_source   = "",
    index_fa          = "",
    index_gff3        = "",
    scpipe_nthreads   = 1,
    max_reads         = 10,
    min_count         = 1,
    number_of_locations = 5,
    bs1 = -1, bl1 = 0, bs2 = 0, bl2 = 16, us = 16, ul = 12
  )
  cfg_file <- write_config(cfg)
  expect_error(
    Run_ST(config = cfg_file, show.config = FALSE),
    "Expected exactly one R1 and one R2 FASTQ file, but found inconsistencies"
  )
})

# Unsupported combination of technology_version and species in probe-based Visium
test_that("Run_ST stops on unsupported probe-based Visium combination", {
  tmp_data <- tempfile("data")
  tmp_out  <- tempfile("out")
  dir.create(tmp_data); dir.create(tmp_out)
  # create one R1 and one R2 to pass earlier checks
  file.create(file.path(tmp_data, "sample_R1.fastq"))
  file.create(file.path(tmp_data, "sample_R2.fastq"))
  cfg <- list(
    data_directory    = tmp_data,
    output_directory  = tmp_out,
    technology_version= "Visium_probe_v2", # v2 & mouse unsupported
    species           = "mouse",
    analysis_source   = "",
    index_fa          = "",
    index_gff3        = "",
    scpipe_nthreads   = 1,
    max_reads         = 10,
    min_count         = 1,
    number_of_locations = 5,
    bs1 = -1, bl1 = 0, bs2 = 0, bl2 = 16, us = 16, ul = 12
  )
  cfg_file <- write_config(cfg)
  expect_error(
    Run_ST(config = cfg_file, show.config = FALSE),
    "Unsupported combination of technology_version and species for probe-based Visium"
  )
})


