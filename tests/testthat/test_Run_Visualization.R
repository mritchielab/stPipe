# Helper: minimal matched.data
make_matched_data <- function(n = 20, seed = 123) {
  set.seed(seed)
  data.frame(
    X_coordinate = runif(n, 0, 100),
    Y_coordinate = runif(n, 0, 100),
    UMI_count    = sample(1:50, n, replace = TRUE),
    spatial_name = paste0("Spot", seq_len(n)),
    stringsAsFactors = FALSE
  )
}

# Test spatial-only visualization
test_that("Run_Visualization returns raw and log UMI ggplot objects when Vis.spatial=TRUE and Vis.read=FALSE", {
  matched <- make_matched_data(30)
  tmp_dir <- tempfile("vis_spatial")
  dir.create(tmp_dir)
  cfg <- list(
    technology_version = "Visium 1.0",
    output_directory   = tmp_dir
  )
  cfg_file <- tempfile(fileext = ".yml")
  yaml::write_yaml(cfg, cfg_file)

  res <- Run_Visualization(
    matched.data = matched,
    config       = cfg_file,
    Vis.spatial  = TRUE,
    Vis.read     = FALSE,
    show.config  = FALSE
  )

  expect_type(res, "list")
  expect_named(res, c("raw_UMI_plot", "log_UMI_plot", "demultiplex_plot", "duplication_plot", "mapping_plot"))

  # raw and log plots are ggplot objects
  expect_s3_class(res$raw_UMI_plot, "ggplot")
  expect_s3_class(res$log_UMI_plot, "ggplot")

  # read-level plots should be NULL
  expect_null(res$demultiplex_plot)
  expect_null(res$duplication_plot)
  expect_null(res$mapping_plot)
})

# Test no visualization when both flags FALSE
test_that("Run_Visualization returns all NULL when Vis.spatial=FALSE and Vis.read=FALSE", {
  matched <- make_matched_data(15)
  tmp_dir <- tempfile("vis_none")
  dir.create(tmp_dir)
  cfg <- list(
    technology_version = "Visium 2.0",
    output_directory   = tmp_dir
  )
  cfg_file <- tempfile(fileext = ".yml")
  yaml::write_yaml(cfg, cfg_file)

  res <- Run_Visualization(
    matched.data = matched,
    config       = cfg_file,
    Vis.spatial  = FALSE,
    Vis.read     = FALSE,
    show.config  = FALSE
  )

  expect_true(all(vapply(res, is.null, logical(1))))
})
