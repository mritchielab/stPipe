make_fake_input <- function(n_spots = 50, n_genes = 20, seed = 42) {
  set.seed(seed)
  gene_ids <- paste0("gene", seq_len(n_genes))
  mat <- matrix(rpois(n_spots * n_genes, lambda = 10), nrow = n_genes, ncol = n_spots)
  colnames(mat) <- paste0("SP", seq_len(n_spots))
  gene_matrix <- data.frame(gene_id = gene_ids, mat, stringsAsFactors = FALSE)

  umi_counts <- colSums(mat)
  md <- data.frame(
    X_coordinate     = runif(n_spots, 0, 100),
    Y_coordinate     = runif(n_spots, 0, 100),
    barcode_sequence = paste0("BC", seq_len(n_spots)),
    spatial_name     = colnames(mat),
    UMI_count        = umi_counts,
    stringsAsFactors = FALSE
  )
  list(gene_matrix = gene_matrix, matched_data = md)
}

test_that("Run_QC throws ggplot aesthetic-length error on slope_max", {
  inp     <- make_fake_input(n_spots = 100, n_genes = 10, seed = 123)
  tmp_dir <- tempfile("qc_slope")
  dir.create(tmp_dir)
  cfg <- list(
    output_directory = tmp_dir,
    qc_filter        = "slope_max",
    qc_per           = "0.2_0.8"
  )
  cfg_file <- tempfile(fileext = ".yml")
  yaml::write_yaml(cfg, cfg_file)

  expect_error(
    Run_QC(
      config       = cfg_file,
      matched.data = inp$matched_data,
      gene.matrix  = inp$gene_matrix,
      show.config  = FALSE
    ),
    "Aesthetics must be either length 1 or the same as the data"
  )
})

test_that("Run_QC throws ggplot aesthetic-length error on EmptyDropletUtils", {
  skip_if_not_installed("DropletUtils")

  inp     <- make_fake_input(n_spots = 80, n_genes = 15, seed = 456)
  tmp_dir <- tempfile("qc_empty")
  dir.create(tmp_dir)
  cfg <- list(
    output_directory = tmp_dir,
    qc_filter        = "EmptyDropletUtils",
    qc_per           = "0.1_0.9"
  )
  cfg_file <- tempfile(fileext = ".yml")
  yaml::write_yaml(cfg, cfg_file)

  expect_error(
    Run_QC(
      config       = cfg_file,
      matched.data = inp$matched_data,
      gene.matrix  = inp$gene_matrix,
      show.config  = FALSE
    ),
    "Aesthetics must be either length 1 or the same as the data"
  )
})


