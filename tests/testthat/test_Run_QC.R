# tests/test_Run_QC.R
# ──────────────────────────────────────────────────────────────────────────────
# This file contains unit tests for the Run_QC function in the stPipe package.
# We verify both slope_max and EmptyDropletUtils filtering modes.

make_fake_input <- function(n_spots = 50, n_genes = 20, seed = 42) {
  set.seed(seed)
  gene_ids <- paste0("gene", seq_len(n_genes))
  mat <- matrix(rpois(n_spots * n_genes, lambda = 10),
                nrow = n_genes, ncol = n_spots)
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

# Test that slope_max mode filters out low-UMI spots
testthat::test_that("slope_max mode removes low-UMI spots", {
  inp <- make_fake_input(n_spots = 100, n_genes = 10, seed = 123)
  tmp_dir <- tempfile("qc_slope"); dir.create(tmp_dir)
  cfg <- list(
    output_directory = tmp_dir,
    qc_filter        = "slope_max",
    qc_per           = "0.2_0.8"
  )
  cfg_file <- tempfile(fileext = ".yml")
  yaml::write_yaml(cfg, cfg_file)

  res <- stPipe::Run_QC(
    config       = cfg_file,
    matched.data = inp$matched_data,
    gene.matrix  = inp$gene_matrix,
    show.config  = FALSE
  )

  # Result must be a list with correct names
  testthat::expect_type(res, "list")
  testthat::expect_named(res, c("filtered_gene_count", "filtered_spatial_coords"))

  # Verify some spots were removed in both outputs
  testthat::expect_lt(ncol(res$filtered_gene_count), ncol(inp$gene_matrix))
  testthat::expect_lt(nrow(res$filtered_spatial_coords), nrow(inp$matched_data))
})

# Test that EmptyDropletUtils mode calls emptyDrops and filters correctly (mocked)
testthat::test_that("EmptyDropletUtils mode applies emptyDrops filter correctly", {
  # Skip if DropletUtils is not installed
  testthat::skip_if_not_installed("DropletUtils")

  inp <- make_fake_input(n_spots = 80, n_genes = 15, seed = 456)
  tmp_dir <- tempfile("qc_empty"); dir.create(tmp_dir)
  cfg <- list(
    output_directory = tmp_dir,
    qc_filter        = "EmptyDropletUtils",
    qc_per           = "0.1_0.9"
  )
  cfg_file <- tempfile(fileext = ".yml")
  yaml::write_yaml(cfg, cfg_file)

  # Create dummy output for emptyDrops: random PValue and FDR columns
  dummy_outs <- data.frame(
    PValue = runif(ncol(inp$gene_matrix) - 1, 0, 0.1),
    FDR    = runif(ncol(inp$gene_matrix) - 1, 0, 0.01)
  )

  # Mock DropletUtils::emptyDrops to return dummy_outs
  testthat::with_mocked_bindings(
    `DropletUtils::emptyDrops` = function(mat, retain) {
      # mat is expected as a data.frame of counts
      testthat::expect_s3_class(mat, "data.frame")
      dummy_outs
    },
    {
      res <- stPipe::Run_QC(
        config       = cfg_file,
        matched.data = inp$matched_data,
        gene.matrix  = inp$gene_matrix,
        show.config  = FALSE
      )

      # Determine which indices should be kept based on dummy_outs
      keep_idx <- which(dummy_outs$FDR <= 0.01 & dummy_outs$PValue <= 0.05)

      # Check that filtered gene count columns match keep_idx length
      testthat::expect_equal(ncol(res$filtered_gene_count), length(keep_idx))

      # Ensure that spatial_names in coords match the filtered columns
      testthat::expect_true(all(res$filtered_spatial_coords$spatial_name %in% colnames(res$filtered_gene_count)))
      testthat::expect_lte(nrow(res$filtered_spatial_coords), ncol(res$filtered_gene_count))
    }
  )
})
