# Helper: generate dummy gene.matrix and matched.data
make_input <- function(n_genes = 50, n_spots = 30, seed = 42) {
  set.seed(seed)
  gene_ids <- paste0("Gene", seq_len(n_genes))
  spot_ids <- paste0("Spot", seq_len(n_spots))
  mat <- matrix(sample(0:100, n_genes * n_spots, replace = TRUE),
                nrow = n_genes, ncol = n_spots,
                dimnames = list(gene_ids, spot_ids))
  gene_matrix <- as.data.frame(mat)
  gene_matrix$gene_id <- rownames(gene_matrix)
  gene_matrix <- gene_matrix[, c("gene_id", spot_ids)]

  md <- data.frame(
    spatial_name     = spot_ids,
    barcode_sequence = paste0("BC", seq_len(n_spots)),
    X_coordinate     = runif(n_spots, 0, 100),
    Y_coordinate     = runif(n_spots, 0, 100),
    UMI_count        = sample(1:500, n_spots, replace = TRUE),
    stringsAsFactors = FALSE
  )
  list(gene_matrix = gene_matrix, matched_data = md)
}

# error when matched.data has unmatched spatial_name
test_that("Run_Create_Obj stops when spatial_name mismatch", {
  inp <- make_input()
  inp$matched_data$spatial_name[1] <- "WrongSpot"
  expect_error(
    Run_Create_Obj(
      gene.matrix   = inp$gene_matrix,
      matched.data  = inp$matched_data,
      obj.type      = "Seurat",
      tech          = "Visium"
    ),
    "Spatial locations in matched.data do not match column names of gene.matrix!"
  )
})

