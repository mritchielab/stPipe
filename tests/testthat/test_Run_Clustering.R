test_that("Run_Clustering returns correct structure on small example", {
  
  # 1) create a tiny gene.count: 6 genes × 4 spots
  gene.count <- matrix(sample(0:100, 200 * 100, replace = TRUE), nrow = 200)
  rownames(gene.count) <- paste0("Gene", seq_len(200))
  colnames(gene.count) <- paste0("Spot", seq_len(100))
  
  # 2) matched.data must have spatial_name matching colnames(gene.count)
  matched.data <- data.frame(
     spatial_name = paste0("Spot", seq_len(100)),
     x_coord = runif(100, 0, 10),
     y_coord = runif(100, 0, 10)
  )
  
  # 3) run the function
  out <- Run_Clustering(
    gene.count    = gene.count,
    matched.data  = matched.data,
    num_clusters  = 3
  )
  
  # 4) check it returns a list with named elements
  expect_type(out, "list")
  expected_names <- c("tsne_df", "tsne_plot", "umap_df", "umap_plot",
                      "tsne_interactive", "umap_interactive")
  expect_named(out, expected_names, ignore.order = TRUE)
  
  # 5) check tsne_df and umap_df are data.frames with correct columns
  expect_s3_class(out$tsne_df, "data.frame")
  expect_true(all(c("TSNE1", "TSNE2", "spot", "cluster") %in% colnames(out$tsne_df)))
  expect_equal(nrow(out$tsne_df), length(unique(rownames(unique(t(gene.count))))))
  
  expect_s3_class(out$umap_df, "data.frame")
  expect_true(all(c("UMAP1", "UMAP2", "spot", "cluster") %in% colnames(out$umap_df)))
  expect_equal(nrow(out$umap_df), nrow(unique(t(gene.count))))
  
  # 6) check clusters are factors and have appropriate levels
  expect_s3_class(out$tsne_df$cluster, "factor")
  expect_s3_class(out$umap_df$cluster, "factor")
  
  # 7) check plots are ggplot objects
  expect_s3_class(out$tsne_plot, "gg")
  expect_s3_class(out$umap_plot, "gg")
  
})

