test_that("pseudo-output has correct structure", {
  
  # 1) matched_data pseudo
  matched_data <- data.frame(
    barcode_sequence = c("B1", "B2", "B3"),
    spatial_name     = c("SPATIAL_1", "SPATIAL_2", "SPATIAL_3"),
    X_coordinate     = c(10, 20, 30),
    Y_coordinate     = c(40, 50, 60),
    UMI_count        = c(1+2+3, 4+5+6, 7+8+9),  # sums of gene counts per 
    stringsAsFactors = FALSE
  )
  
  # 2) gene_count pseudo
  gene_count <- data.frame(
    SPATIAL_1 = c(1, 2, 3),
    SPATIAL_2 = c(4, 5, 6),
    SPATIAL_3 = c(7, 8, 9),
    row.names     = c("gA", "gB", "gC"),
    stringsAsFactors = FALSE,
    check.names   = FALSE
  )
  
  # 3) Tests on matched_data
  expect_s3_class(matched_data, "data.frame")
  expect_equal(nrow(matched_data), 3)
  expect_true(all(c(
    "barcode_sequence", "spatial_name",
    "X_coordinate", "Y_coordinate", "UMI_count"
  ) %in% colnames(matched_data)))
  # Check UMI counts
  expect_equal(matched_data$UMI_count, c(6, 15, 24))
  
  # 4) Tests on gene_count
  expect_s3_class(gene_count, "data.frame")
  expect_equal(nrow(gene_count), 3)
  expect_equal(rownames(gene_count), c("gA", "gB", "gC"))
  expect_true(all(grepl("^SPATIAL_", colnames(gene_count))))
  expect_equal(gene_count["gB", "SPATIAL_2"], 5)
  
})