#' Run_Clustering: Perform Basic Clustering Algorithms
#'
#' This function performs basic clustering using t-SNE and UMAP, either before or after the QC step and outputs visualizations to output directory.
#'
#' @param gene.count A data frame containing gene count data with gene IDs as row names.
#' @param matched.data A data frame containing matched spatial coordinates with raw UMI counts.
#' @export
#' @importFrom ggplot2 ggplot scale_fill_brewer aes geom_bar geom_text theme_minimal theme labs ggsave element_text geom_point scale_color_gradient element_blank element_rect xlim ylim
#' @examples clustering.result <- Run_Clustering(gene.count = gene_count, matched.data =  matching)

Run_Clustering <- function(gene.count, matched.data) {

  # Transpose the gene count data
  gene_count_t <- t(gene.count)
  # Remove duplicates before running t-SNE and UMAP
  gene_count_t <- unique(gene_count_t)
  # Set seed for reproducibility
  set.seed(42)
  # Perform UMAP
  umap_results <- umap::umap(gene_count_t)
  # Perform t-SNE
  tsne_result <- Rtsne::Rtsne(gene_count_t, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 5000, partial_pca = TRUE)

  # Create data frames for UMAP and t-SNE results
  umap_df <- data.frame(
    UMAP1 = umap_results$layout[, 1],
    UMAP2 = umap_results$layout[, 2],
    spot = rownames(gene_count_t)
  )
  tsne_df <- data.frame(
    TSNE1 = tsne_result$Y[, 1],
    TSNE2 = tsne_result$Y[, 2],
    spot = rownames(gene_count_t)
  )

  # Merge UMAP and t-SNE results with matched data
  merged_data_umap <- merge(umap_df, matched.data, by.x = "spot", by.y = "cell_name", all.x = TRUE)
  merged_data_tsne <- merge(tsne_df, matched.data, by.x = "spot", by.y = "cell_name", all.x = TRUE)

  # Generate and save the t-SNE plot
  tsne_plot <- ggplot(merged_data_tsne, aes(x = -TSNE1, y = TSNE2, color = UMI_count)) +
    geom_point(size = 0.5) +
    scale_color_gradient(low = "#FFEEEE", high = "red") +
    labs(title = 't-SNE Plot of UMI Count', x = 't-SNE1', y = 't-SNE2', color = 'UMI Count') +
    theme_minimal(base_size = 15) +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5)
    )

  # Generate and save the UMAP plot
  umap_plot <- ggplot(merged_data_umap, aes(x = UMAP1, y = UMAP2, color = UMI_count)) +
    geom_point(size = 0.5) +
    scale_color_gradient(low = "#FFEEEE", high = "red") +
    labs(title = 'UMAP Plot of UMI Count', x = 'UMAP1', y = 'UMAP2', color = 'UMI Count') +
    theme_minimal(base_size = 15) +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5)
    )

  # return both data frames and plots for TSNE and UMAP
  return(list(tsne_df = tsne_df, tsne_plot = tsne_plot, umap_df = umap_df, umap_plot = umap_plot))

}
