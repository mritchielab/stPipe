#####################################################
# Perform data clustering
#####################################################

#' @name Run_Clustering
#' @title Perform Basic Clustering Algorithms
#' @description This function performs basic clustering using t-SNE and UMAP, either before or after the QC step and outputs visualizations results.
#'
#' @details
#' T-SNE (t-Distributed Stochastic Neighbor Embedding)
#'
#' T-SNE is a popular dimensionality reduction technique used to project high-dimensional data into two or three dimensions.
#' It optimizes for preserving the local structure of the data, meaning that points that are close to each other in high-dimensional space will remain close in the lower-dimensional representation.
#' This method is widely used in single-cell RNA-seq and spatial transcriptomics to explore the heterogeneity of the data and identify cells or spatial regions with similar expression profiles.
#'
#' UMAP (Uniform Manifold Approximation and Projection)
#'
#' UMAP is another popular dimensionality reduction technique that aims to preserve both the local and global structures of the data, with a stronger emphasis on efficiency.
#' It excels at handling large datasets, generating low-dimensional projections in less time compared to T-SNE.
#' Unlike T-SNE, UMAP seeks to retain both local and global data structures, making it more suitable for capturing broad clustering patterns across the datasets.
#' @param gene.count A data frame containing gene count data with gene IDs as row names.
#' @param matched.data A data frame containing matched spatial coordinates with raw UMI counts.
#' @param num_clusters Number of clusters during clustering. Default set to five.
#' @return A list which contains interactive clustering visualization result with related data frames.
#' @examples
#' set.seed(123)
#' gene.count <- matrix(sample(0:100, 200 * 100, replace = TRUE), nrow = 200)
#' rownames(gene.count) <- paste0("Gene", seq_len(200))
#' colnames(gene.count) <- paste0("Spot", seq_len(100))
#' matched.data <- data.frame(
#'   spatial_name = paste0("Spot", seq_len(100)),
#'   x_coord = runif(100, 0, 10),
#'   y_coord = runif(100, 0, 10)
#' )
#' result <- Run_Clustering(
#'   gene.count = gene.count, 
#'   matched.data = matched.data, 
#'   num_clusters = 3
#' )
#' @export
#' @importFrom ggplot2 ggplot scale_fill_brewer aes geom_bar geom_text theme_minimal theme labs ggsave element_text geom_point scale_color_gradient element_blank element_rect xlim ylim scale_color_brewer

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("TSNE1","TSNE2","UMAP1","UMAP2","cluster")
  )
}

Run_Clustering <- function(gene.count, matched.data, num_clusters = 5) {

  gene_count_t <- t(gene.count)
  # Remove duplicates before running t-SNE and UMAP
  gene_count_t <- unique(gene_count_t)

  # UMAP & t-SNE
  umap_results <- umap::umap(gene_count_t)
  tsne_result <- Rtsne::Rtsne(gene_count_t, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 1000, partial_pca = TRUE)

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

  # K-means clustering on UMAP results
  kmeans_result_umap <- stats::kmeans(umap_results$layout, centers = 5)
  kmeans_result_tsne <- stats::kmeans(tsne_result$Y, centers = 5)

  # Assign cluster to UMAP and t-SNE data frames
  umap_df$cluster <- as.factor(kmeans_result_umap$cluster)
  tsne_df$cluster <- as.factor(kmeans_result_tsne$cluster)  # apply clustering results to t-SNE

  # Merge UMAP and t-SNE results with matched data
  merged_data_umap <- merge(umap_df, matched.data, by.x = "spot", by.y = "spatial_name", all.x = TRUE)
  merged_data_tsne <- merge(tsne_df, matched.data, by.x = "spot", by.y = "spatial_name", all.x = TRUE)

  # Generate and save the t-SNE & UMAP plot with cluster labels
  tsne_plot <- ggplot(merged_data_tsne, aes(x = -TSNE1, y = TSNE2, color = cluster)) +
    geom_point(size = 0.5) +
    scale_color_brewer(palette = "Set1") +  # color for different clusters
    labs(title = 't-SNE Plot of Clusters', x = 't-SNE1', y = 't-SNE2', color = 'Cluster') +
    theme_minimal(base_size = 15) +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5)
    )

  umap_plot <- ggplot(merged_data_umap, aes(x = UMAP1, y = UMAP2, color = cluster)) +
    geom_point(size = 0.5) +
    scale_color_brewer(palette = "Set1") +  # color for different clusters
    labs(title = 'UMAP Plot of Clusters', x = 'UMAP1', y = 'UMAP2', color = 'Cluster') +
    theme_minimal(base_size = 15) +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5)
    )

  # interactive plot
  tsne_interactive <- plotly::ggplotly(tsne_plot)
  umap_interactive <- plotly::ggplotly(umap_plot)

  # return data frames and plots for TSNE and UMAP
  return(list(tsne_df = tsne_df, tsne_plot = tsne_plot, umap_df = umap_df, umap_plot = umap_plot, tsne_interactive = tsne_interactive, umap_interactive = umap_interactive))
}
