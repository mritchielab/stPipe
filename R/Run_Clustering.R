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
#' @param n_features Number of features used during clustering. Default set to one thousand.
#' @return A list which contains interactive clustering visualization result with related data frames.
#' @examples
#' set.seed(123)
#' gene.count <- matrix(sample(0:100, 200 * 100, replace = TRUE), nrow = 200)
#' rownames(gene.count) <- paste0("Gene", 1:200)
#' colnames(gene.count) <- paste0("Spot", 1:100)
#' matched.data <- data.frame(
#'   spatial_name = paste0("Spot", 1:100),
#'   x_coord = runif(100, 0, 10),
#'   y_coord = runif(100, 0, 10)
#' )
#' result <- Run_Clustering(gene.count = gene.count, matched.data = matched.data, num_clusters = 3, n_features = 10)
#' @export
#' @importFrom ggplot2 ggplot scale_fill_brewer aes geom_bar geom_text theme_minimal theme labs ggsave element_text geom_point scale_color_gradient element_blank element_rect xlim ylim scale_color_brewer

Run_Clustering <- function(gene.count, matched.data, num_clusters = 5, n_features = 1000) {
  # Feature selection: calculate gene-wise mean, variance, and dispersion
  gene_mean <- rowMeans(gene.count)
  gene_var <- apply(gene.count, 1, var)
  
  # Filter out genes with mean = 0
  valid_genes <- gene_mean > 0
  gene_mean <- gene_mean[valid_genes]
  gene_var <- gene_var[valid_genes]
  
  # Calculate dispersion: variance / mean
  gene_disp <- gene_var / gene_mean
  gene_disp[is.na(gene_disp) | is.infinite(gene_disp)] <- 0
  
  # Select top n_features genes by dispersion
  selected_genes <- names(sort(gene_disp, decreasing = TRUE))[1:min(n_features, length(gene_disp))]
  
  # Subset gene count matrix
  gene.count <- gene.count[selected_genes, ]
  
  # ------------------------
  # Transpose expression matrix so rows represent spots; remove duplicates
  gene_count_t <- t(gene.count)
  gene_count_t <- unique(gene_count_t)
  
  # ------------------------
  # Dimensionality reduction using UMAP and t-SNE
  umap_results <- umap::umap(gene_count_t)
  tsne_result <- Rtsne::Rtsne(gene_count_t, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 1000, pca = FALSE)
  
  # Construct UMAP and t-SNE data frames
  umap_df <- data.frame(
    UMAP1 = umap_results$layout[, 1],
    UMAP2 = umap_results$layout[, 2],
    spot = rownames(gene_count_t),
    stringsAsFactors = FALSE
  )
  
  tsne_df <- data.frame(
    TSNE1 = tsne_result$Y[, 1],
    TSNE2 = tsne_result$Y[, 2],
    spot = rownames(gene_count_t),
    stringsAsFactors = FALSE
  )
  
  # ------------------------
  # K-means clustering on dimensionality reduction results
  kmeans_result_umap <- stats::kmeans(umap_results$layout, centers = num_clusters)
  kmeans_result_tsne <- stats::kmeans(tsne_result$Y, centers = num_clusters)
  
  umap_df$cluster <- as.factor(kmeans_result_umap$cluster)
  tsne_df$cluster <- as.factor(kmeans_result_tsne$cluster)
  
  # ------------------------
  # Merge with spatial coordinate data
  merged_data_umap <- merge(umap_df, matched.data, by.x = "spot", by.y = "spatial_name", all.x = TRUE)
  merged_data_tsne <- merge(tsne_df, matched.data, by.x = "spot", by.y = "spatial_name", all.x = TRUE)
  
  # ------------------------
  # Generate t-SNE plot
  tsne_plot <- ggplot(merged_data_tsne, aes(x = -TSNE1, y = TSNE2, color = cluster)) +
    geom_point(size = 0.5) +
    scale_color_brewer(palette = "Set1") +
    labs(title = 't-SNE Plot of Clusters', x = 't-SNE1', y = 't-SNE2', color = 'Cluster') +
    theme_minimal(base_size = 15) +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          plot.title = element_text(hjust = 0.5))
  
  # Generate UMAP plot
  umap_plot <- ggplot(merged_data_umap, aes(x = UMAP1, y = UMAP2, color = cluster)) +
    geom_point(size = 0.5) +
    scale_color_brewer(palette = "Set1") +
    labs(title = 'UMAP Plot of Clusters', x = 'UMAP1', y = 'UMAP2', color = 'Cluster') +
    theme_minimal(base_size = 15) +
    theme(panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          plot.title = element_text(hjust = 0.5))
  
  # Convert to interactive plots
  tsne_interactive <- plotly::ggplotly(tsne_plot)
  umap_interactive <- plotly::ggplotly(umap_plot)
  
  # Return result list
  return(list(
    tsne_df = tsne_df,
    tsne_plot = tsne_plot,
    umap_df = umap_df,
    umap_plot = umap_plot,
    tsne_interactive = tsne_interactive,
    umap_interactive = umap_interactive
  ))
}
