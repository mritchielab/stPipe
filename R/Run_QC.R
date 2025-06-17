#####################################################
# Perform QC filtering for sST preliminary results
#####################################################

#' @name Run_QC
#' @title QC Control After Upstream Pre-Processing for Sequencing-Based Spatial Transcriptomics
#' @details
#'
#' This function performs QC control on sequencing-based spatial transcriptomics data after upstream pre-processing step such as 'Run_ST' step. Ensure the output directory is the same with the 'Run_ST' one.
#' Filtering is performed either use specific UMI threshold or assign the threshold to 'DropletUtils'.
#'
#' "max_slope"
#'
#' In this approach, the filtering is done based on UMI counts.
#' Spots with counts below a certain threshold are considered low-quality and are filtered out.
#' This method helps retain only the spots with significant transcriptomic signals, reducing noise
#' from sptos with minimal or no meaningful biological information.
#'
#' Threshold Determination:
#' The threshold in this method is computed by analyzing the distribution of UMI counts across spots,
#' and identifying the point of maximum slope in the cumulative UMI distribution curve.
#' This point often corresponds to the transition between background noise and real biological signals.
#'
#' "EmptyDropletUtils"
#'
#' Alternatively, the \code{DropletUtils} package offers a more sophisticated approach
#' by using statistical methods to identify droplets or spots that contain real cells, as opposed
#' to empty droplets or those containing background RNA. This method calculates a
#' false discovery rate (FDR) to assess the likelihood of each droplet containing a real cell. Sptos are retained if they meet the significance criteria for either the p-value or FDR.
#' To learn more details regarding \code{DropletUtils}, visit
#' \href{https://bioconductor.org/packages/release/bioc/html/DropletUtils.html}{this link to its Bioconductor page}.
#'
#' Multiple Thresholds:
#' This method will determine two thresholds based on the config file input parameters.
#' Filtering can be fine-tuned using both p-value and FDR thresholds, offering greater flexibility in distinguishing between noise and meaningful data.
#' Sptos are retained if they meet the significance criteria for either the p-value or FDR.
#' @param config Path to the YAML configuration file.
#' @param matched.data A data frame containing spatial transcriptomics data, including UMI counts and spatial coordinates, this is usually obtained from 'Run_loc_match' function.
#' @param gene.matrix A gene count matrix, this is usually obtained from 'Run_ST' function.
#' @param show.config Logical value indicating whether to print the configuration. Defaults to TRUE.
#' @export
#' @return A list containing filtered gene counts with matched spatial coordinates after QC.
#' @examples
#' output_dir <- tempdir()
#' config_list <- list(
#' output_directory = output_dir,
#' qc_filter = "slope_max",
#' qc_per = "0.4_0.8" 
#' )
#' config_path <- tempfile(fileext = ".yml")
#' yaml::write_yaml(config_list, config_path)
#' set.seed(123)
#' gene_ids <- paste0("gene", seq_len(100))
#' spatial_names <- paste0("SPATIAL_", seq_len(100))
#' count_matrix <- matrix(rpois(100*100, lambda = 20), nrow = 100, ncol = 100)
#' colnames(count_matrix) <- spatial_names
#' gene_matrix <- data.frame(row.names = gene_ids, count_matrix, stringsAsFactors = FALSE)
#' matched_data <- data.frame(
#'  X_coordinate = runif(100, min = 0, max = 1000),
#'  Y_coordinate = runif(100, min = 0, max = 1000),
#'  barcode_sequence = paste0("BC", seq_len(100)),
#'  spatial_name = spatial_names,
#'  stringsAsFactors = FALSE
#' )
#' umi_counts <- colSums(gene_matrix[, -1])
#' matched_data$UMI_count <- umi_counts[match(matched_data$spatial_name, names(umi_counts))]
#' qc_results <- Run_QC(
#'  config = config_path,
#'  matched.data = matched_data,
#'  gene.matrix = gene_matrix,
#'  show.config = FALSE
#' )
#' @importFrom ggplot2 ggplot scale_fill_brewer aes geom_bar geom_text theme_minimal theme labs ggsave element_text geom_point scale_color_gradient element_blank element_rect xlim ylim geom_segment geom_line geom_hline arrow unit
#' @importFrom dplyr %>%
#' @importFrom stats quantile na.omit


Run_QC <- function(config, matched.data, gene.matrix, show.config = TRUE) {

  config <- yaml::read_yaml(config)
  if (show.config) {
    message(config)
  }

  output_directory <- as.character(config$output_directory)
  qc_filter <- as.character(config$qc_filter)
  qc_per <- unlist(strsplit(config$qc_per, "_"))
  qc_per_lower <- as.numeric(qc_per[1])
  qc_per_retain <- as.numeric(qc_per[2])

  # New directory named 'upstream_QC' under output directory
  qc_dir <- file.path(output_directory, "upstream_QC")
  if (!dir.exists(qc_dir)) {
    dir.create(qc_dir)
  }
  gene_count <- gene.matrix
  sum_counts <- colSums(gene_count[, -1])
  spatial_coords <- matched.data[, c("X_coordinate", "Y_coordinate", "barcode_sequence", "spatial_name")]

  # Barcode quality filtering
  quantiles <- quantile(sum_counts, probs = seq(0, 1, 0.01))
  quantile_df <- data.frame(
    Percentile = seq(0, 100, 1),
    Value = as.numeric(quantiles)
  )
  quantile_df$Slope <- c(NA, diff(quantile_df$Value) / diff(quantile_df$Percentile))
  max_slope_index <- which.max(quantile_df$Slope[-1]) + 1  # first slope is NA
  max_slope_interval <- quantile_df[max_slope_index - seq_len(max_slope_index), ]
  max_slope_interval <- max_slope_interval[max_slope_interval$Slope < quantile(na.omit(quantile_df$Slope), qc_per_retain)[[1]], ]
  max_slope_interval <- max_slope_interval[max_slope_interval$Value < quantile(na.omit(quantile_df$Value), qc_per_retain)[[1]], ]
  max_slope_interval2 <- max_slope_interval
  max_slope_interval2 <- max_slope_interval2[max_slope_interval2$Slope < quantile(na.omit(quantile_df$Slope), qc_per_lower)[[1]], ]
  max_slope_interval2 <- max_slope_interval2[max_slope_interval2$Value < quantile(na.omit(quantile_df$Value), qc_per_lower)[[1]], ]
  max_slope_interval$SlopeIncrease <- -c(diff(max_slope_interval$Slope), NA)
  max_slope_interval2$SlopeIncrease <- -c(diff(max_slope_interval2$Slope), NA)
  threshold <- max_slope_interval$Value[which.max(max_slope_interval$SlopeIncrease)]
  ed_lower <- max_slope_interval2$Value[which.max(max_slope_interval2$SlopeIncrease)]

  # Plotting
  suppressWarnings({
    p <- ggplot(quantile_df, aes(x = Percentile, y = Value)) +
      geom_line() +
      geom_point() +
      geom_segment(data = max_slope_interval, 
                   aes(x = Percentile[1], y = Value[1], xend = Percentile[2], yend = Value[2]),
                   color = "red", linewidth = 1, arrow = arrow(type = "closed", length = unit(0.2, "inches"))) +
      geom_text(data = max_slope_interval, 
                aes(x = Percentile[2], y = Value[2], label = paste("EmptyDrop Retain Threshold", round(Value[1], 2))),
                vjust = -1, color = "red") +
      geom_hline(yintercept = ed_lower, linetype = "dashed", color = "blue") +
      geom_text(aes(x = max(Percentile), y = ed_lower, 
                    label = paste("EmptyDrop Lower Threshold", round(ed_lower, 2))),
                vjust = -1, hjust = 1.1, color = "blue") +
      labs(title = "Threshold to filter out low count UMI as potential background noise", 
           x = "UMI Count Percentile", y = "Value") +
      theme_minimal()
    
    ggsave("Threshold.pdf", plot = p, width = 10, height = 8, units = "in", device = "pdf")  
  })
  
  # Filtering based on QC method
  if (qc_filter == "slope_max") {
    keep_cells <- sum_counts[sum_counts > threshold]
    filtered_gene_count <- gene_count[, names(keep_cells)]
  } else if (qc_filter == "EmptyDropletUtils") {
    keep_cells <- sum_counts[sum_counts > threshold]
    gene_count_ED <- gene_count
    rownames(gene_count_ED) <- gene_count_ED$gene_id
    gene_count_ED$gene_id <- NULL
    outs <- DropletUtils::emptyDrops(gene_count_ED, retain = threshold)
    is.cell <- outs$FDR <= 0.01
    is.sign <- outs$PValue <= 0.05
    is.both <- is.cell & is.sign
    filtered_gene_count <- gene_count_ED[, which(is.both), drop = FALSE]
  }

  spatial_coords <- spatial_coords[spatial_coords$spatial_name %in% names(keep_cells) & spatial_coords$spatial_name %in% colnames(filtered_gene_count), ]
  spatial_coords$UMI_count <- sum_counts[spatial_coords$spatial_name]

  return(list(filtered_gene_count = filtered_gene_count, filtered_spatial_coords = spatial_coords))
}
