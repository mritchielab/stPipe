#####################################################
# Visualization for sST upstream results
#####################################################

#' @name Run_Visualization
#' @title Visualize pre-processed sST data
#' @description This function visualizes both spatial-level and read-level information either before or after the 'Run_QC' step. It outputs raw and log-transformed UMI count plots for spatial flag, and demultiplexing and mapping statistics for read flag.
#' @param config Path to the YAML configuration file.
#' @param matched.data A data frame containing spatial transcriptomics data, including UMI counts and spatial coordinates. This is usually obtained from 'Run_loc_match' function.
#' @param Vis.spatial Logical value indicating whether to visualize spatial data. Defaults to TRUE.
#' @param Vis.read Logical value indicating whether to visualize read-level data. Defaults to TRUE.
#' @param show.config Logical value indicating whether to print the configuration. Defaults to TRUE.
#' @return A list contains spatial and read level visualization results.
#' @examples
#' temp_config <- tempfile(fileext = ".yml")
#' writeLines("technology_version: 'Visium 1.0'\noutput_directory: '.'", temp_config)
#' matched.data <- data.frame(
#' X_coordinate = runif(10, 0, 100),
#' Y_coordinate = runif(10, 0, 100),
#' UMI_count = sample(seq_len(100), 10),
#' spatial_name = paste0("Spot", seq_len(10)),
#' stringsAsFactors = FALSE
#' )
#' vis_results <- Run_Visualization(
#'   matched.data = matched.data, 
#'   config = temp_config, 
#'   Vis.spatial = TRUE, 
#'   Vis.read = FALSE, 
#'   show.config = FALSE
#' )
#' @export
#' @importFrom ggplot2 ggplot scale_fill_brewer aes geom_bar geom_text theme_minimal theme labs element_text geom_point scale_color_gradient element_blank element_rect xlim ylim

Run_Visualization <- function(matched.data = NULL, config, Vis.spatial = TRUE, Vis.read = TRUE, show.config = TRUE) {

  config <- yaml::read_yaml(config)
  if (show.config) {
    print(config)
  }
  technology_version <- as.character(config$technology_version)
  output_dir <- as.character(config$output_directory)

  # Visualization of spatial data
  if (Vis.spatial && !is.null(matched.data)) {
    # adjust coordination system for Visium
    if (grepl("Visium", technology_version)) {
      matched.data$Y_coordinate <- -matched.data$Y_coordinate
    }

    # Set plotting limits
    x_range <- range(matched.data$X_coordinate)
    y_range <- range(matched.data$Y_coordinate)
    umi_range <- range(matched.data$UMI_count)

    # Raw UMI count plot
    raw_umi_plot <- ggplot(matched.data, aes(x = X_coordinate, y = Y_coordinate, color = UMI_count)) +
      geom_point(size = 3) +
      scale_color_gradient(low = "white", high = "red") +
      labs(title = 'Spatial Heatmap of UMI Count [raw]', x = 'X Coordinate', y = 'Y Coordinate', color = 'UMI Count') +
      theme_minimal(base_size = 15) +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "gray"),
        plot.background = element_rect(fill = "gray"),
        axis.text = element_blank(),
        axis.title = element_blank()
      )

    # Apply plotting limits if available
    if (!is.null(x_range)) raw_umi_plot <- raw_umi_plot + xlim(x_range)
    if (!is.null(y_range)) raw_umi_plot <- raw_umi_plot + ylim(y_range)

    # Log-transformed UMI count plot
    log_umi_plot <- ggplot(matched.data, aes(x = X_coordinate, y = Y_coordinate, color = log(UMI_count))) +
      geom_point(size = 3) +
      scale_color_gradient(low = "white", high = "red") +
      labs(title = 'Spatial Heatmap of UMI Count [log]', x = 'X Coordinate', y = 'Y Coordinate', color = 'UMI Count [log]') +
      theme_minimal(base_size = 15) +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "gray"),
        plot.background = element_rect(fill = "gray"),
        axis.text = element_blank(),
        axis.title = element_blank()
      )

    # Apply plotting limits if available
    if (!is.null(x_range)) log_umi_plot <- log_umi_plot + xlim(x_range)
    if (!is.null(y_range)) log_umi_plot <- log_umi_plot + ylim(y_range)
  }

  # Visualization of read-level data
  if (Vis.read) {
    # Read in data from the working directory
    sce <- scPipe::create_sce_by_dir(output_dir)
    sce <- scPipe::calculate_QC_metrics(sce)

    # Demultiplexing plot
    demultiplex_info <- scPipe::demultiplex_info(sce)[seq_len(6), ]
    demultiplex_info$count <- as.numeric(demultiplex_info$count)
    total_count <- sum(demultiplex_info$count)
    demultiplex_info$percentage <- (demultiplex_info$count / total_count) * 100

    demultiplex_plot <- ggplot(demultiplex_info, aes(x = reorder(status, -count), y = count, fill = status)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = paste0(round(percentage, 2), "%")), vjust = -0.5, size = 4) +
      scale_fill_brewer(palette = "Set3") +
      labs(title = "Barcode Demultiplexing Information", x = "Status", y = "Percentage") +
      theme_minimal(base_size = 15) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position = "none"
      )

    # UMI duplication rate plot
    duplication_plot <- scPipe::plot_UMI_dup(sce)

    # Mapping statistics plot
    mapping_plot <- scPipe::plot_mapping(sce, percentage = TRUE, dataname = "sample data")
  }

  return(list(
    raw_UMI_plot = if (exists("raw_umi_plot")) raw_umi_plot else NULL,
    log_UMI_plot = if (exists("log_umi_plot")) log_umi_plot else NULL,
    demultiplex_plot = if (exists("demultiplex_plot")) demultiplex_plot else NULL,
    duplication_plot = if (exists("duplication_plot")) duplication_plot else NULL,
    mapping_plot = if (exists("mapping_plot")) mapping_plot else NULL
  ))
}
