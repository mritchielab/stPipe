#' Match spatial location After Upstream Pre-Processing for Sequencing-Based Spatial Transcriptomics
#'
#' This function matches spatial coordinates from sequencing-based spatial transcriptomics data after upstream pre-processing.
#'
#' @param config Path to the YAML configuration file.
#' @param show.config Logical value indicating whether to print the configuration. Defaults to TRUE.
#' @return A data frame contains filtered gene counts and spatial coordinates
#' @export
#' @examples matching <- Run_loc_match(config = "~/Desktop/config_stPipe.yml", show.config = TRUE)

Run_loc_match <- function(config, show.config = TRUE) {

  # Read configuration from the provided YAML file
  config <- yaml::read_yaml(config)
  # Conditionally print config if show.config is TRUE
  if (show.config) {
    print(config)
  }
  output_directory <- as.character(config$output_directory)
  technology_version <- as.character(config$technology_version)
  visium_coordination <- as.character(config$visium_coordination)
  bead_location <- as.character(config$bead_location)

  # Set the working directory to the output directory
  setwd(output_directory)

  # Create a new directory named 'upstream_QC' within the output directory
  qc_dir <- file.path(output_directory, "upstream_QC")
  if (!dir.exists(qc_dir)) {
    dir.create(qc_dir)
  }
  setwd(qc_dir)

  # Function to construct the file path for coordinates based on the technology version
  get_coordinates_file <- function(visium_coordination) {
    file_name <- switch(visium_coordination,
                        "V1" = stPipe::visium_v1_coordinates,
                        "V2" = stPipe::visium_v2_coordinates,
                        "V3" = stPipe::visium_v3_coordinates,
                        "V4" = stPipe::visium_v4_coordinates,
                        "V5" = stPipe::visium_v5_coordinates,
                        stop("Unsupported technology version"))
    return(file_name)
  }

  # Match extracted spatial barcode to x & y coordination
  if (grepl("Visium", technology_version)) {
    coordinates <- get_coordinates_file(visium_coordination)
    print("Coordinates for Visium:")
    print(coordinates)
  } else if (technology_version == "Slideseq" || technology_version == "Curio-seeker") {
    coordinates <- read.csv(bead_location, header = TRUE, stringsAsFactors = FALSE)
    print("Coordinates for Slideseq or Curio-seeker:")
    print(coordinates)
  } else {
    stop("Unknown technology_version in config file.")
  }
  colnames(coordinates) <- c("barcode_sequence", "X_coordinate", "Y_coordinate")

  # Read annotation and gene count data
  annotation_file <- file.path(output_directory, "sample_index.csv")
  annotation <- read.csv(annotation_file, stringsAsFactors = FALSE)
  matched_data <- merge(annotation, coordinates, by = "barcode_sequence", all.x = TRUE)
  matched_data <- matched_data[complete.cases(matched_data), ]

  gene_count_file <- file.path(output_directory, "gene_count.csv")
  gene_count <- read.csv(gene_count_file)
  sum_counts <- colSums(gene_count[, -1])
  matched_data$UMI_count <- as.numeric(sum_counts[match(matched_data$cell_name, names(sum_counts))])
  write.csv(matched_data, file = "matched_spatial_barcode.csv", row.names = FALSE)
  return(matched_data)
}
