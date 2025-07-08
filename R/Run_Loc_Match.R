#####################################################
# Perform spatial location matching for gene count matrix
#####################################################

#' @name Run_Loc_Match
#' @title Match spatial location After Pre-Processing for Sequencing-Based Spatial Transcriptomics
#'
#' @description This function matches spatial coordinates for sST data after upstream pre-processing.
#' 'Run_loc_match' can either map the technology coordination system (such as spot in Visium coordination, bead in Slide-seq coordination)
#' or compute pixel for each spot and map the pixel information back to the image (only for Visium)
#'
#' @param config Path to the YAML configuration file.
#' @param pixel Computing spot pixel or not. If yes, compute pixel for each spot and map back to image; if not, map the Visium coordination system. Defaults to FALSE.
#' @param show.config Logical value indicating whether to print the configuration. Defaults to TRUE.
#' @return A data frame contains gene count matrix with spatial coordinates
#' @examples
#' data_dir <- tempdir()
#' output_dir <- file.path(tempdir(), "Run_Loc_Match_output")
#' if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
#' sample_index <- data.frame(
#'   barcode_sequence = c("BC001", "BC002", "BC003"),
#'   cell_name = c("CELL_1", "CELL_2", "CELL_3"),
#'   stringsAsFactors = FALSE
#' )
#' write.csv(sample_index, file = file.path(output_dir, "sample_index.csv"), row.names = FALSE)
#' gene_count <- data.frame(
#'   gene_id = c("gene1", "gene2"),
#'   CELL_1 = c(10, 5),
#'   CELL_2 = c(20, 10),
#'   CELL_3 = c(30, 15),
#'   stringsAsFactors = FALSE
#' )
#' write.csv(gene_count, file = file.path(output_dir, "gene_count.csv"), row.names = FALSE)
#' config_list <- list(
#'   output_directory = output_dir,
#'   data_directory = data_dir,
#'   technology_version = "Visium_probe_v1",
#'   visium_coordination = "V1"
#' )
#' config_file <- tempfile(fileext = ".yml")
#' yaml::write_yaml(config_list, config_file)
#' result <- Run_Loc_Match(
#'   config = config_file,
#'   pixel = FALSE,
#'   show.config = FALSE
#' )
#' @export

Run_Loc_Match <- function(config, pixel = FALSE, show.config = TRUE) {

  config <- yaml::read_yaml(config)
  if (show.config) {
    message(config)
  }

  output_directory <- as.character(config$output_directory)
  data_directory <- as.character(config$data_directory)
  technology_version <- as.character(config$technology_version)
  visium_coordination <- as.character(config$visium_coordination)
  image_path <- as.character(config$image_path)
  bead_location <- as.character(config$bead_location)

  # New directory named 'upstream_QC' under the output directory
  qc_dir <- file.path(output_directory, "upstream_QC")
  if (!dir.exists(qc_dir)) {
    dir.create(qc_dir)
  }

  # Match extracted spatial barcode to x & y coordination
  if (grepl("Visium", technology_version)) {
    coordinates <- switch(
      visium_coordination,
      "V1" = readRDS(system.file("extdata", "visium_v1_coordinates.rds", package = "stPipe")),
      "V2" = readRDS(system.file("extdata", "visium_v2_coordinates.rds", package = "stPipe")),
      "V3" = readRDS(system.file("extdata", "visium_v3_coordinates.rds", package = "stPipe")),
      "V4" = readRDS(system.file("extdata", "visium_v4_coordinates.rds", package = "stPipe")),
      "V5" = readRDS(system.file("extdata", "visium_v5_coordinates.rds", package = "stPipe")),
      stop("Unsupported technology version")
    )
    message("Heading Coordinates for Visium:")
    message(utils::head(coordinates))
    message("Start mapping coordination...")
    colnames(coordinates) <- c("barcode_sequence", "X_coordinate", "Y_coordinate")
    annotation_file <- file.path(output_directory, "sample_index.csv")
    annotation <- utils::read.csv(annotation_file, stringsAsFactors = FALSE)
    matched_data <- merge(annotation, coordinates, by = "barcode_sequence", all.x = TRUE)
    matched_data <- matched_data[stats::complete.cases(matched_data), ]

    gene_count_file <- file.path(output_directory, "gene_count.csv")
    gene_count <- utils::read.csv(gene_count_file)
    sum_counts <- colSums(gene_count[, -1])
    matched_data$UMI_count <- as.numeric(sum_counts[match(matched_data$cell_name, names(sum_counts))])
  } else if (technology_version == "Slideseq" || technology_version == "Curio-seeker") {
    coordinates <- utils::read.csv(bead_location, header = TRUE, stringsAsFactors = FALSE)
    message("Heading Coordinates for Slideseq or Curio-seeker:")
    message(utils::head(coordinates))
    message("Start mapping coordination...")
    colnames(coordinates) <- c("barcode_sequence", "X_coordinate", "Y_coordinate")

    # Read annotation and gene count data
    annotation_file <- file.path(output_directory, "sample_index.csv")
    annotation <- utils::read.csv(annotation_file, stringsAsFactors = FALSE)
    matched_data <- merge(annotation, coordinates, by = "barcode_sequence", all.x = TRUE)
    matched_data <- matched_data[stats::complete.cases(matched_data), ]

    gene_count_file <- file.path(output_directory, "gene_count.csv")
    gene_count <- utils::read.csv(gene_count_file)
    sum_counts <- colSums(gene_count[, -1])
    matched_data$UMI_count <- as.numeric(sum_counts[match(matched_data$cell_name, names(sum_counts))])

  } else if (technology_version == "Stereoseq") {
    message("Processing Stereoseq technology...")

    # Load sample_index
    sample_index_path <- file.path(output_directory, "sample_index.csv")
    sample_index <- data.table::fread(sample_index_path)

    # Initialize x and y columns
    sample_index$x <- NA
    sample_index$y <- NA

    # Process R1 file
    r1_path <- list.files(data_directory, pattern = "R1", full.names = TRUE)
    decompressed_r1 <- tempfile(fileext = ".fq")
    system(paste("gzip -dc", shQuote(r1_path), ">", decompressed_r1))

    r1_content <- data.table::fread(decompressed_r1, sep = "\n", header = FALSE)
    r1_headers <- r1_content[seq(1, nrow(r1_content), 4), ]
    r1_barcodes <- r1_content[seq(2, nrow(r1_content), 4), ]

    r1_table <- data.frame(
      read_id = sub("^@(.*)/.*$", "\\1", r1_headers[[1]]),
      barcode = substr(r1_barcodes[[1]], 1, 25)
    )

    # Process demultiplexed file
    demultiplexed_path <- file.path(output_directory, "demultiplexed.fq")
    demultiplexed <- data.table::fread(demultiplexed_path, sep = "\n", header = FALSE)
    demultiplexed_headers <- demultiplexed[seq(1, nrow(demultiplexed), 4), ]

    demultiplexed_table <- do.call(rbind, lapply(demultiplexed_headers[[1]], function(line) {
      read_id <- sub("^@(.*)[:_].*$", "\\1", line)
      coords <- regmatches(line, gregexpr("x\\d+|y\\d+", line, perl = TRUE))[[1]]
      if (length(coords) == 2) {
        return(data.frame(
          read_id = read_id,
          x = as.integer(sub("x", "", coords[1])),
          y = as.integer(sub("y", "", coords[2]))
        ))
      } else {
        return(data.frame(read_id = read_id, x = NA, y = NA))
      }
    }))

    # Merge R1 and demultiplexed
    matched_table <- merge(
      r1_table,
      demultiplexed_table,
      by = "read_id",
      all.x = TRUE
    )

    # Use read_id to match barcode_sequence in sample_index
    final_table <- merge(
      sample_index,
      matched_table[, c("barcode", "x", "y")],
      by.x = "barcode_sequence",
      by.y = "barcode",
      all.x = TRUE
    )

    # Process gene count
    gene_count_file <- file.path(output_directory, "gene_count.csv")
    gene_count <- data.table::fread(gene_count_file)
    umi_sums <- colSums(gene_count[, -1])
    final_table$UMI_count <- umi_sums[match(final_table$cell_name, names(umi_sums))]

    # Update column names
    colnames(final_table)[2] <- "spatial_name"
    final_table$spatial_name <- gsub("^CELL_", "SPATIAL_", final_table$spatial_name)
    colnames(gene_count) <- gsub("^CELL_", "SPATIAL_", colnames(gene_count))

    return(list(matched_data = final_table, gene_count = gene_count))
  }

  else {
    stop("Unknown technology_version in config file, please check again.")
  }

  # Compute pixel information if pixel = T
  if (pixel) {
    message("Start computing pixel information...")

    # set R_PACKAGE_DIR environment variable
    r_package_dir <- system.file(package = "stPipe")
    #Sys.setenv(R_PACKAGE_DIR = system.file(package = "stPipe"))
    pixel_output_csv <- file.path(qc_dir, "mapped_pixel.csv")

    # Call Python script
    basilisk::basiliskRun(
      env = stPipe:::stPipe_env,
      fun = function(test_img_fn, output_dir) {
        Visium_pixel <- reticulate::import_from_path("Visium_pixel", system.file("python", package = "stPipe"))
        Visium_pixel$main(test_img_fn, output_dir)
      },
      test_img_fn = image_path,
      output_dir = qc_dir
    )

    # Load pixel coordinates
    pixel_data <- utils::read.csv(pixel_output_csv, stringsAsFactors = FALSE)

    # Read ref_pos.csv internally
    ref_pos_path <- system.file("extdata", "ref_pos.csv", package = "stPipe")
    ref_pos <- utils::read.csv(ref_pos_path, header = FALSE, stringsAsFactors = FALSE)
    colnames(ref_pos) <- c("barcode_sequence", "in_tissue", "X_coordinate", "Y_coordinate", "x", "y")
    ref_pos <- ref_pos[-1, ]  # Remove the first row

    # Remove "-1" from barcode_sequence
    ref_pos$barcode_sequence <- gsub("-1$", "", ref_pos$barcode_sequence)
    # Add barcode_sequence to pixel_data
    pixel_data$barcode_sequence <- ref_pos$barcode_sequence
    # Merge pixel_data with matched_data
    matched_data <- merge(matched_data, pixel_data, by = "barcode_sequence", all.x = TRUE)
    colnames(matched_data)[7:8] <- c("Y_pixel", "X_pixel")
  }

  colnames(matched_data)[2] <- c("spatial_name")
  matched_data$spatial_name <- gsub("^CELL_", "SPATIAL_", matched_data$spatial_name)

  gene_count <- utils::read.csv(gene_count_file)
  rownames(gene_count) <- gene_count$gene_id
  gene_count$gene_id <- NULL
  colnames(gene_count) <- gsub("^CELL_", "SPATIAL_", colnames(gene_count))
  return(list(matched_data = matched_data, gene_count = gene_count))
}
