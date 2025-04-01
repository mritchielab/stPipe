#####################################################
# Generate HTML Report for upstream results
#####################################################

#' @name Run_HTML
#' @title HTML report generation in RMarkDown format
#'
#' @description This function generates a HTML report of stPipe upstream processing steps in Rmarkdown format. It extracts plots in the result directory and outputs them in the whole report with corresponding explanation.
#' @examples
#' ## @example
#' temp_dir <- tempdir()
#' dummy_gene_count <- data.frame(GeneA = c(10, 20), GeneB = c(5, 15))
#' rownames(dummy_gene_count) <- c("spot1", "spot2")
#' write.csv(dummy_gene_count, file = file.path(temp_dir, "gene_count.csv"), row.names = TRUE)
#' png_files <- c("Mapping_statistics_plot.png", "UMI_duplication_plot.png", 
#'              "Barcode_demultiplexing_plot.png", "Spatial_Heatmap_of_UMI_Count_raw.png", 
#'                          "Spatial_Heatmap_of_UMI_Count_log.png", "Threshold.png")
#'                          for (f in png_files) {
#'                            fileConn <- file(file.path(temp_dir, f))
#'                              writeLines("dummy image content", fileConn)
#'                                close(fileConn)}
#' html_files <- c("tsne_interactive.html", "umap_interactive.html")
#' for (f in html_files) {
#'  fileConn <- file(file.path(temp_dir, f))
#'  writeLines("<html>dummy interactive html</html>", fileConn)
#'  close(fileConn)}
#' dummy_rmd <- paste(
#'  "Report for stPipe",
#'  "Mapping Statistics: <<Mapping_statistics_plot>>",
#'  "UMI Duplication: <<UMI_duplication_plot>>",
#'  "Barcode Demultiplexing: <<Barcode_demultiplexing_plot>>",
#'  "Spatial Heatmap Raw: <<Spatial_Heatmap_of_UMI_Count_raw>>",
#'  "Spatial Heatmap Log: <<Spatial_Heatmap_of_UMI_Count_log>>",
#'  "Threshold: <<Threshold>>",
#'  "tSNE Interactive: <<tsne_interactive>>",
#'  "UMAP Interactive: <<umap_interactive>>",
#'  sep = "\n"
#')
#' dummy_template_path <- file.path(temp_dir, "stPipe_Report_Skeleton.Rmd")
#' writeLines(dummy_rmd, dummy_template_path)
#' Run_HTML_test <- function(path) {
#'  rmd_template_path <- dummy_template_path
#'  if (rmd_template_path == "") {
#'    stop("Template file not found. Please check if it's properly installed in your package.")
#'  }
#' rmd_template <- readLines(rmd_template_path)
#'  gene_count <- read.csv(file.path(path, "gene_count.csv"), row.names = 1)
#'    Mapping_statistics_plot <- normalizePath(file.path(path, "Mapping_statistics_plot.png"))
#'  UMI_duplication_plot <- normalizePath(file.path(path, "UMI_duplication_plot.png"))
#'  Barcode_demultiplexing_plot <- normalizePath(file.path(path, "Barcode_demultiplexing_plot.png"))
#'  Spatial_Heatmap_of_UMI_Count_raw <- normalizePath(file.path(path, "Spatial_Heatmap_of_UMI_Count_raw.png"))
#'  Spatial_Heatmap_of_UMI_Count_log <- normalizePath(file.path(path, "Spatial_Heatmap_of_UMI_Count_log.png"))
#'  Threshold <- normalizePath(file.path(path, "Threshold.png"))
#'  tsne_interactive <- normalizePath(file.path(path, "tsne_interactive.html"))
#'  umap_interactive <- normalizePath(file.path(path, "umap_interactive.html"))
#'  rmd_content <- gsub("<<path>>", path, rmd_template)
#'  rmd_content <- gsub("<<Mapping_statistics_plot>>", Mapping_statistics_plot, rmd_content)
#'  rmd_content <- gsub("<<UMI_duplication_plot>>", UMI_duplication_plot, rmd_content)
#'  rmd_content <- gsub("<<Barcode_demultiplexing_plot>>", Barcode_demultiplexing_plot, rmd_content)
#'  rmd_content <- gsub("<<Spatial_Heatmap_of_UMI_Count_raw>>", Spatial_Heatmap_of_UMI_Count_raw, rmd_content)
#'  rmd_content <- gsub("<<Spatial_Heatmap_of_UMI_Count_log>>", Spatial_Heatmap_of_UMI_Count_log, rmd_content)
#'  rmd_content <- gsub("<<Threshold>>", Threshold, rmd_content)
#'  rmd_content <- gsub("<<tsne_interactive>>", tsne_interactive, rmd_content)
#'  rmd_content <- gsub("<<umap_interactive>>", umap_interactive, rmd_content)
#'  rmd_file <- file.path(path, "report.Rmd")
#'  writeLines(rmd_content, rmd_file)
#'  rmarkdown::render(rmd_file, output_file = file.path(path, "report.html"))}
#' Run_HTML_test(temp_dir)
#' @param path Path to where results are stored, including the ones obtained from 'Run_ST', 'Run_loc_match', 'Run_Vis', and 'Run_Clustering'.
#' @return Generated HTML report in RMarkDown format.
#' @importFrom rmarkdown render
#' @export

Run_HTML <- function(path) {

  rmd_template_path <- system.file("rmarkdown/templates/stpipe-report/skeleton/stPipe_Report_Skeleton.Rmd", package = "stPipe")
  if (rmd_template_path == "") {
    stop("Template file not found. Please check if it's properly installed in your package.")
  }

  rmd_template <- readLines(rmd_template_path)

  # gene count information
  gene_count <- read.csv(file.path(path, "gene_count.csv"), row.names = 1)
  #filtered_gene_count <- read.csv(file.path(path, "upstream_QC/filtered_gene_count.csv"), row.names = 1)

  # barplot for QC filtering
  #total_gene_count <- length(colnames(gene_count))
  #filtered_gene_count_len <- length(colnames(filtered_gene_count))
  #filtered_out <- total_gene_count - filtered_gene_count_len

  # obtain plot path
  #bar_plot_path <- normalizePath(file.path(path, "gene_count_diff_plot.png"))
  Mapping_statistics_plot <- normalizePath(file.path(path, "Mapping_statistics_plot.png"))
  UMI_duplication_plot <- normalizePath(file.path(path, "UMI_duplication_plot.png"))
  Barcode_demultiplexing_plot <- normalizePath(file.path(path, "Barcode_demultiplexing_plot.png"))
  Spatial_Heatmap_of_UMI_Count_raw <- normalizePath(file.path(path, "Spatial_Heatmap_of_UMI_Count_raw.png"))
  Spatial_Heatmap_of_UMI_Count_log <- normalizePath(file.path(path, "Spatial_Heatmap_of_UMI_Count_log.png"))
  Threshold <- normalizePath(file.path(path, "Threshold.png"))
  tsne_interactive <- normalizePath(file.path(path, "tsne_interactive.html"))
  umap_interactive <- normalizePath(file.path(path, "umap_interactive.html"))

  #png(bar_plot_path)
  #barplot(c(total_gene_count, filtered_gene_count_len),
  #        names.arg = c("Original", "Filtered"),
  #        main = "Number of Spatial location Before and After QC",
  #        col = c("red", "green"))
  #dev.off()
  #description_text <- sprintf("After the quality control step, %d spots/bins/beads have been filtered out due to their low quality and %d are remained.",
  #                            filtered_out, filtered_gene_count_len)

  # replace placeholder
  #rmd_content <- gsub("<<bar_plot_path>>", bar_plot_path, rmd_content)
  #rmd_content <- gsub("<<filtered_out>>", filtered_out, rmd_content)
  #rmd_content <- gsub("<<filtered_gene_count_len>>", filtered_gene_count_len, rmd_content)
  #rmd_content <- gsub("description_text", description_text, rmd_content
  rmd_content <- gsub("<<path>>", path, rmd_template)
  rmd_content <- gsub("<<Mapping_statistics_plot>>", Mapping_statistics_plot, rmd_content)
  rmd_content <- gsub("<<UMI_duplication_plot>>", UMI_duplication_plot, rmd_content)
  rmd_content <- gsub("<<Barcode_demultiplexing_plot>>", Barcode_demultiplexing_plot, rmd_content)
  rmd_content <- gsub("<<Spatial_Heatmap_of_UMI_Count_raw>>", Spatial_Heatmap_of_UMI_Count_raw, rmd_content)
  rmd_content <- gsub("<<Spatial_Heatmap_of_UMI_Count_log>>", Spatial_Heatmap_of_UMI_Count_log, rmd_content)
  rmd_content <- gsub("<<Threshold>>", Threshold, rmd_content)
  rmd_content <- gsub("<<tsne_interactive>>", tsne_interactive, rmd_content)
  rmd_content <- gsub("<<umap_interactive>>", umap_interactive, rmd_content)

  # write into new RMarkdown file
  rmd_file <- file.path(path, "report.Rmd")
  writeLines(rmd_content, rmd_file)

  # generate HTML report under specified path
  rmarkdown::render(rmd_file, output_file = file.path(path, "report.html"))
}
