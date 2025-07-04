#' Paired-End FASTQ Demo Files
#'
#' A pair of gzipped FASTQ files for demo (FFPE, probe-based mouse spleen sample).
#'
#' @docType data
#' @name demo_data
#' @aliases demo_R1.fq.gz demo_R2.fq.gz
#' @keywords datasets
#'
#' @usage demo_data
#'
#' @format A character vector of length 2:
#' \describe{
#'   \item{demo\_R1.fq.gz}{Read1 FASTQ file}
#'   \item{demo\_R2.fq.gz}{Read2 FASTQ file}
#' }
#'
#' @source
#' <https://zenodo.org/records/12788291>
NULL

demo_data <- function() {
  system.file(
    "extdata",
    c("demo_R1.fq.gz", "demo_R2.fq.gz"),
    package = "stPipe"
  )
}


#' Visium Spatial Barcode Coordinates (Versions 1–5)
#'
#' Five separate datasets containing spot barcodes and their
#' X/Y coordinates on the Visium capture area. Each object is
#' a data.frame with three columns which stand for: `barcode`, `x`, and `y`.
#'
#' @name visium_v1_coordinates
#' @aliases visium_v1_coordinates
#' @docType data
#' @title Visium v1 Spatial Barcode Coordinates
#'
#' @format A data.frame with 4992 rows (14336 rows for v5) and 3 columns:
#' \describe{
#'   \item{barcode}{Spot barcode (e.g., "AAACGCTGAAAGCGTC-1")}
#'   \item{x}{X-coordinate on the capture area}
#'   \item{y}{Y-coordinate on the capture area}
#' }
#'
#' @usage
#' visium_v1_coordinates
#'
#'
#' @source
#' 10x Genomics Support: Visium Spatial Barcode Inclusion List
#' (<https://kb.10xgenomics.com/hc/en-us/articles/360041426992-Where-can-I-find-the-Space-Ranger-barcode-inclusion-list-formerly-barcode-whitelist-and-their-coordinates-on-the-slide>)
NULL

visium_coordinates <- function(version = 1L) {
  ver <- match.arg(as.character(version), as.character(1:5))
  fname <- sprintf("visium-v%s_coordinates.rds", ver)
  path <- system.file("extdata", fname, package = "stPipe")
  if (nzchar(path)) {
    readRDS(path)
  } else {
    stop("Unsupported Visium version: ", version)
  }
}

visium_v1_coordinates <- function() {
  system.file("extdata", "visium_v1_coordinates.rds", package = "stPipe")
}

#' @rdname visium_v1_coordinates
#' @aliases visium_v2_coordinates
#' @title Visium v2 Spatial Barcode Coordinates
#'
#' @name visium_v2_coordinates
#' @docType data
NULL

#' @rdname visium_v1_coordinates
#' @aliases visium_v3_coordinates
#' @title Visium v3 Spatial Barcode Coordinates
#'
#' @name visium_v3_coordinates
#' @docType data
NULL

#' @rdname visium_v1_coordinates
#' @aliases visium_v4_coordinates
#' @title Visium v4 Spatial Barcode Coordinates
#'
#' @name visium_v4_coordinates
#' @docType data
NULL

#' @rdname visium_v1_coordinates
#' @aliases visium_v5_coordinates
#' @title Visium v5 Spatial Barcode Coordinates
#'
#' @name visium_v5_coordinates
#' @docType data
NULL

#' Probe set for Human Transcriptome v1
#'
#' @docType data
#' @name Visium_Human_Transcriptome_Probe_Set_v1_0_GRCh38_2020_A
#' @aliases Visium_Human_Transcriptome_Probe_Set_v1_0_GRCh38_2020_A.rds
#' @keywords datasets
#'
#' @usage Visium_Human_Transcriptome_Probe_Set_v1_0_GRCh38_2020_A
#'
#' @format A data.frame with 18681 rows and 2 columns:
#' \describe{
#'   \item{gene_id}{Ensembl gene ID}
#'   \item{probe_seq}{Sequence for designed probe set}
#' }
#'
#' @source
#' <https://www.10xgenomics.com/support/software/space-ranger/downloads#probe-set-downloads>
NULL

Visium_Human_Transcriptome_Probe_Set_v1_0_GRCh38_2020_A <- function() {
  system.file("extdata", "Visium_Human_Transcriptome_Probe_Set_v1_0_GRCh38_2020_A.rds", package = "stPipe")
}

#' Probe set for Human Transcriptome v2
#'
#' @docType data
#' @name Visium_Human_Transcriptome_Probe_Set_v2_0_GRCh38_2020_A
#' @aliases Visium_Human_Transcriptome_Probe_Set_v2_0_GRCh38_2020_A.rds
#' @keywords datasets
#'
#' @usage Visium_Human_Transcriptome_Probe_Set_v2_0_GRCh38_2020_A
#'
#' @format A data.frame with 53504 rows and 2 columns:
#' \describe{
#'   \item{gene_id}{Ensembl gene ID}
#'   \item{probe_seq}{Sequence for designed probe set}
#' }
#'
#' @source
#' <https://www.10xgenomics.com/support/software/space-ranger/downloads#probe-set-downloads>
NULL

Visium_Human_Transcriptome_Probe_Set_v2_0_GRCh38_2020_A <- function() {
  system.file("extdata", "Visium_Human_Transcriptome_Probe_Set_v2_0_GRCh38_2020_A.rds", package = "stPipe")
}

#' Probe set for Mouse Transcriptome v1
#'
#' @docType data
#' @name Visium_Mouse_Transcriptome_Probe_Set_v1_0_mm10_2020_A
#' @aliases Visium_Mouse_Transcriptome_Probe_Set_v1_0_mm10_2020_A.rds
#' @keywords datasets
#'
#' @usage Visium_Mouse_Transcriptome_Probe_Set_v1_0_mm10_2020_A
#'
#' @format A data.frame with 20057 rows and 2 columns:
#' \describe{
#'   \item{gene_id}{Ensembl gene ID}
#'   \item{probe_seq}{Sequence for designed probe set}
#' }
#'
#' @source
#' <https://www.10xgenomics.com/support/software/space-ranger/downloads#probe-set-downloads>
NULL

Visium_Mouse_Transcriptome_Probe_Set_v1_0_mm10_2020_A <- function() {
  system.file("extdata", "Visium_Mouse_Transcriptome_Probe_Set_v1_0_mm10_2020_A.rds", package = "stPipe")
}

#' H&E Image for mouse olfactory bulb
#'
#' @docType data
#' @name tissue_hires_image
#' @aliases tissue_hires_image.png
#' @keywords datasets
#'
#' @usage tissue_hires_image
#'
#' @format A high-resolution image file
#'
#' @source
#' <https://github.com/ucagenomix/SiT>
NULL

tissue_hires_image <- function() {
  system.file("extdata", "tissue_hires_image.png", package = "stPipe")
}


#' Coordinates for mouse olfactory bulb H&E Image
#'
#' @docType data
#' @name ref_pos
#' @aliases ref_pos.csv
#' @keywords datasets
#'
#' @usage ref_pos
#'
#' @format A data.frame with 4992 rows and 5 columns:
#' \describe{
#'   \item{barcode}{Spot barcode (e.g., "AAACGCTGAAAGCGTC-1")}
#'   \item{x}{X-coordinate on the capture area}
#'   \item{y}{Y-coordinate on the capture area}
#'   \item{x-pixel}{Y-pixel on the capture area}
#'   \item{y-pixel}{Y-pixel on the capture area}
#'}
#'
#' @source
#' <https://github.com/ucagenomix/SiT>
NULL

ref_pos <- function() {
  path <- system.file("extdata", "ref_pos.csv", package = "stPipe")
  read.csv(path, stringsAsFactors = FALSE)
}
