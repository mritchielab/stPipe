#' stPipe: A package for pre-processing sequencing-based spatial transcriptomics data.
#'
#' @description The stPipe will do spatial barcode demultiplexing, UMI deduplication, spatial location matching and quality control on
#' fastq data generated from all mainstream protocols
#'
#' @section stPipe functions:
#' The stPipe functions Run_ST, Run_loc_match, Run_QC, Run_Vis, Run_Clustering, Run_create_obj, Run_HTML, Run_Interactive
#'
#'
#' @author Yang Xu <xu.ya@wehi.edu.au>
#' @docType package
#' @name stPipe
#' @import Rhtslib
#' @import SingleCellExperiment
#' @import SpatialExperiment
#' @import Rcpp
#' @importFrom methods as
#' @useDynLib stPipe, .registration = TRUE
#' @aliases stPipe stPipe-package
NULL
