#####################################################
# Generate downstream objects for mainstream tools for further analysis
#####################################################

#' @name Run_create_obj
#' @title This function creates specified spatial transcriptomics data object for further personalized downstream analysis
#' @details Current mainstream analytic tools for downstream includes: Seurat, SpatialExperiment supported tools, and Squidpy.
#' This function can help creates corresponding objects for further downstream analysis.
#'
#' @param gene.matrix Gene count matrix. This is usually obtained from 'Run_ST' function.
#' @param matched.data Data frame which contains spatial localization matched with gene count matrix. This is usually obtained from 'Run_loc_match' function.
#' @param obj.type Type of spatial transcriptomics object being created, can be 'Seurat', 'SpatialExperiment', or 'AnnData'.
#' @param ss.radius Optional. Radius for filtering SlideSeq or Curio-seeker spots for Seurat object. Default is 3000.
#' @param tech Type of spatial transcriptomics sequencing technology, can be "Visium", "Slideseq", "Curio-seeker", or "Stereoseq".
#' @return Created spatial transcriptomics data object as required by 'obj.type'.
#' #' @examples
#' \dontrun{
#' my_ST_Obj <- Run_create_obj(gene.matrix = gene_count, matched.data = matched_spatial_barcode, obj.type = "SpatialExperiment", tech = "Visium")
#' }
#' @export
#' @importFrom ggplot2 theme

Run_create_obj <- function(gene.matrix, matched.data, obj.type, tech, ss.radius = 3000) {

  if (!all(matched.data$cell_name %in% colnames(gene_count))) {
    stop("Spatial locations in matched.data do not match column names of gene.matrix!")
  }

  # replace spot names with spatial barcode sequences
  gene_columns <- colnames(gene.matrix)
  mapping <- setNames(matched.data$barcode_sequence, matched.data$cell_name)
  new_columns <- mapping[gene_columns]
  if (all(!is.na(new_columns))) {
    colnames(gene.matrix) <- new_columns
  } else {
    warning("Some spatial location names did not have a matching barcode sequence.")
  }

  # get matched meta.data with gene.matrix and reformat
  matched.data <- obj@meta.data[,4:5]
  rownames(matched.data) <- matched.data$barcode_sequence
  matched.data <- matched.data[, 4:5]
  colnames(matched.data) <- c("xcoord", "ycoord")

  # Create Seurat object
  if (obj.type == "Seurat") {

    # Create a Seurat object
    obj <- SeuratObject::CreateSeuratObject(counts = gene.matrix, meta.data = matched.data, assay = 'Spatial')

    # add spatial information to created Seurat object based on different technology
    if (tech == "Visium"){

      matched.data$ycoord <- -matched.data$ycoord
      assay_name <- "Spatial"
      fov_name <- "fov"
      fov <- SeuratObject::CreateFOV(
        matched.data,
        type = "centroids",
        assay = assay_name,
        key = Key(fov_name, quiet = TRUE))
      fov <- fov[SeuratObject::Cells(obj)]
      obj[[fov_name]] <- fov

    } else if (tech == "Slideseq" | tech == "Curio-seeker") {

    obj[['image']] <- new(
      Class = 'SlideSeq',
      assay = "Spatial",
      coordinates = matched.data
    )
    # remove stray beads which fall outside the main Slideseq puck area
    obj <- Seurat::FilterSlideSeq(object = obj, radius = ss.radius, do.plot = T)
    }


  } else if (obj.type == "SpatialExperiment") {

    obj <- SpatialExperiment::SpatialExperiment(
      assays = list(counts = gene.matrix),
      spatialCoords = as.matrix(matched.data)
    )
    rowData(obj)$gene_name <- rownames(rowData(obj))


  } else if (obj.type == "AnnData") {

    # Create an AnnData object
    anndata <- reticulate::import("anndata")
    obj <- anndata$AnnData(X = gene.matrix)
    obj$obs <- matched.data

    py_create_anndata <- "
    import anndata
    adata = anndata.AnnData(X = r.gene_matrix)
    adata.obs = r.matched_data
    "
    reticulate::py_run_string(py_create_anndata)

    obj <- reticulate::py$adata

  } else {

    stop("Unsupported obj.type! Please choose from 'Seurat', 'SpatialExperiment', or 'AnnData'.")
  }

  return(obj)



}
