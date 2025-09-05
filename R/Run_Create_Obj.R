#####################################################
# Generate downstream objects for mainstream tools for further analysis
#####################################################

#' @name Run_Create_Obj
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
#' @examples
#' set.seed(123)
#' gene.count <- matrix(sample(0:100, 200*100, TRUE), nrow=200)
#' rownames(gene.count) <- paste0('Gene', 1:200)
#' colnames(gene.count) <- paste0('Spot', 1:100)
#' matched.data <- data.frame(
#'   spatial_name = colnames(gene.count),
#'   barcode_sequence = paste0('BC', 1:100),
#'   X_coordinate = runif(100,0,10),
#'   Y_coordinate = runif(100,0,10),
#'   UMI_count = sample(50:200,100,TRUE)
#' )
#' seu <- Run_Create_Obj(
#'   gene.matrix = gene.count,
#'   matched.data = matched.data,
#'   obj.type = 'Seurat',
#'   tech = 'Slideseq'
#' )
#' @export
#' @importFrom ggplot2 theme
#' @importFrom methods new
#' @importFrom SummarizedExperiment rowData<-

Run_Create_Obj <- function(gene.matrix, matched.data, obj.type, tech, ss.radius = 3000) {
  
  if (!all(matched.data$spatial_name %in% colnames(gene.matrix))) {
    stop("Spatial locations in matched.data do not match column names of gene.matrix!")
  }
  
  counts_mat <- as.matrix(gene.matrix)
  matched.data <- matched.data[match(colnames(counts_mat), matched.data$spatial_name), ]
  meta_df <- data.frame(
    X   = matched.data$X_coordinate,
    Y   = matched.data$Y_coordinate,
    UMI = matched.data$UMI_count,
    row.names = matched.data$barcode_sequence
  )
  colnames(counts_mat) <- matched.data$barcode_sequence
  
  # Create Seurat object
  if (obj.type == "Seurat") {
    
    # Create a Seurat object
    obj <- SeuratObject::CreateSeuratObject(counts = counts_mat, meta.data = meta_df, assay = 'Spatial')
    coords <- meta_df[, c("X","Y"), drop = FALSE]
    
    obj[['image']] <- new(
      Class = 'SlideSeq',
      assay = "Spatial",
      coordinates = coords)
    
    # add spatial information to created Seurat object based on different technology
    if (tech == "Slideseq" | tech == "Curio-seeker"){
      
      # remove stray beads which fall outside the main Slideseq puck area
      obj <- Seurat::FilterSlideSeq(object = obj, radius = ss.radius, do.plot = TRUE)
      
    }
  } else if (obj.type == "SpatialExperiment") {
    
    obj <- SpatialExperiment::SpatialExperiment(
      assays = list(counts = counts_mat),
      colData = meta_df,
      spatialCoords = as.matrix(meta_df[, c('X','Y')])
    )
    rowData(obj)$gene_name <- rownames(gene.matrix)
    
    
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
    
    stop("Unknown obj.type! Please choose class from 'Seurat', 'SpatialExperiment', or 'AnnData'.")
  }
  
  return(obj)
  
}
