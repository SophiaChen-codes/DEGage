#Alicia Petrany, draft as of 4.25.2023
#This file contains functions for DEGage_Preprocess

#' init.prepprocess.dependencies
#'
#' initializes dependencies for DEGage's preprocess module
#' @return Nothing
init.preprocess.dependencies <- function(){
  if(!suppressPackageStartupMessages(require("SingleR"))){
    suppressPackageStartupMessages(install.packages("SingleR"))
    library("SingleR")
  }
  if(!suppressPackageStartupMessages(require("SingleCellExperiment"))){
    suppressPackageStartupMessages(install.packages("SingleCellExperiment"))
    library("SingleCellExperiment")
  }
  if(!suppressPackageStartupMessages(require("Seurat"))){
    suppressPackageStartupMessages(install.packages("Seurat"))
    library("Seurat")
  }
}

#'preprocess.QC
#'
#'performs QC on user input
preprocess.QC <- function(input, dir.type, output){
  if(class(input) == "character"){
    if(!dir.exists(input)){
      stop("The provided input directory does not exist")
    }
  }

  allowed.dir.types <- c("mtx")
  if(!(dir.type %in% allowed.dir.types)){
    stop("dir.type is not an accepted format. See ?DEGage_preprocess for supported data types")
  }

  allowed.output.types <- c("Seurat", "df")
  if(!(output %in% allowed.output.types)){
    stop("output is not an accepted format. See ?DEGage_preprocess for supported formats")
  }

}


#'parse.input
#'
#'parses user input into a seurat object
#'@param input User input of whatever format
#'@param type type of directory that a user input path leads to, if relevant.
parse.input <- function(input, type){
  dtype <- class(input)
  if(dtype == "data.frame"){
    return(CreateSeuratObject(input))
  }
  else if(dtype == "Seurat"){
    return(input)
  }
  else if(dtype == "character" && type == "mtx"){
    output <- Read10X(data.dir = input)
    return(CreateSeuratObject(output))
  }
  else{
    stop("Input is not of a supported data type. See ?DEGage_Preprocess for supported data types.")
  }
}

#'filter.and.normalize
#'
#'fiters and normalizes counts according to lightweight Seurat workflow
#'@param sobj A Seurat object containing raw counts
#'@returns A seurat object containing raw counts
filter.and.normalize <- function(sobj, min.nFeatureRNA, max.nFeatureRNA, mt.percent){
  #TODO: add parameters for min and max nFeature_RNA
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  sobj <- subset(sobj, subset = nFeature_RNA > min.nFeatureRNA &
                   nFeature_RNA < max.nFeatureRNA &
                   percent.mt < (mt.percent*100))
  sobj <- NormalizeData(sobj)
  sobj <- FindVariableFeatures(object = sobj)
  return(sobj)
}

#'celltype.annotation
#'
#'performs celltype annotations with single R, if possible. Can be limited if dataset contains a small number of cells
#'@param sobj A seurat object containing normalized counts
#'@param celltype.min Filters out cells of celltypes that are present in too few cells.
#'
celltype.annotation <- function(sobj, celltype.min, output){

  #generating singleR predictions
  sce <- as.SingleCellExperiment(sobj)
  hpca.sce <- suppressMessages(celldex::HumanPrimaryCellAtlasData())

  #its a trycatch because singleR is finicky and honestly I have no idea why
  preds = tryCatch({preds<- SingleR(test = sce, ref = hpca.sce,
                                      labels = hpca.sce$label.main)},
                error = function(e){return(NULL)})
  if(is.null(preds)){
    stop("Could not generate cell type annotations with SingleR. It's possible that there are too few cells in your dataset. Try DEGage() instead.")
  }
  #Filtering out junk cells
  pred.table <- table(preds$labels)
  cell.list <- levels(factor(preds$labels))
  too.few <- vector()
  for(i in 1:length(pred.table)){
    if (pred.table[[i]] < celltype.min){
      too.few <- c(too.few, cell.list[i])
    }
  }

  i <- which(preds$labels %in% too.few  )
  sobj$celltype <- preds$labels
  sobj <- sobj[,-i]

  sobj <- clustering.w.singleR(sobj, output)
  return(sobj)
}

#'cluster.w.singleR
#'
#'generates tsne coordinates for cells, then outputs the plot with cell types annotated
#'
#'@param sobj A seurat object containing processed counts
#'@param output The type of output the user specified
#'@return nothing, calls output
clustering.w.singleR <- function(sobj, output){
  sobj <- ScaleData(sobj)
  sobj <- RunPCA(sobj)
  sobj <- RunTSNE(sobj)
  return(output.format.w.celltypes(sobj, output))
}

#'output
#'
#'formats output to the users preference
#'@param sobj Seurat objecct containing processed counts and cell type annotations
#'@param output format of output as specific by user
#'@return filtered counts and cell type annotations as specified by user
output.format.w.celltypes <- function(sobj, output){
  if(output == "Seurat"){
    return(sobj)
  }
  if(output == "df"){
    return(list(sobj@assays$RNA@counts,
                 sobj$celltype))
  }
}


