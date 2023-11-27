#' @title Importing GEO object and files
#'
#' @param GEO_ID GEO ID for the dataset
#'
#' @return A list of phenotype data, expression data, annotation data and object file
#' @export
#'
#' @examples GEO_ID()
import_GEO <- function(GEO_ID){

  gse <- getGEO(GEO_ID)
  gse <- gse[[1]]

  subj <- pData(gse)
  expr <- exprs(gse)
  annot <- fData(gse)

  out_list <- list('Phenotype_Data' = subj,
                   'Expression_Data' = expr,
                   'Annotation_Data' = annot,
                   'Object_file' = gse,
                   'gse' = gse)

  return(out_list)

}
