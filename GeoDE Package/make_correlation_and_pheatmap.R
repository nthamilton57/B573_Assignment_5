#' @title Correlation between Expression Data and Phenotypic Data
#' @description Does correlation between Expression Data and phenotypic data and displays the corresponding pheatmap.
#' @param Expression_Data Expression data obtained from GSE using exprs function
#' @param Phenotype_Data Phenotype data obtained from GSE using pData function
#'
#' @return correlation matrix and pheatmap
#' @export
#'
#' @examples make_correlation_and_pheatmap(Expression_Data, Phenotype_Data)
make_correlation_and_pheatmap <- function(Expression_Data, Phenotype_Data){

  corMatrix <- cor(Expression_Data,use="c")
  rownames(Phenotype_Data) <- colnames(corMatrix)
  out_graph <- pheatmap(corMatrix, annotation_col=Phenotype_Data)

  return(list('Correlation_Matrix' = corMatrix,
              'Graph' = out_graph))
}
