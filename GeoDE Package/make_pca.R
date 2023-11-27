#' @title Create PCA Scatterplot
#' @description Create a scatter plot of your data utilizing Principle Component Analysis (PCA) vectors
#' @param Expression_Data Expression data obtained from GSE using exprs function
#' @param Phenotype_Data Phenotype data obtained from GSE using pData function
#'
#' @return plot scatterplot of the data using PCA
#' @export
#'
#' @examples make_pca(Expression_Data,Phenotype_Data)

make_pca <- function(Expression_Data,Phenotype_Data){

  pca <- prcomp(t(Expression_Data))
  ## Join the PCs to the sample information
  plot <- cbind(Phenotype_Data, pca$x) %>%
    ggplot(aes(x = PC1, y=PC2, col=group,label=paste("Patient", patient))) + geom_point() + geom_text_repel()
  return(plot)
}
