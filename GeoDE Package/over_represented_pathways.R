#' @title Use the identified DEGs to determine over-represented Reactome pathways
#' @description Assess whether the number of selected genes associated with a reactome pathway is larger than expected
#' @param top_gene_dataframe all results of differentially expressed genes
#' @param p_value_cutoff cutoff value for p-value
#' @param logfc_cutoff cutoff value for log fold change
#'
#' @return dataframe of over-represented reactome pathways
#' @export
#'
#' @examples over_represented_pathways(top_gene_dataframe, 0.05, 1.5)

over_represented_pathways <- function(top_gene_dataframe, p_value_cutoff, logfc_cutoff){

  de <- (top_gene_dataframe$Entrez_Gene_ID)[abs(top_gene_dataframe$logFC) > logfc_cutoff]

  x <- enrichPathway(gene=de, pvalueCutoff = p_value_cutoff, readable=TRUE)
  x <- as.data.frame(x)

  return(x)

}
