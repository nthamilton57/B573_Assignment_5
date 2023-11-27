#' @title Conducts disease enrichment analysis
#' @description Disease enrichment analysis is a valuable approach in genomics that helps identify associations between a set of genes and specific diseases or medical conditions. This function takes a dataframe of top genes, usually resulting from a differential expression analysis, and performs disease enrichment analysis on these genes.
#' @param top_gene_dataframe all results of differentially expressed genes
#'
#' @return the result of the disease enrichment analysis that includes information about diseases or medical conditions that are enriched among the provided top genes, along with statistical measures such as p-values and adjusted p-values.
#' @export
#'
#' @examples do_disease_enrichment_analysis(top_gene_dataframe)

do_disease_enrichment_analysis <- function(top_gene_dataframe){

  genelist <- as.vector(top_gene_dataframe$logFC)
  genelist <- sort(genelist, decreasing = TRUE)
  genelist <- setNames(genelist, top_gene_dataframe$Entrez_Gene_ID)

  x <- gseDO(genelist,
             minGSSize     = 120,
             pvalueCutoff  = 0.2,
             pAdjustMethod = "BH",
             verbose       = FALSE)
  x <- as.data.frame(x)

  return(x)

}
