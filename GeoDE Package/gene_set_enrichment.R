#' @title Performs gene set enrichment analysis (GSEA) using the gsePathway function
#' @description The function returns the result of the gene set enrichment analysis. This result includes information about pathways or gene sets that are enriched among the provided top genes, along with statistical measures such as p-values and adjusted p-values.
#' @param top_gene_dataframe all results of differentially expressed genes
#'
#' @return information about pathways or gene sets that are enriched among the provided top genes
#' @export
#'
#' @examples gene_set_enrichment(top_genes)

gene_set_enrichment <- function(top_gene_dataframe){

  genelist <- as.vector(top_gene_dataframe$logFC)
  genelist <- sort(genelist, decreasing = TRUE)
  genelist <- setNames(genelist, top_gene_dataframe$Entrez_Gene_ID)

  y <- gsePathway(genelist,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  verbose = FALSE)

  return(y)

}
