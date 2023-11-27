#' @title Get annotated differentially expressed data
#' @description Get important information for genes that are differentially expressed
#' @param Annotation_Data Annotation Data
#' @param Fitted_File all results of differentially expressed genes
#' @param top_genes Number of top genes wanted to be contained in table
#'
#' @return dataframe with certain number of top genes
#' @export
#'
#' @examples make_annotations(Annotation_Data, Fitted_File, top_genes = 50)

make_annotations <- function(Object_File, Fitted_File, top_genes = 50){

  anno <- fData(Object_File)
  anno <- dplyr::select(anno, Symbol, Entrez_Gene_ID, Chromosome, Cytoband)  #selects the information we want to show for each gene
  Fitted_File$genes <- anno
  #top_gene_dataframe <- topTable(Fitted_File, number = top_genes)  #shows the topTable with the additional information

  #return(top_gene_dataframe)
  return(Fitted_File)

}
