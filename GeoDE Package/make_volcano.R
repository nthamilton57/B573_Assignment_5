#' @title Creating a volcano plot
#' @description The function produces a volcano plot, a widely-used visualization in genomics. In this plot, each point represents a gene, with the x-axis indicating the log-fold change and the y-axis representing the negative logarithm of the adjusted p-value. Points are color-coded to distinguish between significant and non-significant results.
#' @param Fitted_File all results of differentially expressed genes
#' @param p_cutoff the p-value cutoff for the graph
#' @param fc_cutoff the foldchange cutoff for the graph
#'
#' @return a volcano graph
#' @export
#'
#' @examples make_volcano(Fitted_File, p_cutoff = 0.05, fc_cutoff = 1)

make_volcano <- function(Fitted_File, p_cutoff = 0.05, fc_cutoff = 1){

  full_results <- topTable(Fitted_File, number=Inf)    #include all results
  full_results <- tibble::rownames_to_column(full_results,"ID")


  ## change according to your needs
  p_cutoff <- p_cutoff
  fc_cutoff <- fc_cutoff

  out_graph <- full_results %>%
    mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>%
    ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()

  out <- list(out_graph, full_results)
  return(out)
}
