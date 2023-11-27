#' @title Visualize Reactome pathway
#' @description implemented the viewPathway() to visualized selected reactome pathway
#' @param pathway the selected reactome pathway to visualize
#'
#' @return nothing
#' @export
#'
#' @examples visualize_pathway_network("Nitric oxide stimulates guanylate cyclase")

visualize_pathway_network <- function(pathway){

  viewPathway(pathway,
              readable = TRUE)


}
