#' @title Log Normalize the Expression Data
#' @description This function takes expression data and returns a the data with a log base 2 transformation applied.
#'
#' @param Expression_Data Expression data obtained from GSE using exprs function
#'
#' @return Normalized Expression Data
#' @export
#'
#' @examples log_normalize(Expression_Data)
log_normalize <- function(Expression_Data){

  Normal_Expression_Data <- log2(Expression_Data)    #transform the data to smaller values/more manageable
  return(Normal_Expression_Data)

}
