#' @title Check Values and Plot a Box Plot
#' @description This function will take in expression data, create a boxplot of the data as well as return a summary of the data.
#'
#' @param Expression_Data Expression data obtained from GSE using exprs function
#'
#' @return A summary of the data and a boxplot of the data
#' @export
#' @examples
#' check_values(Expression_Data)

check_values <- function(Expression_Data){

  boxplot(Expression_Data)
  return(summary(Expression_Data))

}
