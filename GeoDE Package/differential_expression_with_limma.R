#' @title Perform Differential Expression with Limma
#' @description This function is designed for conducting differential expression analysis using the Limma package in R. The function returns a Limma object containing the results of the differential expression analysis. The Limma object encapsulates information such as log-fold changes, p-values, and adjusted p-values for each gene.
#'
#' @param Phenotype_Data Phenotype data obtained from GSE using pData function
#' @param Expression_Data Expression data obtained from GSE using exprs function
#' @param Object_File GSE object for the current analysis
#' @param threshold_genes Minimum number of subjects to have these genes
#'
#' @return Returns a fitted object with logFC, pvalue and un-annotated names
#' @export
#'
#' @examples
#' differential_expression_with_limma(Phenotype_Data, Expression_Data, Object_File, threshold_genes = 2)


# differential_expression_with_limma <- function(Phenotype_Data, Expression_Data, Object_File, threshold_genes = 2){
#
#   design <- model.matrix(~0+Phenotype_Data$group)
#   colnames(design) <- c("Normal","Tumour")
#   cutoff <- median(Expression_Data)
#   is_expressed <- (Expression_Data) > cutoff
#   keep <- rowSums(is_expressed) > threshold_genes
#   Object_File <- Object_File[keep,]
#
#   fit <- lmFit(exprs(Object_File), design)
#   contrasts <- makeContrasts(Tumour - Normal, levels=design)
#   fit2 <- contrasts.fit(fit, contrasts)
#   fit2 <- eBayes(fit2)
#
#   return(fit2)
#
# }


differential_expression_with_limma <- function(sampleInfo,gse){

  design <- model.matrix(~0+sampleInfo$group)   #the ~0+sampleInfo$group is like the equivalent of the linear regression formula where 0 is the intercept for each subject (baseline for each subject)

  ## the column names are a bit ugly, so we will rename
  colnames(design) <- c("Normal","Tumour")

  ## calculate median expression level
  cutoff <- median(exprs(gse))

  ## TRUE or FALSE for whether each gene is "expressed" in each sample. this will get rid of genes that are expressed lower than the median
  is_expressed <- exprs(gse) > cutoff

  ## Identify genes expressed in more than 2 samples
  keep <- rowSums(is_expressed) > 2

  ## subset to just those expressed genes
  gse <- gse[keep,]

  #here we are fitting the model to see if the group affects the gene expression
  fit <- lmFit(exprs(gse), design)

  #n order to perform the differential analysis, we have to define the contrast that we are interested in. In our case we only have two groups and one contrast of interest. Multiple contrasts can be defined in the makeContrasts function.
  #Here we are looking at differences between Tumour and Normal
  contrasts <- makeContrasts(Tumour - Normal, levels=design)    #shows that -1 is Normal and 1 is the Tumour cells

  fit2 <- contrasts.fit(fit, contrasts) #runs the differential expression with the contrasts now
  fit2 <- eBayes(fit2)

  return(fit2)

}
