#' @title Get phenotypic data with important factors
#' @description Get factors that may be needed for analysis.
#' @param Phenotype_Data Phenotype data obtained from GSE using pData function
#'
#' @return Phenotypic data with simplified factors that will be useful for differential expression analysis
#' @export
#'
#' @examples make_sample_simple(Phenotype_Data)
make_sample_simple <- function(Phenotype_Data){

  # Phenotype_Data <- select(sampleInfo, source_name_ch1,characteristics_ch1.1)
  # Phenotype_Data <- rename(sampleInfo,group = source_name_ch1, patient=characteristics_ch1.1)

  library(dplyr)

  # Use select and rename from dplyr to transform the data
  Phenotype_Data <- Phenotype_Data %>%
    dplyr::select(source_name_ch1, characteristics_ch1.1) %>%
    dplyr::rename(group = source_name_ch1, patient = characteristics_ch1.1)

  return(Phenotype_Data)
}
