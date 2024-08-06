#' Normalization 0,1
#'
#' This function normalizes input data to the range (0,1).
#' @param mat Matrix with the first column containing measurement times and the second column
#' containing annotations (organs). The third and subsequent columns contain measurement data.
#' @param center The metric used to compute the center, options include "mean" and "median"
#' @return The data matrix with scaled measurements measurements. Scaled center measurements are provided for
#' each measurement/annotation (e.g., gene/organ) combination. This is the data format required for HMF().
#' @export
#' @examples datHMF = scale_zeroOne(matrix, center="mean")
scale_zeroOne <- function(mat,center="mean") {
  
  ###normalization 0,1
  hmfmelt <- melt(mat, id.vars = c("Organ", "Age")) %>% as.data.frame()
  
  hmfmean <- summarise_all(group_by(hmfmelt,Organ, Age, variable), list(function(x)(mean(x, na.rm = T))))
  hmfmean2 <- data.frame(Organ_gene = paste0(hmfmean$Organ,"-", hmfmean$variable), Age = hmfmean$Age, exp = hmfmean$value)
  
  hmfnorm <- apply(dcast(hmfmean2, Age~Organ_gene) %>% tibble::column_to_rownames("Age"),2, function(x)( x-min(x, na.rm = T) ) / ( max(x, na.rm = T) - min(x, na.rm = T) ))
  
  hmfnormmelt <- melt(hmfnorm, value.name = "Exp", varnames = c("Age","Organ-Gene")) %>% separate(`Organ-Gene`, into = c("Organ","Gene"), sep = "-")
  hmfnormmat <- dcast(hmfnormmelt, Age+Organ ~ Gene) %>% arrange(Organ); #hmfnormmat <- hmfnormmat[-which(apply(hmfnormmat[,3:ncol(hmfnormmat)], 1, function(x)(all(is.na(x))))),]
  
  return(hmfnormmat)
  
}
