#' Select M ranges
#'
#' Select a set of values for M = max(m)
#' @param datHMF The data matrix with scaled measurements measurements. Scaled center measurements are provided for
#' each measurement/annotation (e.g., gene/organ) combination. This format is produced by scale_zeroOne(). Defaults to datHMF.
#' @param Mall A set of m-values. Defaults to NULL.
#' @param Msc A scale factor that determines the M value (M = max(m)) corresponding to the largest range of m-values.
#' Mmax <- ceiling(Msc * Np / 2) where Np = number of parameters.
#' Defaults to 2.
#' @return A List with the minimum number of modulating functions (Mmin), the maximum number of modulating functions (Mmax),
#' and the entire set of modulating function m-values (Mall).
#' @export
#' @examples mvals = get_M_vals(datHMF, Msc=2)
##### select M ranges
get_M_vals <- function(datHMF=datHMF,Mall=NULL,Msc=2){
    if(Msc != 2){Msc = Msc}   # scale factor, user defined input
    Np <- (ncol(datHMF)-2)*length(unique(datHMF[,2]))      # number of features times number of annotations
    Mmin <- ceiling(Np/2)     # minimum number of modulating functions
    Mmax <- ceiling(Msc*Np/2) # maximum number of modulating functions
    if(Mmax < Mmin) {Mmax = Mmin + 10; print("Mmax set to default of Mmin + 10") } # default
    Mpos <- seq(1,Mmax,1)
    Mneg <- -Mpos
    if(is.null(Mall) == F){Mall = Mall}
    if(is.null(Mall) == T){Mall <- c(0,Mpos,Mneg)}
  out = list(Mmin=Mmin, Mmax=Mmax, Mall=Mall)
  return(out)
}


#' m-value information
#'
#' Get indices in the vector of m-values and select which m-value ranges to use
#' @param N_mRanges The number of m-value ranges to use.
#' @param mvals A set of m-values, the object generated by get_M_vals(). Defaults to NULL.
#' @param Msc A scale factor that determines the M value (M = max(m)) corresponding to the largest range of m-values.
#' Mmax <- ceiling(Msc * Np / 2) where Np = number of parameters.
#' Defaults to 2.
#' @param mr Indices for a specific m-value range.
#' @return A List with m-value information.
#' @examples mRange_best = get_M_indices(N_mRanges,mvals,Msc,indM_best)$mRange
get_M_indices = function(N_mRanges,mvals=NULL,Msc=2,mr=NULL, datHMF = datHMF){

  if(is.null(mvals) == F){mvals = mvals}
  if(is.null(mvals) == T){
    mvals = get_M_vals(datHMF=datHMF,Mall=NULL,Msc=2)
  }
  Mall = mvals$Mall
  Mmin = mvals$Mmin
  Mmax = mvals$Mmax

  ind0 <- which(Mall==0)
  indp1 <- which(Mall==1)
  indn1 <- which(Mall==-1)
  indMin <- which(Mall == Mmin)
  indMax <- which(Mall == Mmax)
  Mspectrum <- c(1:(Mmax - Mmin+ 1))
  indM <- floor( seq(1, length(Mspectrum), length.out=N_mRanges) )

  if(is.null(mr)==T){
    mRange=NULL
    mRangeLen=NULL
  }
  if(is.null(mr)==F){
    mRange=c(ind0, c(indp1:(indMin+mr-1)), c(indn1:(indn1+indMin+mr-3)))
    mRangeLen <- c()
    for (mm in 1:length(indM)) {
      # m-indices for a given range of m-values
      mrang <- indM[mm]
      mRange2 <- c(ind0, c(indp1:(indMin+mrang-1)), c(indn1:(indn1+indMin+mrang-3)))
      mRangeLen = c(mRangeLen, length(mRange2))
    }
  }

  out = list(Mspectrum=Mspectrum, indM=indM, ind0=ind0, mRangeLen=mRangeLen, mRange=mRange)
  return(out)
} # end get_M_indices



#' Select M ranges
#'
#' Select a set of values for M = max(m)
#' @param datHMF The data matrix with scaled measurements measurements. Scaled center measurements are provided for
#' each measurement/annotation (e.g., gene/organ) combination. This format is produced by scale_zeroOne(). Defaults to datHMF.
#' @param Mall A set of m-values. Defaults to NULL.
#' @param Msc A scale factor that determines the M value (M = max(m)) corresponding to the largest range of m-values.
#' Mmax <- ceiling(Msc * Np / 2) where Np = number of parameters.
#' Defaults to 2.
#' @param N_mRanges The number of m-value ranges to use.
#' @param mvals A set of m-values, the object generated by get_M_vals(). Defaults to NULL.
#' @param mr Indices for a specific m-value range.
#' @return A List with m-value information.
#' @export
#' @examples mvals = get_M_vals(datHMF, Msc=2)
##### select M ranges
get_M_info <- function(datHMF=datHMF,Mall=NULL,Msc=2,N_mRanges,minfo=NULL,mr=NULL){
  
  if(is.null(minfo) == F){minfo = minfo}
  if(is.null(minfo) == T){
    if(Msc != 2){Msc = Msc}   # scale factor, user defined input
    Np <- (ncol(datHMF)-2)*length(unique(datHMF[,2]))      # number of features times number of annotations
    Mmin <- ceiling(Np/2)     # minimum number of modulating functions
    Mmax <- ceiling(Msc*Np/2) # maximum number of modulating functions
    if(Mmax < Mmin) {Mmax = Mmin + 10; print("Mmax set to default of Mmin + 10") } # default
    Mpos <- seq(1,Mmax,1)
    Mneg <- -Mpos
    # if(is.null(Mall) == F){Mall = Mall}
    # if(is.null(Mall) == T){Mall <- c(0,Mpos,Mneg)}
    Mall <- c(0,Mpos,Mneg)
    
    ind0 <- which(Mall==0)
    indp1 <- which(Mall==1)
    indn1 <- which(Mall==-1)
    indMin <- which(Mall == Mmin)
    indMax <- which(Mall == Mmax)
    Mspectrum <- c(1:(Mmax - Mmin+ 1))
    indM <- floor( seq(1, length(Mspectrum), length.out=N_mRanges) )  
  }else{
      Mmin = minfo$Mmin
      Mmax = minfo$Mmax
      Mall = minfo$Mall
      Mspectrum = minfo$Mspectrum
      indM = minfo$indM
      ind0 = minfo$ind0
      indp1 <- which(Mall==1)
      indn1 <- which(Mall==-1)
      indMin <- which(Mall == Mmin)
      indMax <- which(Mall == Mmax)
    }
  
  
  
  if(is.null(mr)==T){
    mRange=NULL
    mRangeLen=NULL
  }
  if(is.null(mr)==F){
    mRange=c(ind0, c(indp1:(indMin+mr-1)), c(indn1:(indn1+indMin+mr-3)))
    mRangeLen <- c()
    for (mm in 1:length(indM)) {
      # m-indices for a given range of m-values
      mrang <- indM[mm]
      mRange2 <- c(ind0, c(indp1:(indMin+mrang-1)), c(indn1:(indn1+indMin+mrang-3)))
      mRangeLen = c(mRangeLen, length(mRange2))
    }
  }
  
  
  out = list(Mmin=Mmin, Mmax=Mmax, Mall=Mall,Mspectrum=Mspectrum, indM=indM, ind0=ind0, mRangeLen=mRangeLen, mRange=mRange)
  return(out)
}
