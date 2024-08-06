#' The continuous Hartley transform
#'
#' This function performs linear interpolation between data points and implements the continuous Hartley transform
#' @param w This is omega_0 (Equation 4, Anderson et al., 2017). Defaults to NULL.
#' @param m Vector of slope values. Defaults to NULL.
#' @param b Vector of intercept values. Defaults to NULL.
#' @param tf Vector of final time points. Defaults to NULL.
#' @param ti Vector of initial time points. Defaults to NULL.
#' @param Harg0 If omega_0 is equal to zero (x=0), set omega_0 to Harg0 (w=Harg0). Defaults to NULL.
#' @return A vector with the transformation of the measurement data into the frequency domain.
#' @export
#' @examples This function is called by HMF1compute() and HMFcompute(). This function is not recommended to be used in isolation.
#' Hcompute(omega, slope, intercept, final time, initial time, H0)
#' @references \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005627}
Hcompute <- function(w=NULL,m=NULL,b=NULL,tf = NULL, ti= NULL,Harg0=NULL) { # continuous Hartley transform approximation
  H <- 0 # initiallized CHT, to be updated in loop over all intervals


  w[w==0] = Harg0
  #if (w==0) {w = Harg0}
  wtf  = tf %o% w
  wti  = ti %o% w
  mw2 = outer(m,w^2, "/")
  bw   = outer(b,w, "/")
  # Imtcoswtdt <- (mm/(w^2)) * ( (cos(w*tf)+w*tf*sin(w*tf)) - (cos(w*ti)+w*ti*sin(w*ti)) )
  #   Ibcoswtdt  <- (b/w) * ( sin(w*tf) - sin(w*ti) )
  #   Imtsinwtdt <- (mm/(w^2)) * ( (sin(w*tf)-w*tf*cos(w*tf)) - (sin(w*ti)-w*ti*cos(w*ti)) )
  #   Ibsinwtdt  <- (-b/w) * ( cos(w*tf) - cos(w*ti) )
  # browser()
  Imtcoswtdt <- (mw2) * ( (cos(wtf)+wtf*sin(wtf)) - (cos(wti)+wti*sin(wti)) )
  Ibcoswtdt  <- (bw)  * ( sin(wtf) - sin(wti) )
  Imtsinwtdt <- (mw2) * ( (sin(wtf)-wtf*cos(wtf)) - (sin(wti)-wti*cos(wti)) )
  Ibsinwtdt  <- (-bw) * ( cos(wtf) - cos(wti) )

  H <- colSums( rbind(colSums(Imtcoswtdt),colSums(Ibcoswtdt),colSums(Imtsinwtdt),colSums(Ibsinwtdt)))
  #H <- H + Imtcoswtdt + Ibcoswtdt + Imtsinwtdt + Ibsinwtdt
  return(H)
}


#' HMF spectral component
#'
#' This function computes the mth HMF spectral component of the measurement profile (Equation 5, Anderson et al., 2017).
#' @param meds Matrix with the first column containing measurement times and the second column containing the corresponding measurements. Defaults to NULL.
#' @param m The m value.
#' @param W0 The omega_0 term.
#' @param n Highest order of the system (n=1 for the first degree system that this analysis was designed for).
#' @param Harg0 If omega_0 is equal to zero (x=0), set omega_0 to Harg0 (w=Harg0).
#' @return The spectral component data.
#' @export
#' @examples This function is not recommended to be used in isolation.
#' HMFcompute(data, m, omega, n, H0)
#' @references \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005627}
HMFcompute <- function(meds,m,W0,n,Harg0) {  # spectral component, i=1
  #browser()
  # t <- meds[,1] # time
  # x <- meds[,2] # expression data, -ddCt (imputed)
  i <- 1        # order of the derivative, applied for dE/dt (i=1)
  HhatX = HhatY = vector(length = length(m))     # initialize Hhat, to be updated in loop

  xff = meds[-1         ,2]
  xii = meds[-nrow(meds),2]
  tff = meds[-1         ,1]
  tii = meds[-nrow(meds),1]
  slope = (xff-xii)/(tff-tii)
  intercept = xii - slope*tii

  for (j in c(0:n)) {
    # casargX <- (i*pi/2) * ((n+m-j)^i)
    # casX <- cos(casargX) - sin(casargX)
    HargX <- (n+m-j) * W0
    HX <- Hcompute(HargX,slope,intercept,tff,tii,Harg0)
    HhatX <- HhatX + ((-1)^j) * choose(n,j) * HX

    casargY <- (i*pi/2)
    casY <- cos(casargY) - sin(casargY)
    HargY <- ((-1)^i) * (n+m-j) * W0
    HY <- Hcompute(HargY,slope,intercept,tff,tii, Harg0)

    HhatY = HhatY + ((-1)^j) * choose(n,j) * casY * ((n+m-j)^i) * (W0^i) * HY

  }
  Hhat = data.frame(X = HhatX, Y= HhatY)
  return(Hhat)
}

#' X and Y data
#'
#' This function for getting the X and Y data (Equation 10, Anderson et al., 2017).
#' @param datHMF The data matrix with scaled measurements measurements. Scaled center measurements are provided for
#' each measurement/annotation (e.g., gene/organ) combination. This format is produced by scale_zeroOne(). Defaults to datHMF.
#' @param Mall A set of m-values. Defaults to NULL.
#' @param Msc A scale factor that determines the M value (M = max(m)) corresponding to the largest range of m-values.
#' Mmax <- ceiling(Msc * Np / 2) where Np = number of parameters.
#' @param n Highest order of the system (n=1 for the first degree system that this analysis was designed for).
#' @param Harg0 If omega_0 is equal to zero (x=0), set omega_0 to Harg0 (w=Harg0).
#' @return The X and Y transformations obtained using Equations 4-6 for all m-value ranges
#' @export
#' @examples XY = getXY(datHMF,Mall=Mall,Msc=Msc,n=n,Harg0=Harg0)
#' @references \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005627}
getXY = function(datHMF=datHMF,Mall=NULL,Msc=2,n=1,Harg0=NULL){

  W0 = 2*pi / (max(unique(datHMF[,1]),na.rm = T)-min(unique(datHMF[,1]),na.rm = T))

  if(is.null(Mall) == F){Mall = Mall}
  if(is.null(Mall) == T){
    mvals = get_M_vals(datHMF=datHMF,Mall=NULL,Msc=2)
    Mall = mvals$Mall
  }

  #X = Y = array( c(0), c(ncol(datHMF)-2, length(unique(datHMF[,2])), length(Mall) ) )

  dat.long = melt(datHMF ,id=c(2,1),variable.name = "Gene",value.name = "Exp",factorsAsStrings = TRUE) %>% dplyr::select(2,1,3,4) %>% arrange_at(.vars=2)

  grouped_XY = dat.long %>% group_by_at(.vars = c(2,3)) %>% do({ #Use multidplyr parallel do #group by Organ+Gene

    # if (length(.$Age)<4){paste("not enough ages","_",.$Strain[1],"_",.$Organ[1],"_",.$Gene[1])} #hardcoded <4
    meds = data.frame(.[,1],.[,4])
    meds_XY = HMFcompute(meds, Mall, W0, n, Harg0)
  }) %>% ungroup() %>% as.data.frame()   #Order: Mall, Gene, Organ
  #browser()
  dims.out = c(length(Mall),length(unique(grouped_XY[,2])),length(unique(grouped_XY[,1]))) 
  X = array(grouped_XY$X, dim = dims.out) #Mall,Gene,Organ
  Y = array(grouped_XY$Y, dim = dims.out) #Mall,Gene,Organ
  return(list(X=X,Y=Y))
} # end getXY()
