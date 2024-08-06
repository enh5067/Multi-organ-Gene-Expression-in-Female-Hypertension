#' System identification
#'
#' This function performs system identification to infer the connectivity and dynamics of a
#' linear network model (Anderson et al., 2017).
#' @param datHMF The data matrix with scaled measurements measurements. Scaled center measurements are provided for
#' each measurement/annotation (e.g., gene/organ) combination. This format is produced by scale_zeroOne(). Defaults to datHMF.
#' @param N_mRanges The number of m-value sets to evaluate.
#' @param alpha The L1 regularization terms for constraining the regression. Defaults to seq(0,1,0.2).
#' @param nlambda The number of L2 regularization terms for consideration. Defaults to 10.
#' @param Mall A set of m-values. Defaults to NULL.
#' @param Msc A scale factor that determines the M value (M = max(m)) corresponding to the largest range of m-values.
#' Mmax <- ceiling(Msc * Np / 2) where Np = number of parameters.
#' @param n Highest order of the system (n=1 for the first degree system that this analysis was designed for).
#' @param dt Simulation time step
#' @param epsilon Term for computing simulation error. Defaults to 10^(-10).
#' @param Harg0 If omega_0 is equal to zero (x=0), set omega_0 to Harg0 (w=Harg0).
#' @return Error data for all simulations, simulation data for the best simulation, a parameter matrix for the best simulation.
#' @export
#' @examples err0 = HMF_fit(datHMF,N_mRanges=2,dt=0.1,alpha=c(0.2,0.6))
#' @references \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005627}
HMF_fit = function(datHMF=NULL,N_mRanges=10,alpha=NULL,nlambda=NULL,Mall=NULL,Msc=2,n=1,dt=0.1,epsilon=10^(-10),Harg0=0.000001, path = NULL){
  #setwd(path)
  library(magrittr)
  library(dplyr)
  library(reshape2)
  library(deSolve)
  library(Matrix)
  
  # regularization terms
  if(is.null(alpha)==F){alph = alpha}
  if(is.null(alpha)==T){alph  = seq(0,1,0.2)}
  if(is.null(nlambda)==F){nlambda = nlambda}
  if(is.null(nlambda)==T){nlambda  = 10}
  
  
  #datHMF = datHMF %>% group_by_at(.vars = c(2,1))  %>% summarise_at(.vars = -c(1,2), funs(mean)) %>% arrange_at(.vars = c(1,2)) %>% ungroup() %>% select_at(c(2,1,3:ncol(datHMF))) %>% as.data.frame()
  
  #datHMF = datHMF %>% group_by_at(.vars = c(2,1))  %>% summarise_at(.vars = -c(1,2), funs(mean)) %>% ungroup() %>% select_at(c(2,1,3:ncol(datHMF))) %>% as.data.frame()
  
  # X and Y data, dims: gene, organ, m value (in Mall)
  # mvals = get_M_vals(datHMF=datHMF,Msc=Msc, Mall = Mall)
  # Mall = mvals$Mall
  minfo = get_M_info(datHMF=datHMF,Msc=Msc, Mall = Mall, N_mRanges = N_mRanges)
  Mall = minfo$Mall
  Mspectrum = minfo$Mspectrum # indices in the vector of m-values
  indM = minfo$indM # select which m-value ranges to plot
  
  XY = getXY(datHMF,Mall=Mall,Msc=Msc,n=n,Harg0=Harg0)
  
  # mspecs = get_M_indices(N_mRanges, mvals, datHMF) # m-value indices
  # Mspectrum = mspecs$Mspectrum # indices in the vector of m-values
  # indM = mspecs$indM # select which m-value ranges to plot
  

  # genes  = c(1,2,4,7,10,13,16,19,22,3,6,9,12,15,18,21,24,11,14,17,20,23)
  # datHMF.old = datHMF
  # datHMF = datHMF[,c(1,2,2+genes)]
  
  # model initial conditions for simulation
  ic = get_init_condits(datHMF)   #Move later with getXY
  
  # total number of model parameters
  numPar <- ((ncol(datHMF)-2)^2)*(length(unique(datHMF[,2]))^2)
  
  
    # loop through alpha
  simulation_error = list()
  parameters = list()
  for (alp in 1:length(alph)) {
    
    message(paste("alpha = ", alph[alp]))


    # space for saving rate constant data
    #fitCoefs <- array(c(0),c(numPar,nlambda,length(indM))) # fit coefficients
    MLAM <- matrix(c(0),nlambda*length(indM),5)
    
    # loop through every set of m-value ranges
    # implement parameter estimation--regression for each set of m-values
    for (mm in 1:length(indM)) {
      
      message(paste("M-ind = ", mm ))
      
      
      # m-indices for a given range of m-values
      #mRange = get_M_indices(N_mRanges = N_mRanges,mvals = mvals,Msc = Msc,mr = indM[mm],datHMF = datHMF)$mRange
      mRange = get_M_info(N_mRanges = N_mRanges,Mall = Mall,Msc = Msc,mr = indM[mm],datHMF = datHMF)$mRange
      Yreg = get_Yreg(datHMF,XY$Y,numPar,mRange)
      Xreg = get_Xreg(datHMF,XY$X,numPar,mRange)
      
      # fit the regression model using regularization
      # alpha <- alph[alp]
      fit <- glmnet::glmnet( Xreg,Yreg, family="gaussian",
                             alpha=alph[alp], intercept=F, standardize=F, nlambda=nlambda)
      
      fit = get_rateConst(fit,nlambda)
      rateConst = fit$beta
     
      indLam <- ((mm-1)*nlambda+1) : (mm*nlambda)
      MLAM[indLam,1] <- mm # m-value index (1,2,3,...)
      MLAM[indLam,2] <- c(1:nlambda) # lambda index (1,2,3,...)
      MLAM[indLam,3] <- indM[mm] # which value from indM
      MLAM[indLam,4] <- fit$lambda # numeric lambda value
      colnames(MLAM) = c("m-value_index","lambda_index","indM_value","lambda_value","error")
    
      # model simulation - get errors
      tlen = max(unique(datHMF[,1])) - min(unique(datHMF[,1]))
      Nr = length(unique(datHMF[,2])) # number of organs (annotations)
      
      
      MLAM[indLam,5] <- sim_error(datHMF,Nr,tlen,rateConst,ic,dt,epsilon,numPar)
      
      rm(mRange, Xreg, Yreg, fit, rateConst)
   
    } # end mm loop (m-value ranges)
    
    # paraOut = list()
    #fitCoefs[,1:ncol(rateConst),mm] <- rateConst
    # for (mv in 1:length(indM)) {
    #   paraOut[[mv]] <- fitCoefs[,,mv]
    # }
    # parameters[[alp]] = paraOut
    simulation_error[[alp]] = MLAM
    
  } # end alp loop
  
  # select best simulation
  
  simulate_error = format_error(simulation_error,alph)
  best_params = select_best(simulate_error)
  
  # simulate the best model
  # best_params: col1=alpha, col2=mm, col4=indM[mm], col5=lambda
  indM_best = best_params$indM_value
  # mRange_best = get_M_indices(N_mRanges,mvals,Msc,indM_best, datHMF)$mRange
  mRange_best = get_M_info(datHMF = datHMF, minfo = minfo, mr = indM_best)$mRange
  tsim <- seq(0,tlen,by=dt)
  #tsim <- seq(0,tlen,length.out = 100)
  #dt = tsim[2]-tsim[1]
  
  best_sims = sim_best(best_params,mRange_best,tsim,numPar,XY,datHMF,ic, dt, nlambda = 10)
  best_sims$simulations[,1] = best_sims$simulations[,1] + min(datHMF[,1])
  
  return(list(errors = simulate_error, best_sim = best_sims$simulations, best_params = best_sims$param_matrix, 
              best_settings = best_params))
  
} # end HMF_fit()
