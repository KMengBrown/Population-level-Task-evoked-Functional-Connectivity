###############################################################
###############################################################
##############  ptFCE and AMUSE-ptFCE Algorithms  #############
###############################################################
###############################################################



###############################################################
##############            Section 1              ##############
##############         ptFCE Algorithm           ############## 
###############################################################



############################
## 1.1 Preparations
############################


## The Fourier transform returning 
## (1/n) * sum_{k=1}^n z[k]*exp(-2*pi*1i*(k-1)*(h-1)/n)
my.fft=function(x, TR.input){ return( fft(x)/length(x) ) }


## The following function provides the time-wise demean procedure.
center_scale <- function(x) { scale(x, scale = FALSE) }


## Project integers to a modulo group.
## This projection is prepared for the function periodic extension.
## e.g., 
## x=(-7):7
## for (i in 1:length(x)) { print(modulo.group.function(x[i], 5)) }
## Then [ -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7 ] becomes
## [3 4 5 1 2 3 4 5 1 2 3 4 5 1 2].
modulo.group.function=function(a, n){
  if(a%%n==0){ resp=n }else{ resp=a%%n }
  return(resp)
} 


## For a given sequence (x[1], x[2],...,x[n]), we first periodically extend it
## to be a periodic function and then shift its time-variate by a given "shift.size."
periodically.shift=function(x, shift.size){
  resp=vector()
  for (t in 1:length(x)) { resp[t]=x[modulo.group.function(t-shift.size, length(x))] }
  return(resp)
}


## Auto-Covariance Function: it returns the auto-covariance
## of a pair of stochastic processes.
Auto_covariance=function(Signal.1, Signal.2){
  
  # Both 'Signal.1' and 'Signal.2' are matrices.
  # Each row of these matrices is a realization of a stochastic process.
  # Each column corresponds to observations at a specific time point.
  # Each row dimension is the number of participants.
  # Each column dimension is the number fo time points.
  
  resp=vector() 
  for (iter in 0:(dim(Signal.1)[2]-1)) {
    vec.storage=vector()
    for (iter1 in 1:dim(Signal.1)[2]) { 
      vec.storage[iter1]=mean(Signal.1[, iter1]*Signal.2[, modulo.group.function(iter1+iter, dim(Signal.1)[2])])
    }
    resp[iter+1]=mean(vec.storage) 
  }
  return(resp)
}



############################
## 1.2 Semi-finished ptFCE Algorithm
##     (when underlying reference
##     signals are available)
############################


Semi_finished_ptFCE=function(Y_k, Y_l, 
               R_k, R_l,
               TR,
               cut_off=0.15,
               freq_plot=TRUE){
  
  # Y_k and Y_l: They present task-evoked BOLD signals.
  #            * They are matrices.
  #            * Each row is a time series, which is a realization of a stochastic process.
  #            * The row dimension is the number of participants, and the column dimension
  #              is the number of time points.
  # R_k and R_l: They present reference signals.
  #            * They are matrices.
  #            * Each row is a time series, which is a realization of a stochastic process.
  # TR:          Repeatition time.
  # cut_off:     We are interested in the frequencies between 0 and cut_off. Default is 0.15
  #              as HRFs act as band--pass filters (typically 0-0.15 Hz).
  # freq_plot:   If TRUE, the ptFCE() function plots estimates across all Fourier frequencies 
  #              in ( 0, 1/(2*TR) ). Default = TRUE.
  # Output:      An estimation of ptFC corr(beta_k, beta_l).
  
  # Centralization for zero-mean
  Y_k=apply(Y_k, 2, center_scale)
  Y_l=apply(Y_l, 2, center_scale)
  R_k=apply(R_k, 2, center_scale)
  R_l=apply(R_l, 2, center_scale)
  
  if( (dim(Y_k)[1]!=dim(Y_l)[1]) | (dim(Y_k)[2]!=dim(Y_l)[2])){
    print("The dimensionalities of input matrices are not compatible.")
  }else{
    
    ###### The key chunk starts #####
    
    A_k.shifted=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
    A_l.shifted=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
    R_k.shifted=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
    R_l.shifted=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
    for (iter in 1:dim(Y_k)[1]) {  
      # Introduce uniformly distributed time-shift to each row of A_k/A_l.
      # Seed is forbidden!!!
      shift=sample(1:dim(Y_k)[2], 1)        
      A_k.shifted[iter,]=periodically.shift(Y_k[iter,], shift.size = shift)
      A_l.shifted[iter,]=periodically.shift(Y_l[iter,], shift.size = shift)
      R_k.shifted[iter,]=periodically.shift(R_k[iter,], shift.size = shift)
      R_l.shifted[iter,]=periodically.shift(R_l[iter,], shift.size = shift)
    }
    
    Auto.cov.kl=Auto_covariance(A_k.shifted, A_l.shifted) 
    Auto.cov.kk=Auto_covariance(A_k.shifted, A_k.shifted) 
    Auto.cov.ll=Auto_covariance(A_l.shifted, A_l.shifted) 
    Auto.cov.kl.baseline=Auto_covariance(R_k.shifted, R_l.shifted)
    Auto.cov.kk.baseline=Auto_covariance(R_k.shifted, R_k.shifted)
    Auto.cov.ll.baseline=Auto_covariance(R_l.shifted, R_l.shifted)
    
  }
  ###### The key chunk ends #####
  
  auto_cov=list(A_kl=Auto.cov.kl, A_kk=Auto.cov.kk, A_ll=Auto.cov.ll)
  
  numerator=Mod(my.fft(Auto.cov.kl-Auto.cov.kl.baseline))
  denominator.kk=sqrt(Mod(my.fft(Auto.cov.kk-Auto.cov.kk.baseline)))
  denominator.ll=sqrt(Mod(my.fft(Auto.cov.ll-Auto.cov.ll.baseline)))
  
  vector.all.freq=numerator/(denominator.kk*denominator.ll)
  vector.all.freq=pmin(vector.all.freq, 1)
  vector.all.freq=pmax(vector.all.freq, 0)
  frequencies=(0:(dim(Y_k)[2]-1))/(TR*dim(Y_k)[2])
  ind.storage=which(frequencies<cut_off)
  # The following line restricts the estimation result in [0,1].
  opt.est=max(c(0, min(c(median(vector.all.freq[ind.storage]), 1))))
  
  if(freq_plot==TRUE){
    index.storage=which(frequencies<=1/(2*TR))
    plot(frequencies[index.storage], vector.all.freq[index.storage], 
         type = "l", col="blue",
         xlab = "Fourier Frequencies", ylab = "ptFC Estimates",
         main = " ")
    lines(frequencies[ind.storage], vector.all.freq[ind.storage],
          col="orange", type = "l", lty=1, lwd=3)
    abline(h=opt.est, lty=3, col="red", lwd=3)
    legend("topright", c("Est(Freq)", "Opt Est"), col=c("blue", "red"), lty=c(1,3))
    
  }
  
  resp=list(est=opt.est,
            est_all_freq=vector.all.freq, 
            auto_cov_factors=auto_cov)
  
  # output: "est" gives the desired output, i.e., an estimate of ptFC.
  
  return(resp)
  
}




###############################################################
##############            Section 2              ##############
##############      AMUSE-ptFCE Algorithm        ############## 
###############################################################


############################
## 2.1 Packages
############################

library("MASS")
library("JADE")
library("neuRosim")


############################
## 2.2 ptFCE Algorithm
############################

ptFCE=function(Y_k, Y_l, N, TR,
                     freq_plot=TRUE, 
                     cut_off=0.15){
  
  # Y_k and Y_l: They present task-evoked BOLD signals.
  #            * They are matrices.
  #            * Each row is a time series, which is a realization of a stochastic process.
  # N：          The stimulus signal. 
  #              Specifically, N herein is a vector, whose each entry is the value of the 
  #              stimulus signal at the corresponding time point.
  # TR:          Repeatition time.
  # freq_plot and cut_off are explained in ptFCE().
  
  # Observation time points
  t.variate=((1:(dim(Y_k)[2]))-1)*TR
  
  # Canonical HRF values at observation time points
  cano_HRF=canonicalHRF(t.variate, verbose = FALSE)
  
  # Convolution
  N_conv_h=convolve(N, cano_HRF)
  
  N_conv_h.shifted=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
  Y_k.shifted=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
  Y_l.shifted=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
  for (iter in 1:dim(Y_k)[1]) {
    shift=sample(1:dim(Y_k)[2], 1)
    N_conv_h.shifted[iter,]=periodically.shift(N_conv_h, shift.size = shift)
    Y_k.shifted[iter,]=periodically.shift(Y_k[iter,], shift.size = shift)
    Y_l.shifted[iter,]=periodically.shift(Y_l[iter,], shift.size = shift)
  }
  
  Y_k_store=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
  for(iter in 1:dim(Y_k)[1]){
    X_store=matrix(NA, nrow = dim(Y_k)[2], ncol = 2)
    X_store[,1]=Y_k.shifted[iter, ]
    X_store[,2]=N_conv_h.shifted[iter,]-mean(N_conv_h)
    # The "mean(N_conv_h)" herein is the constant "C_k" in our paper.
    res1<-AMUSE(X_store)
    r_correlation=vector()
    r_correlation[1]=cor(res1$S[,1], N_conv_h.shifted[iter,]-mean(N_conv_h))
    r_correlation[2]=cor(res1$S[,2], N_conv_h.shifted[iter,]-mean(N_conv_h))
    A_hat_k=ginv(res1$W)
    Y_k_store[iter,]=A_hat_k[1, max(which(max(r_correlation)==r_correlation))]*res1$S[, max(which(max(r_correlation)==r_correlation))]
  }
  
  Y_l_store=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
  for(iter in 1:dim(Y_k)[1]){
    X_store=matrix(NA, nrow = dim(Y_k)[2], ncol = 2)
    X_store[,1]=Y_l.shifted[iter, ]
    X_store[,2]=N_conv_h.shifted[iter,]-mean(N_conv_h)
    res1<-AMUSE(X_store)
    r_correlation=vector()
    r_correlation[1]=cor(res1$S[,1], N_conv_h.shifted[iter,]-mean(N_conv_h))
    r_correlation[2]=cor(res1$S[,2], N_conv_h.shifted[iter,]-mean(N_conv_h))
    A_hat_k=ginv(res1$W)
    Y_l_store[iter,]=A_hat_k[1, max(which(max(r_correlation)==r_correlation))]*res1$S[, max(which(max(r_correlation)==r_correlation))]
  }
  
  R_k=matrix(0, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
  R_l=R_k
  
  resp=Semi_finished_ptFCE(Y_k_store, Y_l_store, R_k, R_l, TR, freq_plot = freq_plot, cut_off = cut_off)
  return(resp)
  
}


############################
## End
############################

