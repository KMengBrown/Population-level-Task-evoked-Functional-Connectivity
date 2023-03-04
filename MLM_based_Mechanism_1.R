###############################################################
##############        Simulation Study for        #############
##############     the MLM-based Approaches       #############
##############        Using Mechanism 1           #############
###############################################################

library("mgcv")

###################################
## 1.1 Basic parameters for data
##     and simulations
###################################

# task_FC correlation
task_FC_seq=c(0.4, 0.6)

# Repetition time
TR=0.72

#Observation times
t.obs=(0:283)*TR

## Number of participants
N.participants=308

## Number of simulations
n_simulations=500



###################################
## 1.2 Stimulus signals
###################################

N.lf=function(t){
  resp=0
  if(t>=71.35 & t<83.35){ resp=1 }
  if(t>=177.125 & t<189.125){ resp=1 }
  return(resp)
}

N.lh=function(t){
  resp=0
  if(t>=11 & t<23){ resp=1 }
  if(t>=116.63 & t<128.63){ resp=1 }
  return(resp)
}

N.rf=function(t){
  resp=0
  if(t>=26.13 & t<38.13){ resp=1 }
  if(t>=146.88 & t<158.88){ resp=1 }
  return(resp)
}

N.rh=function(t){
  resp=0
  if(t>=86.5 & t<98.5){ resp=1 }
  if(t>=162 & t<174){ resp=1 }
  return(resp)
}

N.t=function(t){
  resp=0
  if(t>=56.26 & t<68.26){ resp=1 }
  if(t>=131.75 & t<143.75){ resp=1 }
  return(resp)
}

## Stimulus signal of interest
N=function(t){ return(N.rh(t)) }

###################################
## 1.3 HRF
###################################

# (Biased) HRFs
parameter.list_1=list(a1=4, a2=10, b1=0.8, b2=0.8, c=0.4)
h_1_function=function(t){ return(canonicalHRF(t, parameter.list_1, verbose = FALSE)) }
parameter.list_2=list(a1=8, a2=14, b1=1, b2=1, c=0.3)
h_2_function=function(t){ return(canonicalHRF(t, parameter.list_2, verbose = FALSE)) }
parameter.list_1=list(a1=4, a2=10, b1=0.8, b2=0.8, c=0.4)
h_3_function=function(t){ return(canonicalHRF(t, parameter.list_1, verbose = FALSE)) }
h_canonical_function=function(t){ return(canonicalHRF(t, verbose = FALSE)) }

###################################
## 1.4 Reaction delay times
###################################

t.0.k=0  
t.0.l=0
# Latency times are usually much smaller than the TR in experiments of interest.


###################################
## 1.5 Convolutions N * HRF
###################################

N.values=vector()
N.lf.values=vector()
N.lh.values=vector()
N.rf.values=vector()
N.t.values=vector()

h.1.values=vector()
h.2.values=vector()
h.3.values=vector()
h.canonical.values=vector()
for (i in 1:length(t.obs)) {
  N.values[i]=N(t.obs[i])
  N.lf.values[i]=N.lf(t.obs[i])
  N.lh.values[i]=N.lh(t.obs[i])
  N.rf.values[i]=N.rf(t.obs[i])
  N.t.values[i]=N.t(t.obs[i])
  
  h.1.values[i]=h_1_function(t.obs[i])
  h.2.values[i]=h_2_function(t.obs[i])
  h.3.values[i]=h_3_function(t.obs[i])
  h.canonical.values[i]=h_canonical_function(t.obs[i])
}

h.1.conv.N=convolve(N.values, h.1.values)
h.2.conv.N=convolve(N.values, h.2.values)
h.3.conv.N=convolve(N.values, h.3.values)
h.canonical.conv.N.t=convolve(N.values, h.canonical.values)
h.1.conv.N.lf=convolve(N.lf.values, h.1.values)
h.2.conv.N.lf=convolve(N.lf.values, h.2.values)
h.3.conv.N.lf=convolve(N.lf.values, h.3.values)
h.1.conv.N.lh=convolve(N.lh.values, h.1.values)
h.2.conv.N.lh=convolve(N.lh.values, h.2.values)
h.3.conv.N.lh=convolve(N.lh.values, h.3.values)
h.1.conv.N.rf=convolve(N.rf.values, h.1.values)
h.2.conv.N.rf=convolve(N.rf.values, h.2.values)
h.3.conv.N.rf=convolve(N.rf.values, h.3.values)
h.1.conv.N.t=convolve(N.t.values, h.1.values)
h.2.conv.N.t=convolve(N.t.values, h.2.values)
h.3.conv.N.t=convolve(N.t.values, h.3.values)


###################################
## 2 Simulation iteration 
##   (for n_simulation = 500 runs)
###################################

V.beta=diag(c(2,3,2))
V.beta[1,2]=task_FC_seq[1]*sqrt(2*3)
V.beta[2,1]=V.beta[1,2]
V.beta[2,3]=task_FC_seq[2]*sqrt(2*3)
V.beta[3,2]=V.beta[2,3]

corr.other=0.3
V.beta.other=diag(c(2,3,2))
V.beta.other[1,2]=corr.other*sqrt(2*3)
V.beta.other[2,1]=corr.other*sqrt(2*3)
V.beta.other[2,3]=corr.other*sqrt(2*3)
V.beta.other[3,2]=corr.other*sqrt(2*3)

var.noise=diag(rep(30, 3))

correctness=vector()

for (P in 1:n_simulations) {
  
  print(P)
  
  
  Y.1.mat=matrix(NA, ncol = 284, nrow = N.participants)
  Y.2.mat=matrix(NA, ncol = 284, nrow = N.participants)
  Y.3.mat=matrix(NA, ncol = 284, nrow = N.participants)
  
  # The following for-chunk generates data.
  for (iter in 1:N.participants) {
    
    beta.interest=rmvn(N.participants, mu=rep(0, 3), V=V.beta)
    
    beta.lf=rmvn(N.participants, mu=rep(0, 3), V=V.beta.other)
    beta.lh=rmvn(N.participants, mu=rep(0, 3), V=V.beta.other)
    beta.rf=rmvn(N.participants, mu=rep(0, 3), V=V.beta.other)
    beta.t=rmvn(N.participants, mu=rep(0, 3), V=V.beta.other)
    
    Noise=rmvn(length(t.obs), mu=rep(0, 3), V=var.noise)
    
    Y.1.mat[iter,] = 9000 + beta.interest[iter,1]*h.1.conv.N + beta.lf[iter,1]*h.1.conv.N.lf+beta.lh[iter,1]*h.1.conv.N.lh+beta.rf[iter,1]*h.1.conv.N.rf+beta.t[iter,1]*h.1.conv.N.t + Noise[,1]
    Y.2.mat[iter,] = 9000 + beta.interest[iter,2]*h.2.conv.N + beta.lf[iter,2]*h.2.conv.N.lf+beta.lh[iter,2]*h.2.conv.N.lh+beta.rf[iter,2]*h.2.conv.N.rf+beta.t[iter,2]*h.2.conv.N.t + Noise[,2]
    Y.3.mat[iter,] = 9000 + beta.interest[iter,3]*h.3.conv.N + beta.lf[iter,3]*h.3.conv.N.lf+beta.lh[iter,3]*h.3.conv.N.lh+beta.rf[iter,3]*h.3.conv.N.rf+beta.t[iter,3]*h.3.conv.N.t + Noise[,3]
    
  }
  
  N=N.values
  
  Y_k=Y.1.mat
  Y_l=Y.2.mat
  Y_m=Y.3.mat
  
  t.variate=((1:(dim(Y_k)[2]))-1)*TR
  
  # Canonical HRF values at observation time points
  cano_HRF=canonicalHRF(t.variate, verbose = FALSE)
  
  # Convolution
  N_conv_h=convolve(N, cano_HRF)
  
  N_conv_h.shifted=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
  Y_k.shifted=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
  Y_l.shifted=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
  Y_m.shifted=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
  for (iter in 1:dim(Y_k)[1]) {
    shift=sample(1:dim(Y_k)[2], 1)
    N_conv_h.shifted[iter,]=periodically.shift(N_conv_h, shift.size = shift)
    Y_k.shifted[iter,]=periodically.shift(Y_k[iter,], shift.size = shift)
    Y_l.shifted[iter,]=periodically.shift(Y_l[iter,], shift.size = shift)
    Y_m.shifted[iter,]=periodically.shift(Y_m[iter,], shift.size = shift)
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
    A_hat_l=ginv(res1$W)
    Y_l_store[iter,]=A_hat_l[1, max(which(max(r_correlation)==r_correlation))]*res1$S[, max(which(max(r_correlation)==r_correlation))]
  }
  
  Y_m_store=matrix(NA, nrow = dim(Y_k)[1], ncol = dim(Y_k)[2])
  for(iter in 1:dim(Y_k)[1]){
    X_store=matrix(NA, nrow = dim(Y_k)[2], ncol = 2)
    X_store[,1]=Y_m.shifted[iter, ]
    X_store[,2]=N_conv_h.shifted[iter,]-mean(N_conv_h)
    res1<-AMUSE(X_store)
    r_correlation=vector()
    r_correlation[1]=cor(res1$S[,1], N_conv_h.shifted[iter,]-mean(N_conv_h))
    r_correlation[2]=cor(res1$S[,2], N_conv_h.shifted[iter,]-mean(N_conv_h))
    A_hat_m=ginv(res1$W)
    Y_m_store[iter,]=A_hat_m[1, max(which(max(r_correlation)==r_correlation))]*res1$S[, max(which(max(r_correlation)==r_correlation))]
  }

  Y_k_vector=as.vector(t(Y_k_store))
  Y_l_vector=as.vector(t(Y_l_store))
  Y_m_vector=as.vector(t(Y_m_store))
  
  independent_variables=c(NA)
  subject_id=c(NA)
  for (subjects_index in 1:N.participants) {
    independent_variables=c(independent_variables, h.canonical.conv.N.t)
    subject_id=c(subject_id, rep(subjects_index, 284))
  }
  subject_id=subject_id[-1]
  independent_variables=independent_variables[-1]
  
  MLM_based=cbind(Y_k_vector, Y_l_vector, Y_m_vector, subject_id, independent_variables)
  colnames(MLM_based)=c("Y_k_lmer", "Y_l_lmer", "Y_m_lmer", "IDs", "N_conv_HRF")
  MLM_based=as.data.frame(MLM_based)
  
  results_k=lmer(Y_k_lmer~N_conv_HRF+(N_conv_HRF|IDs)-1, data = MLM_based)
  beta_k=coef(results_k)
  beta_k=beta_k$IDs[,2]
  
  results_l=lmer(Y_l_lmer~N_conv_HRF+(N_conv_HRF|IDs)-1, data = MLM_based)
  beta_l=coef(results_l)
  beta_l=beta_l$IDs[,2]
  
  results_m=lmer(Y_m_lmer~N_conv_HRF+(N_conv_HRF|IDs)-1, data = MLM_based)
  beta_m=coef(results_m)
  beta_m=beta_m$IDs[,2]
  
  correctness[P]=cor(beta_k, beta_l)<cor(beta_l, beta_m)
  
}


mean(correctness)
1.96*sqrt( (mean(correctness)*(1-mean(correctness)))/n_simulations )
