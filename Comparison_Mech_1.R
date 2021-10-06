###############################################################
##############       Simulation Study for The     #############
##############         Comparison of ptFCE        #############
##############       with Existing Approaches     #############
##############          Using Mechanism 1         #############
###############################################################

library("mgcv")

###################################
## 1.1 Basic parameters for data
##     and simulations
###################################

# task_FC correlation
task_FC_seq=c(0.4, 0.6)

# Repeatition time
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
for (i in 1:length(t.obs)) {
  N.values[i]=N(t.obs[i])
  N.lf.values[i]=N.lf(t.obs[i])
  N.lh.values[i]=N.lh(t.obs[i])
  N.rf.values[i]=N.rf(t.obs[i])
  N.t.values[i]=N.t(t.obs[i])
  
  h.1.values[i]=h_1_function(t.obs[i])
  h.2.values[i]=h_2_function(t.obs[i])
  h.3.values[i]=h_3_function(t.obs[i])
}

h.1.conv.N=convolve(N.values, h.1.values)
h.2.conv.N=convolve(N.values, h.2.values)
h.3.conv.N=convolve(N.values, h.3.values)
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
## 2 Simulation iteraton 
##   (for n_simulation = 500 runs)
###################################

est_results_matrix=matrix(NA, ncol = length(task_FC_seq),
                          nrow = n_simulations)
naive_Pearson_mean=matrix(NA, ncol = length(task_FC_seq),
                          nrow = n_simulations)
naive_Pearson_median=matrix(NA, ncol = length(task_FC_seq),
                            nrow = n_simulations)
task_Pearson_mean=matrix(NA, ncol = length(task_FC_seq),
                         nrow = n_simulations)
task_Pearson_median=matrix(NA, ncol = length(task_FC_seq),
                           nrow = n_simulations)
beta_series_mean=matrix(NA, ncol = length(task_FC_seq),
                        nrow = n_simulations)
beta_series_median=matrix(NA, ncol = length(task_FC_seq),
                          nrow = n_simulations)
coherence_mean=matrix(NA, ncol = length(task_FC_seq),
                      nrow = n_simulations)
coherence_median=matrix(NA, ncol = length(task_FC_seq),
                        nrow = n_simulations)

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
  
  ## The ptFCE algorithm
  est_ptFCE_12=ptFCE(Y_k=Y.1.mat, Y_l = Y.2.mat, N=N.values, TR=0.72, freq_plot=FALSE)
  est_results_matrix[P, 1]=est_ptFCE_12$est
  est_ptFCE_23=ptFCE(Y_k=Y.2.mat, Y_l = Y.3.mat, N=N.values, TR=0.72, freq_plot=FALSE)
  est_results_matrix[P, 2]=est_ptFCE_23$est
  
  
  ## Naive Pearson correlation
  
  corr.vector_12=vector()
  corr.vector_23=vector()
  for (j in 1:N.participants) {
    corr.vector_12[j]=abs(cor(Y.1.mat[j,], Y.2.mat[j,]))
    corr.vector_23[j]=abs(cor(Y.2.mat[j,], Y.3.mat[j,]))
  }
  naive_Pearson_median[P, ]=c(median(corr.vector_12), median(corr.vector_23))
  naive_Pearson_mean[P, ]=c(mean(corr.vector_12), mean(corr.vector_23))
  
  
  ## Task Pearson correlation
  
  time.section.1=which(t.obs>=86.5 & t.obs<98.5)   
  time.section.2=which(t.obs>=162 & t.obs<174)
  time.section=c(time.section.1, time.section.2)
  Data.mat.1=Y.1.mat[, time.section]
  Data.mat.2=Y.2.mat[, time.section]
  Data.mat.3=Y.3.mat[, time.section]
  corr.vector_12=vector()
  corr.vector_23=vector()
  for (j in 1:N.participants) {
    corr.vector_12[j]=abs(cor(Data.mat.1[j,], Data.mat.2[j,]))
    corr.vector_23[j]=abs(cor(Data.mat.2[j,], Data.mat.3[j,]))
  }
  task_Pearson_median[P, ]=c(median(corr.vector_12), median(corr.vector_23))
  task_Pearson_mean[P, ]=c(mean(corr.vector_12), mean(corr.vector_23))
  
  
  ######################################################
  ## Beta-series regression
  ######################################################
  
  ## We apply the procedure described in Chapter 9.2 of Ashby (2019) 
  ## to estimate functional connectivity using beta-series regression, 
  ## except we ignore the “nuisance term” therein.
  
  N.lf.values=vector()
  N.lh.values=vector()
  N.rf.values=vector()
  N.values=vector()
  N.t.values=vector()
  HRF.values=vector()
  for (i in 1:length(t.obs)) {
    HRF.values=canonicalHRF(t.obs, verbose = FALSE)
    N.lf.values[i]=N.lf(t.obs[i])
    N.lh.values[i]=N.lh(t.obs[i])
    N.rf.values[i]=N.rf(t.obs[i])
    N.values[i]=N(t.obs[i])
    N.t.values[i]=N.t(t.obs[i])
  }
  
  indicator.N=as.numeric(N.values==1)
  indicator.N.lf=as.numeric(N.lf.values==1)
  indicator.N.lh=as.numeric(N.lh.values==1)
  indicator.N.rf=as.numeric(N.rf.values==1)
  indicator.N.t=as.numeric(N.t.values==1)
  
  zeros=rep(0, length(t.obs))
  
  design.matrix.N=matrix(NA, ncol = 1, nrow = length(t.obs))
  for (i in which(N.values==1)) {
    replacement=zeros
    replacement[i]=1
    replacement=convolve(replacement, HRF.values)
    design.matrix.N=cbind(design.matrix.N, replacement)
  }
  design.matrix.N=design.matrix.N[,-1]
  
  design.matrix.N.lf=matrix(NA, ncol = 1, nrow = length(t.obs))
  for (i in which(N.lf.values==1)) {
    replacement=zeros
    replacement[i]=1
    replacement=convolve(replacement, HRF.values)
    design.matrix.N.lf=cbind(design.matrix.N.lf, replacement)
  }
  design.matrix.N.lf=design.matrix.N.lf[,-1]
  
  design.matrix.N.lh=matrix(NA, ncol = 1, nrow = length(t.obs))
  for (i in which(N.lh.values==1)) {
    replacement=zeros
    replacement[i]=1
    replacement=convolve(replacement, HRF.values)
    design.matrix.N.lh=cbind(design.matrix.N.lh, replacement)
  }
  design.matrix.N.lh=design.matrix.N.lh[,-1]
  
  design.matrix.N.rf=matrix(NA, ncol = 1, nrow = length(t.obs))
  for (i in which(N.rf.values==1)) {
    replacement=zeros
    replacement[i]=1
    replacement=convolve(replacement, HRF.values)
    design.matrix.N.rf=cbind(design.matrix.N.rf, replacement)
  }
  design.matrix.N.rf=design.matrix.N.rf[,-1]
  
  design.matrix.N.t=matrix(NA, ncol = 1, nrow = length(t.obs))
  for (i in which(N.t.values==1)) {
    replacement=zeros
    replacement[i]=1
    replacement=convolve(replacement, HRF.values)
    design.matrix.N.t=cbind(design.matrix.N.t, replacement)
  }
  design.matrix.N.t=design.matrix.N.t[,-1]
  
  DESIGN.MATRIX=cbind(design.matrix.N, design.matrix.N.lf, design.matrix.N.lh, design.matrix.N.rf, design.matrix.N.t)
  
  
  Data.mat.1=Y.1.mat
  Data.mat.2=Y.2.mat
  Data.mat.3=Y.3.mat
  
  est.store_12=vector()
  est.store_23=vector()
  for (iter in 1:N.participants) {
    Y.1.regression=Data.mat.1[iter,]
    Y.2.regression=Data.mat.2[iter,]
    Y.3.regression=Data.mat.3[iter,]
    reg.1=lm(Y.1.regression~DESIGN.MATRIX)
    beta.est.1=as.vector(reg.1$coefficients[2:(dim(design.matrix.N)[2]+1)])
    reg.2=lm(Y.2.regression~DESIGN.MATRIX)
    beta.est.2=as.vector(reg.2$coefficients[2:(dim(design.matrix.N)[2]+1)])
    reg.3=lm(Y.3.regression~DESIGN.MATRIX)
    beta.est.3=as.vector(reg.3$coefficients[2:(dim(design.matrix.N)[2]+1)])
    est.store_12[iter]=abs(cor(beta.est.1, beta.est.2))
    est.store_23[iter]=abs(cor(beta.est.2, beta.est.3))
  }
  
  beta_series_median[P, ]=c(median(est.store_12), median(est.store_23))
  beta_series_mean[P, ]=c(mean(est.store_12), mean(est.store_23))
  
  
  ######################################################
  ## Coherence analysis
  ######################################################
  
  library("seewave")
  
  Data.mat.1=Y.1.mat
  Data.mat.2=Y.2.mat
  Data.mat.3=Y.3.mat
  est.store_12=vector()
  est.store_23=vector()
  for (iter in 1:N.participants) {
    coh.results=coh(Data.mat.1[iter, ], Data.mat.2[iter, ], f=1/TR, plot=FALSE)
    est.store_12[iter]=abs(median(coh.results[coh.results[,1]<0.15,2]))
    coh.results=coh(Data.mat.2[iter, ], Data.mat.3[iter, ], f=1/TR, plot=FALSE)
    est.store_23[iter]=abs(median(coh.results[coh.results[,1]<0.15,2]))
  }
  coherence_median[P, ]=c(median(est.store_12), median(est.store_23))
  coherence_mean[P, ]=c(mean(est.store_12), mean(est.store_23))
  
}



par(mfrow=c(3,3))
boxplot(est_results_matrix, main="ptFC")
abline(h=0.4, col="red", lty=2)
abline(h=0.6, col="red", lty=2)
boxplot(naive_Pearson_mean, main="naive_Pearson_mean")
boxplot(naive_Pearson_median, main="naive_Pearson_median")
boxplot(task_Pearson_mean, main="task_Pearson_mean")
boxplot(task_Pearson_median, main="task_Pearson_median")
boxplot(beta_series_mean, main="beta_series_mean")
boxplot(beta_series_median, main="beta_series_median")
boxplot(coherence_mean, main="coherence_mean")
boxplot(coherence_median, main="coherence_median")

write.csv(est_results_matrix, file = "Mech1_est_results_matrix.csv")
write.csv(naive_Pearson_mean, file = "Mech1_naive_Pearson_mean.csv")
write.csv(naive_Pearson_median, file = "Mech1_naive_Pearson_median.csv")
write.csv(task_Pearson_mean, file = "Mech1_task_Pearson_mean.csv")
write.csv(task_Pearson_median, file = "Mech1_task_Pearson_median.csv")
write.csv(beta_series_mean, file = "Mech1_beta_series_mean.csv")
write.csv(beta_series_median, file = "Mech1_beta_series_median.csv")
write.csv(coherence_mean, file = "Mech1_coherence_mean.csv")
write.csv(coherence_median, file = "Mech1_coherence_median.csv")


effectiveness=function(v){
  s=0
  if((v[1]<v[2])){ s=1 }
  return(s)
}

ptFC_rate=0
naive_Pearson_mean_rate=0
naive_Pearson_median_rate=0
task_Pearson_mean_rate=0
task_Pearson_median_rate=0
beta_series_mean_rate=0
beta_series_median_rate=0
coherence_mean_rate=0
coherence_median_rate=0
for (i in 1:n_simulations) {
  ptFC_rate=ptFC_rate+effectiveness(est_results_matrix[i,])
  naive_Pearson_mean_rate=naive_Pearson_mean_rate+effectiveness(naive_Pearson_mean[i,])
  naive_Pearson_median_rate=naive_Pearson_median_rate+effectiveness(naive_Pearson_median[i,])
  task_Pearson_mean_rate=task_Pearson_mean_rate+effectiveness(task_Pearson_mean[i,])
  task_Pearson_median_rate=task_Pearson_median_rate+effectiveness(task_Pearson_median[i,])
  beta_series_mean_rate=beta_series_mean_rate+effectiveness(beta_series_mean[i,])
  beta_series_median_rate=beta_series_median_rate+effectiveness(beta_series_median[i,])
  coherence_mean_rate=coherence_mean_rate+effectiveness(coherence_mean[i,])
  coherence_median_rate=coherence_median_rate+effectiveness(coherence_median[i,])
}
ptFC_rate=ptFC_rate/n_simulations
naive_Pearson_mean_rate=naive_Pearson_mean_rate/n_simulations
naive_Pearson_median_rate=naive_Pearson_median_rate/n_simulations
task_Pearson_mean_rate=task_Pearson_mean_rate/n_simulations
task_Pearson_median_rate=task_Pearson_median_rate/n_simulations
beta_series_mean_rate=beta_series_mean_rate/n_simulations
beta_series_median_rate=beta_series_median_rate/n_simulations
coherence_mean_rate=coherence_mean_rate/n_simulations
coherence_median_rate=coherence_median_rate/n_simulations

RATE_vector=c(ptFC_rate, naive_Pearson_mean_rate, naive_Pearson_median_rate, task_Pearson_mean_rate, task_Pearson_median_rate, beta_series_mean_rate, beta_series_median_rate, coherence_mean_rate, coherence_median_rate)
write.csv(RATE_vector, file = "Mech1_RATE_vector.csv")
