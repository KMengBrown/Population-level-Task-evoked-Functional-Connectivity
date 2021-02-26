###############################################################
##############  Simulation Study for ptFCE and    #############
##############     AMUSE-ptFCE Algorithms in      #############
##############              Table 2               #############
###############################################################


###################################
## 1 Basic parameters for data
###################################

## The true underlying ptFC
corr.vector=c(0, 0.25, 0.5, 0.75, 1)

# Repeatition time
TR=0.72

#Observation times
t.obs=(0:283)*TR

## Number of participants
N.participants=308


###################################
## 2 Stimulus signals
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
## 3 HRF
###################################

# HRF at the k-th node is of the canonical form
h_k_function=function(t){ return(canonicalHRF(t, verbose = FALSE)) }

# HRF at the l-th node is noncanonical.
parameter.list=list(a1=10, a2=15, b1=0.9, b2=0.9, c=0.35)
h_l_function=function(t){ return(canonicalHRF(t, parameter.list, verbose = FALSE)) }


###################################
## 4 Reaction delay times
###################################

t.0.k=0  
t.0.l=0 


###################################
## 5 Convolutions N * HRF
###################################

N.values=vector()
N.lf.values=vector()
N.lh.values=vector()
N.rf.values=vector()
N.t.values=vector()

h.k.values=vector()
h.l.values=vector()
for (i in 1:length(t.obs)) {
  N.values[i]=N(t.obs[i])
  N.lf.values[i]=N.lf(t.obs[i])
  N.lh.values[i]=N.lh(t.obs[i])
  N.rf.values[i]=N.rf(t.obs[i])
  N.t.values[i]=N.t(t.obs[i])
  
  h.k.values[i]=h_k_function(t.obs[i])
  h.l.values[i]=h_l_function(t.obs[i])
}

h.k.conv.N=convolve(N.values, h.k.values)
h.l.conv.N=convolve(N.values, h.l.values)
h.k.conv.N.lf=convolve(N.lf.values, h.k.values)
h.l.conv.N.lf=convolve(N.lf.values, h.l.values)
h.k.conv.N.lh=convolve(N.lh.values, h.k.values)
h.l.conv.N.lh=convolve(N.lh.values, h.l.values)
h.k.conv.N.rf=convolve(N.rf.values, h.k.values)
h.l.conv.N.rf=convolve(N.rf.values, h.l.values)
h.k.conv.N.t=convolve(N.t.values, h.k.values)
h.l.conv.N.t=convolve(N.t.values, h.l.values)
h.k.conv.N.shifted=vector()
h.l.conv.N.shifted=vector()
for (i in 1:length(t.obs)) {
  h.k.conv.N.shifted[i]=h.k.conv.N[modulo.group.function(i-t.0.k, length(t.obs))]
  h.l.conv.N.shifted[i]=h.l.conv.N[modulo.group.function(i-t.0.l, length(t.obs))]
}



###################################
## 6 Simulate for I = 100 times
###################################

est.with.R.vector=vector()
est.without.R.vector=vector()
sample.cor.vector=vector()
order.est.with.R=matrix(NA, ncol = 5, nrow = 1)
order.est.without.R=matrix(NA, ncol = 5, nrow = 1)
order.sample.cor=matrix(NA, ncol = 5, nrow = 1)


for (I in 1:100) {
  
  print(I)
  
  for (ITER in 1:5) {
    
    ###################################
    ## Random effect 
    ## coefficients beta's
    ###################################
    
    corr.true=corr.vector[ITER]
    V1=2
    V2=3
    V12=sqrt(V1)*sqrt(V2)*corr.true
    V.beta=cbind(c(V1, V12), c(V12, V2))
    beta.interest=rmvn(N.participants, mu=c(0,0), V=V.beta)
    
    corr.other=0.3
    V1=2
    V2=3
    V12=sqrt(V1)*sqrt(V2)*corr.other
    V.beta.other=cbind(c(V1, V12), c(V12, V2))
    beta.lf=rmvn(N.participants, mu=c(0,0), V=V.beta.other)
    beta.lh=rmvn(N.participants, mu=c(0,0), V=V.beta.other)
    beta.rf=rmvn(N.participants, mu=c(0,0), V=V.beta.other)
    beta.t=rmvn(N.participants, mu=c(0,0), V=V.beta.other)
    
    
    ###################################
    ## Random noise
    ###################################
    
    v1.noise=1
    v2.noise=1
    cor.noise=0.2
    var.noise=cbind(c(v1.noise, sqrt(v1.noise)*sqrt(v2.noise)*cor.noise), c(sqrt(v1.noise)*sqrt(v2.noise)*cor.noise, v2.noise))

    
    ###################################
    ## Generating task-evoked and
    ## reference signals
    ###################################
    
    Y.k.mat=matrix(NA, nrow = N.participants, ncol = length(t.obs))
    Y.l.mat=matrix(NA, nrow = N.participants, ncol = length(t.obs))
    R.k.mat=matrix(NA, nrow = N.participants, ncol = length(t.obs))
    R.l.mat=matrix(NA, nrow = N.participants, ncol = length(t.obs))
    for (iter in 1:N.participants) {
      
      Noise=rmvn(length(t.obs), mu=rep(0, 2), V=var.noise)
      Y.k.mat[iter,] = 9000 + beta.interest[iter,1]*h.k.conv.N.shifted + beta.lf[iter,1]*h.k.conv.N.lf+beta.lh[iter,1]*h.k.conv.N.lh+beta.rf[iter,1]*h.k.conv.N.rf+beta.t[iter,1]*h.k.conv.N.t + Noise[,1]
      Y.l.mat[iter,] = 9000 + beta.interest[iter,2]*h.l.conv.N.shifted + beta.lf[iter,2]*h.l.conv.N.lf+beta.lh[iter,2]*h.l.conv.N.lh+beta.rf[iter,2]*h.l.conv.N.rf+beta.t[iter,2]*h.l.conv.N.t + Noise[,2]
      
      Noise1=rmvn(length(t.obs), mu=rep(0, 2), V=var.noise)
      R.k.mat[iter,] = 9000 + beta.lf[iter,1]*h.k.conv.N.lf+beta.lh[iter,1]*h.k.conv.N.lh+beta.rf[iter,1]*h.k.conv.N.rf+beta.t[iter,1]*h.k.conv.N.t + Noise1[,1]
      R.l.mat[iter,] = 9000 + beta.lf[iter,2]*h.l.conv.N.lf+beta.lh[iter,2]*h.l.conv.N.lh+beta.rf[iter,2]*h.l.conv.N.rf+beta.t[iter,2]*h.l.conv.N.t + Noise1[,2]
      
    }
    
    
    ###################################
    ## Estimation using the
    ## ptFCE algorithm
    ###################################
    
    est.result=ptFCE(Y_k = Y.k.mat, Y_l = Y.l.mat, R_k = R.k.mat, R_l = R.l.mat, TR=TR, freq_plot = FALSE)
    est.with.R.vector[ITER]=est.result$est
    
    
    ###################################
    ## Estimation using the
    ## AMUSE-ptFCE algorithm
    ###################################
    
    est.without.R=AMUSE_ptFCE(Y_k=Y.k.mat, Y_l = Y.l.mat, N=N.values, TR=0.72, freq_plot=FALSE)
    est.without.R.vector[ITER]=est.without.R$est
    
    
    sample.cor.vector[ITER]=abs(cor(beta.interest[,1], beta.interest[,2]))
    
    
  }
  
  order.est.with.R=rbind(order.est.with.R, order(est.with.R.vector))
  order.est.without.R=rbind(order.est.without.R, order(est.without.R.vector))
  order.sample.cor=rbind(order.sample.cor, order(sample.cor.vector))
  
}


order.est.with.R=order.est.with.R[-1,]
order.est.without.R=order.est.without.R[-1,]
order.sample.cor=order.sample.cor[-1,]


###################################
## Simulation results for the 
## ptFCE and AMUSE-ptFCE algorithms
## in Table 2
###################################

## Correct grading rate of the ptFCE algorithm
for (i in 1:5) {print(mean(order.est.with.R[,i]==i))}

## Correct grading rate of the AMUSE-ptFCE algorithm
for (i in 1:5) {print(mean(order.est.without.R[,i]==i))}

for (i in 1:5) {print(mean(order.sample.cor[,i]==i))}


############################
## End
############################