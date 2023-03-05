
library("mgcv")
library("plot.matrix")

N_nodes=50
Connectivity_true=matrix(NA, ncol = N_nodes, nrow = N_nodes)
for (i in 1:N_nodes) {
  for (j in 1:i) {
    if(i==j){
      Connectivity_true[i,j]=1
    }else{
      if(j<=floor(N_nodes/2)){
        Connectivity_true[i,j]=runif(1, min = 0, max = 0.4)
        Connectivity_true[j,i]=Connectivity_true[i,j]
      }else{
        Connectivity_true[i,j]=runif(1, min = 0.6, max = 1)
        Connectivity_true[j,i]=Connectivity_true[i,j]
      }
    }
  }
}


###############################################################
#### Mechanism 1
###############################################################

Connectivity_estimated_1=matrix(NA, ncol = N_nodes, nrow = N_nodes)
for (node_k in 1:N_nodes) {
  for (node_l in 1:node_k) {
    print(paste(as.character(node_k), ", ", as.character(node_l)))
    if(node_k==node_l){
      Connectivity_estimated_1[node_k, node_l]=1
    }else{
      
      ###################################
      ## 1.1 Basic parameters for data
      ##     and simulations
      ###################################
      
      # task_FC correlation
      task_FC_seq=Connectivity_true[node_k, node_l]
      
      # Repeatition time
      TR=0.72
      
      #Observation times
      t.obs=(0:283)*TR
      
      ## Number of participants
      N.participants=308
      
      
      
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
      h.canonical.values=vector()
      for (i in 1:length(t.obs)) {
        N.values[i]=N(t.obs[i])
        N.lf.values[i]=N.lf(t.obs[i])
        N.lh.values[i]=N.lh(t.obs[i])
        N.rf.values[i]=N.rf(t.obs[i])
        N.t.values[i]=N.t(t.obs[i])
        
        h.1.values[i]=h_1_function(t.obs[i])
        h.2.values[i]=h_2_function(t.obs[i])
        h.canonical.values[i]=h_canonical_function(t.obs[i])
      }
      
      h.1.conv.N=convolve(N.values, h.1.values)
      h.2.conv.N=convolve(N.values, h.2.values)
      h.canonical.conv.N.t=convolve(N.values, h.canonical.values)
      h.1.conv.N.lf=convolve(N.lf.values, h.1.values)
      h.2.conv.N.lf=convolve(N.lf.values, h.2.values)
      h.1.conv.N.lh=convolve(N.lh.values, h.1.values)
      h.2.conv.N.lh=convolve(N.lh.values, h.2.values)
      h.1.conv.N.rf=convolve(N.rf.values, h.1.values)
      h.2.conv.N.rf=convolve(N.rf.values, h.2.values)
      h.1.conv.N.t=convolve(N.t.values, h.1.values)
      h.2.conv.N.t=convolve(N.t.values, h.2.values)
      
      
      ###################################
      ## 2 Simulation iteraton 
      ##   (for n_simulation = 500 runs)
      ###################################
      
      V.beta=diag(c(2,3))
      V.beta[1,2]=task_FC_seq[1]*sqrt(2*3)
      V.beta[2,1]=V.beta[1,2]
      
      corr.other=0.3
      V.beta.other=diag(c(2,3))
      V.beta.other[1,2]=corr.other*sqrt(2*3)
      V.beta.other[2,1]=corr.other*sqrt(2*3)
      
      var.noise=diag(rep(30, 2))
      
      
      Y.1.mat=matrix(NA, ncol = 284, nrow = N.participants)
      Y.2.mat=matrix(NA, ncol = 284, nrow = N.participants)
      
      # The following for-chunk generates data.
      for (iter in 1:N.participants) {
        
        beta.interest=rmvn(N.participants, mu=rep(0, 2), V=V.beta)
        
        beta.lf=rmvn(N.participants, mu=rep(0, 2), V=V.beta.other)
        beta.lh=rmvn(N.participants, mu=rep(0, 2), V=V.beta.other)
        beta.rf=rmvn(N.participants, mu=rep(0, 2), V=V.beta.other)
        beta.t=rmvn(N.participants, mu=rep(0, 2), V=V.beta.other)
        
        Noise=rmvn(length(t.obs), mu=rep(0, 2), V=var.noise)
        
        Y.1.mat[iter,] = 9000 + beta.interest[iter,1]*h.1.conv.N + beta.lf[iter,1]*h.1.conv.N.lf+beta.lh[iter,1]*h.1.conv.N.lh+beta.rf[iter,1]*h.1.conv.N.rf+beta.t[iter,1]*h.1.conv.N.t + Noise[,1]
        Y.2.mat[iter,] = 9000 + beta.interest[iter,2]*h.2.conv.N + beta.lf[iter,2]*h.2.conv.N.lf+beta.lh[iter,2]*h.2.conv.N.lh+beta.rf[iter,2]*h.2.conv.N.rf+beta.t[iter,2]*h.2.conv.N.t + Noise[,2]

      }
      
      est_ptFCE=ptFCE(Y_k=Y.1.mat, Y_l = Y.2.mat, N=N.values, TR=0.72, freq_plot=FALSE)
      Connectivity_estimated_1[node_k, node_l]=est_ptFCE$est
      Connectivity_estimated_1[node_l, node_k]=Connectivity_estimated_1[node_k, node_l]
      
    }
    
  }
}



###############################################################
#### Mechanism 2
###############################################################


Connectivity_estimated_2=matrix(NA, ncol = N_nodes, nrow = N_nodes)
for (node_k in 1:N_nodes) {
  for (node_l in 1:node_k) {
    print(paste(as.character(node_k), ", ", as.character(node_l)))
    if(node_k==node_l){
      Connectivity_estimated_2[node_k, node_l]=1
    }else{
      
      ###################################
      ## 1.1 Basic parameters for data
      ##     and simulations
      ###################################
      
      # Resting-state FC correlation
      resting_state_FC=0
      
      # task_FC correlation
      task_FC_seq=Connectivity_true[node_k, node_l]
      
      # Resting random noise
      V_resting=diag(rep(30, 2))
      
      # Task random noise
      V_task=diag(rep(30, 2))
      V_task[1,2]=task_FC_seq*30
      V_task[2,1]=V_task[1,2]
      
      # Repeatition time
      TR=0.72
      
      #Observation times
      t.obs=(0:283)*TR
      
      ## Number of participants
      N.participants=308
      
      
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
      h.canonical.values=vector()
      for (i in 1:length(t.obs)) {
        N.values[i]=N(t.obs[i])
        N.lf.values[i]=N.lf(t.obs[i])
        N.lh.values[i]=N.lh(t.obs[i])
        N.rf.values[i]=N.rf(t.obs[i])
        N.t.values[i]=N.t(t.obs[i])
        
        h.1.values[i]=h_1_function(t.obs[i])
        h.2.values[i]=h_2_function(t.obs[i])
        h.canonical.values[i]=h_canonical_function(t.obs[i])
      }
      
      h.1.conv.N=convolve(N.values, h.1.values)
      h.2.conv.N=convolve(N.values, h.2.values)
      h.canonical.conv.N.t=convolve(N.values, h.canonical.values)
      h.1.conv.N.lf=convolve(N.lf.values, h.1.values)
      h.2.conv.N.lf=convolve(N.lf.values, h.2.values)
      h.1.conv.N.lh=convolve(N.lh.values, h.1.values)
      h.2.conv.N.lh=convolve(N.lh.values, h.2.values)
      h.1.conv.N.rf=convolve(N.rf.values, h.1.values)
      h.2.conv.N.rf=convolve(N.rf.values, h.2.values)
      h.1.conv.N.t=convolve(N.t.values, h.1.values)
      h.2.conv.N.t=convolve(N.t.values, h.2.values)
      
      
      ###################################
      ## 2 Simulation iteraton 
      ##   (for n_simulation = 500 runs)
      ###################################
      
      Y.1.mat=matrix(NA, ncol = 284, nrow = N.participants)
      Y.2.mat=matrix(NA, ncol = 284, nrow = N.participants)
      
      # The following for-chunk generates data.
      for (iter in 1:N.participants) {
        
        epsilon_1=rmvn(121, mu=rep(0, 2), V=V_resting)
        epsilon_2=rmvn(16, mu=rep(0, 2), V=V_task)
        epsilon_3=rmvn(88, mu=rep(0, 2), V=V_resting)
        epsilon_4=rmvn(17, mu=rep(0, 2), V=V_task)
        epsilon_5=rmvn(42, mu=rep(0, 2), V=V_resting)
        epsilon_noise=rbind(epsilon_1, epsilon_2, epsilon_3, epsilon_4, epsilon_5)
        
        Y.1.mat[iter,]=9000+h.1.conv.N+h.1.conv.N.lf+h.1.conv.N.lh+h.1.conv.N.rf+h.1.conv.N.t+epsilon_noise[,1]
        Y.2.mat[iter,]=9000+h.2.conv.N+h.2.conv.N.lf+h.2.conv.N.lh+h.2.conv.N.rf+h.2.conv.N.t+epsilon_noise[,2]
        
      }
      
      est_ptFCE=ptFCE(Y_k=Y.1.mat, Y_l = Y.2.mat, N=N.values, TR=0.72, freq_plot=FALSE)
      Connectivity_estimated_2[node_k, node_l]=est_ptFCE$est
      Connectivity_estimated_2[node_l, node_k]=Connectivity_estimated_2[node_k, node_l]
      
    }
  }
}


par(mfrow=c(1,3))
plot(Connectivity_true, main="True Connectivity", xlab="", ylab="", col = topo.colors)
plot(Connectivity_estimated_1, main="Estimated Connectivity for Data from Mechanism 1",
     xlab="", ylab="", key = NULL, col = topo.colors)
plot(Connectivity_estimated_2, main="Estimated Connectivity for Data from Mechanism 2",
     xlab="", ylab="", key = NULL, col = topo.colors)

