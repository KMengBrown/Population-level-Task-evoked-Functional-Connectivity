###############################################################
##############     HCP Data Set Analysis Using    #############
##############   the AMUSE-ptFCE Algorithm and    #############
##############         Existing Methods           #############
###############################################################

### We present results obtained by the estimation of ptFC using a data set 
### in the task-evoked functional MRI component of HCP.

### The visulization of estimation results here is Figures 6 and 7 in 
### the paper "Population-level Task-evoked Functional Connectivity" by 
### Kun Meng and Ani Eloyan (hereafter ME).

#install.packages("brainGraph")
library("brainGraph")
#install.packages("scatterplot3d")
library("scatterplot3d")
#install.packages("network")
library("network")

###################################
## 1 Data
###################################

## Put the "Subjects_IDs.csv" file and all the 308 ".rda" data files in the same folder 
## and set the corresponding working directory, where the 308 ".rda" files are the output 
## from running the "HCPdata.R" code.

###################################
## 1.1 Read data and orgnize them 
##     into an appliable form
###################################

## Read data
IDs=read.csv("Subjects_IDs.csv")
file_names=vector()
for (i in 1:dim(IDs)[1]) { file_names[i]=paste("tmmatrx_", as.character(IDs[i,1]), ".rda", sep = "") }

Data_time_courses=matrix(NA, ncol = 284, nrow = 1)
for (i in 1:dim(IDs)[1]) {
  print(i)
  load(file = file_names[i])
  Data_time_courses=rbind(Data_time_courses, as.matrix(tm.matrx))
}
Data_time_courses=Data_time_courses[-1,]

## This is a 36960-by-284 matrix.
## There are 308 participants in this data set, and each participant
## has 120 regions of interest in his/her brain. Therefore, there are
## 36960=308*120 rows in this matrix. 
## Each BOLD signal is observed at 284 time points (T=0,1,...,283).
## Therefore, there are 284 columns.
dim(Data_time_courses)


###################################
## 1.2 Organize the entire data 
##     matrix into 120 matrices
###################################

## Each of the 120 matrices corresponds to a region of interest.
Data_list=list()
for (i in 1:120) {
  Data_list[[i]]=matrix(NA, nrow = 1, ncol = 284)
  for (j in 1:308) { Data_list[[i]]=rbind(Data_list[[i]], Data_time_courses[(j-1)*120+i,]) }
  Data_list[[i]]=Data_list[[i]][-1,]
}
# Data_list[[i]][j,] denotes the task-evoked BOLD signals in the 
# i-th region of the j-th participant.


## The 1st matrix Data_list[[1]] corresponds to the seed region.


###################################
## 2 Experimental Parameters
###################################

## Repeatition time
TR=0.72

## Observation time points
t.obs=(1:284)*TR

## The task stimulus signal of interest
N=function(t){
  resp=0
  if(t>=86.5 & t<98.5){ resp=1 }
  if(t>=162 & t<174){ resp=1 }
  return(resp)
}

## Task stimulus signal values at observation time points
N.values=vector()
for (i in 1:284) { N.values[i]=N(t.obs[i]) }



######################################################
## 3 Apply the AMUSE-ptFCE Algorithm
######################################################

## There are 120 regions of interest in the AAL atlas.
iteration.index=1:120
## Because of numerical issues, we omit the regions CB7.L, CB7.R, and CB10.L 
## and investigate the rest 117 regions.
iteration.index=iteration.index[-c(105, 106, 111)]

est_AMUSE_ptFCE=vector()

for (i in iteration.index) {
  print(i)
  results=AMUSE_ptFCE(Y_k = Data_list[[1]], Y_l = Data_list[[i]], N=N.values, TR=TR, freq_plot = FALSE)
  est_AMUSE_ptFCE[i]=results$est
}



######################################################
## 4 Naive Pearson correlation
######################################################

est_naive_Pearson_correlation=vector()

for (i in iteration.index) {
  
  print(i)
  
  Data.mat.1=Data_list[[1]]
  Data.mat.i=Data_list[[i]]
  
  corr.vector=vector()
  for (j in 1:308) {
    corr.vector[j]=cor(Data.mat.1[j,], Data.mat.i[j,])
  }
  est_naive_Pearson_correlation[i]=median(corr.vector)
  
}


######################################################
## 5 Task Pearson correlation
######################################################

time.section.1=which(t.obs>=86.5 & t.obs<98.5)   
time.section.2=which(t.obs>=162 & t.obs<174)
time.section=c(time.section.1, time.section.2)

est_task_Pearson_correlation=vector()
for (i in iteration.index) {
  
  print(i)
  
  Data.mat.1=Data_list[[1]][, time.section]
  Data.mat.i=Data_list[[i]][, time.section]
  
  corr.vector=vector()
  for (j in 1:308) {
    corr.vector[j]=cor(Data.mat.1[j,], Data.mat.i[j,])
  }
  est_task_Pearson_correlation[i]=median(corr.vector)
  
}



######################################################
## 6 Beta-series regression
######################################################

## We apply the procedure described in Chapter 9.2 of Ashby (2019) 
## to estimate functional connectivity using beta-series regression, 
## except we ignore the “nuisance term” therein.

N.participants=308

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
  if(t>=56.244 & t<68.244){ resp=1 }
  if(t>=131.75 & t<143.75){ resp=1 }
  return(resp)
}

N=function(t){ return(N.rh(t)) }


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

est_beta_series=vector()
for (i in iteration.index) {
  
  print(i)
  
  Data.mat.1=Data_list[[1]]
  Data.mat.i=Data_list[[i]]
  
  est.store=vector()
  for (iter in 1:N.participants) {
    Y.k.regression=Data.mat.1[iter,]
    Y.l.regression=Data.mat.i[iter,]
    reg.k=lm(Y.k.regression~DESIGN.MATRIX)
    beta.est.k=as.vector(reg.k$coefficients[2:(dim(design.matrix.N)[2]+1)])
    reg.l=lm(Y.l.regression~DESIGN.MATRIX)
    beta.est.l=as.vector(reg.l$coefficients[2:(dim(design.matrix.N)[2]+1)])
    est.store[iter]=cor(beta.est.k, beta.est.l)
  }
  
  est_beta_series[i]=median(est.store)
  
}



######################################################
## 7 Coherence analysis
######################################################

library("seewave")

est_coherence_analysis=vector()
for (i in iteration.index) {
  
  print(i)
  
  Data.mat.1=Data_list[[1]]
  Data.mat.i=Data_list[[i]]
  
  est.store=vector()
  for (iter in 1:N.participants) {
    coh.results=coh(Data.mat.1[iter, ], Data.mat.i[iter, ], f=1/TR, plot=FALSE)
    est.store[iter]=median(coh.results[coh.results[,1]<0.15,2])
  }
  est_coherence_analysis[i]=median(est.store)
  
}



######################################################
## 8 Visualization in Figure 6
######################################################

iteration.index=1:120
iteration.index=iteration.index[-c(105, 106, 111)]
AAA=iteration.index[which(iteration.index%%2==0)]
BBB=iteration.index[which(iteration.index%%2==1)]
iteration.index=c(BBB, AAA)
iteration.index=iteration.index[-1]

naive_Pearson_correlation_plot=(est_naive_Pearson_correlation[iteration.index]-min(est_naive_Pearson_correlation[iteration.index]))/(max(est_naive_Pearson_correlation[iteration.index])-min(est_naive_Pearson_correlation[iteration.index]))
task_Pearson_correlation_plot=(est_task_Pearson_correlation[iteration.index]-min(est_task_Pearson_correlation[iteration.index]))/(max(est_task_Pearson_correlation[iteration.index])-min(est_task_Pearson_correlation[iteration.index]))
beta_series_plot=(est_beta_series[iteration.index]-min(est_beta_series[iteration.index]))/(max(est_beta_series[iteration.index])-min(est_beta_series[iteration.index]))
coherence_analysis_plot=(est_coherence_analysis[iteration.index]-min(est_coherence_analysis[iteration.index]))/(max(est_coherence_analysis[iteration.index])-min(est_coherence_analysis[iteration.index]))
AMUSE_ptFCE_plot=(est_AMUSE_ptFCE[iteration.index]-min(est_AMUSE_ptFCE[iteration.index]))/(max(est_AMUSE_ptFCE[iteration.index])-min(est_AMUSE_ptFCE[iteration.index]))

par(mfrow=c(1, 1), mar=c(4.5, 4, 0.5, 0.5), xaxt="n")

ddd=0.25
plot((1:116)+ddd, naive_Pearson_correlation_plot, col="red", 
     type = "l", ylim = c(0,1.3), 
     xlab = " ", ylab = "Connectivity", lwd=2,
     xlim = c(1,117))
for (i in 1:116) {
  abline(v=i+ddd, lty=3)
}
lines((1:116)+ddd, naive_Pearson_correlation_plot, col="red", 
      type = "l", lwd=2)
text(x = seq(1, 116, by=1), par("usr")[3]-0.15, labels = aal2.120$name[iteration.index], srt = 90, pos = 1, xpd = TRUE,
     cex=0.7)
lines((1:116)+ddd, task_Pearson_correlation_plot, col="orange", lwd=2)
lines((1:116)+ddd, beta_series_plot, col="darkgreen", lwd=2)
lines((1:116)+ddd, coherence_analysis_plot, col="blue", lwd=2)
lines((1:116)+ddd, AMUSE_ptFCE_plot, lwd=2)

lines((1:116)+ddd, naive_Pearson_correlation_plot, col="red", 
      pch=16, type = "p", cex=1.2,
      xlab = "Region Indices", ylab = "Connectivity", lwd=1)
lines((1:116)+ddd, task_Pearson_correlation_plot, col="orange", lwd=1, pch=16, type = "p", cex=1.2)
lines((1:116)+ddd, beta_series_plot, col="darkgreen", lwd=1, pch=16, type = "p", cex=1.2)
lines((1:116)+ddd, coherence_analysis_plot, col="blue", lwd=1, pch=16, type = "p", cex=1.2)
lines((1:116)+ddd, AMUSE_ptFCE_plot, lwd=1, pch=16, type = "p", cex=1.2)
legend("topright", c("AMUSE-ptFCE algorithm", "Naive Pearson correlation", "Task Pearson correlation", "Beta-series regression", "Coherence analysis"),
       col = c("black", "red", "orange", "darkgreen", "blue"), 
       pch = rep(16, 5))


######################################################
## 8 Visualization in Figure 7
######################################################

iteration.index=1:120
iteration.index=iteration.index[-c(105, 106, 111)]
AAA=iteration.index[which(iteration.index%%2==0)]
BBB=iteration.index[which(iteration.index%%2==1)]
iteration.index=c(BBB, AAA)

RESULTS=cbind(aal2.120$x.mni[iteration.index]-min(aal2.120$x.mni[iteration.index]), aal2.120$y.mni[iteration.index]-min(aal2.120$y.mni[iteration.index]), aal2.120$z.mni[iteration.index]-min(aal2.120$z.mni[iteration.index]))

par(mfrow=c(1,2), mar=rep(0,4))

s3d=scatterplot3d(RESULTS,
                  xlab = " ", ylab = " ", zlab = " ", 
                  pch=20, type = "p", highlight.3d = TRUE,
                  box = FALSE, tick.marks = FALSE, cex.symbols = 3,
                  col.grid="lightblue", 
                  mar = c(0,0,0,0)
)

iteration.index=iteration.index[-1]
plot.scale.para=(est_AMUSE_ptFCE[iteration.index]-min(est_AMUSE_ptFCE[iteration.index]))/(max(est_AMUSE_ptFCE[iteration.index])-min(est_AMUSE_ptFCE[iteration.index]))
colors=as.color(1-plot.scale.para)
for (i in order(plot.scale.para, decreasing = FALSE)) {
  PPP=rbind(RESULTS[1,], RESULTS[i,])
  s3d$points3d(PPP, type = "l", lwd=2*(plot.scale.para), col=colors[i])
}
RN=rep(" ", 117)
RN[1]=aal2.120$name[1]
text(s3d$xyz.convert(RESULTS), labels = RN,
     cex= 1, col = "darkgreen", adj = 1.2)
RN=aal2.120$name[iteration.index]
RN[1]=" "
text(s3d$xyz.convert(RESULTS), labels = RN,
     cex= 0.5, col = "darkgreen", adj = 1.6)


s3d=scatterplot3d(RESULTS,
                  xlab = " ", ylab = " ", zlab = " ", 
                  pch=20, type = "p", highlight.3d = TRUE,
                  box = FALSE, tick.marks = FALSE, cex.symbols = 3,
                  col.grid="lightblue", 
                  mar = c(0,0,0,0)
)
N=30
alg2.plot=order(est_AMUSE_ptFCE[iteration.index], decreasing = TRUE)[1:N]
for (i in alg2.plot) {
  PPP=rbind(RESULTS[1,], RESULTS[i,])
  s3d$points3d(PPP, type = "l", col="blue", lwd=2)
}

names.store=aal2.120$name[iteration.index]

RN=rep(" ", 117)
RN[1]=aal2.120$name[1]
text(s3d$xyz.convert(RESULTS), labels = RN,
     cex= 1, col = "darkgreen", adj = 1.2)

RN=rep(" " , 117)
for (i in alg2.plot) {
  RN[i]=names.store[i]
}
RN[1]=" "
text(s3d$xyz.convert(RESULTS), labels = RN,
     cex= 0.5, col = "darkgreen", adj = 1.6)


############################
## End
############################
