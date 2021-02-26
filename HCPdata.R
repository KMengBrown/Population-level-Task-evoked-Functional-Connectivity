library(oro.nifti)

### Obtain the file paths for all subjects, mask and anatomical image
### Assuming the data folders are saved in a folder called HCP
dat.dir = "~HCP"

### Specify the directory for saving the time courses
dat.out = "~/TimeCourses"

### Paths of subject specific zip foldrs
file.dir = dir(dat.dir, pattern = "*", full.names = TRUE)

### Path to the AAL atlas
file.parc = "~/atlas/AAL2.nii"

### Read in the parcellation and anatomical image
parc = readNIfTI(file.parc)
my.mean = function(x){mean(x, na.rm = TRUE)}

### total number of subjects
n.subj = length(file.dir)

### Extract the unique parcel indicators
parcels = unique(as.vector(parc))
parcels = parcels[order(parcels)][2:length(parcels)]

### The number of parcels
n.clust = length(parcels)

setwd(dat.out)

for(subj in 1:n.subj){
  
  ### Unzip the subject specific zip folder
  unzip(file.dir[subj], exdir = dat.dir.unzip)
  
  ### Extract the file of interest
  dir.fmri = dir(dat.dir.unzip, pattern = "*", full.names = TRUE)
  subject = dir(dat.dir.unzip, pattern = "*", full.names = FALSE)
  file.fmri = paste(dir.fmri,"/MNINonLinear/Results/tfMRI_MOTOR_RL/tfMRI_MOTOR_RL.nii.gz",sep="")
  tm.matrx = c()
  
  ### Read in the fMRI object
  fmri1 = readNIfTI(file.fmri)
  
  ### obtain the parcel specific time courses
  tm = dim(fmri1)[4]
  for(c in 1:n.clust){
    
    ind.clust = which(parc == parcels[c], arr.ind = TRUE)
    time.matrx = c()
    for(t in 1:tm){
      temp.img = fmri1[,,,t]
      time.matrx = cbind(time.matrx, temp.img[ind.clust])
    }
    temp.vec = apply(time.matrx, 2, my.mean)
    tm.matrx = rbind(tm.matrx, temp.vec)
  }
  
  ### Save the resulting matrix of time courses
  setwd(dat.out)
  save(tm.matrx, file = paste(paste("tmmatrx_", subject, sep = ""),".rda", sep = ""))
  
  ### Delete the unzipped folder
  system(paste("rm -r", dir.fmri))
  
  print(subj)
  
}

### Unzip each subject specific folder and save the onset folders
for(subj in 1:n.subj){
  
  unzip(file.dir[subj], exdir = dat.dir.unzip)
  
  dir.fmri = dir(dat.dir.unzip, pattern = "*", full.names = TRUE)
  subject = dir(dat.dir.unzip, pattern = "*", full.names = FALSE)
  file.fmri = paste(dir.fmri,"/MNINonLinear/Results/tfMRI_MOTOR_RL/EVs",sep="")
  
  setwd(dat.out)
  
  system(paste(paste("mv", file.fmri), dat.out))
  
  system(paste("mv EVs", paste(subject, "EVs", sep = "")))
  system(paste("rm -r", dir.fmri))
  
  print(subj)
  
}
