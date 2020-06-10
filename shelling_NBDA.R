#######################################################################################################
######################################  R Code to Wild et al (2020):   ################################
############################################ Current Biology ##########################################
### Integrating genetic, environmental and social networks to reveal transmission pathways of a dolphin foraging innovation ###

## load NBDA package available at https://rdrr.io/github/whoppitt/NBDA/

install.packages("devtools")
library(devtools)
install_github("whoppitt/NBDA")
library(NBDA)

## read in empirical data (shelling)

# set working directory
setwd("..../data")

## save tabs from excel file Data S1 as separate csv files, then load
## (available in the supplementary electronic material) 

# read association network
SRI_all <- read.csv("Association Matrix.csv", row.names=1, header=TRUE)
SRI_all <- as.matrix(SRI_all)

# read ecological network
ecol_all <- read.csv("Home Range Overlaps.csv", row.names=1, header=TRUE)
ecol_all <- as.matrix(ecol_all)  

# read relatedness network
relate_all <- read.csv("Relatedness Matrix.csv", row.names=1, header=TRUE)
relate_all <- as.matrix(relate_all)

# read ILVs
ILV_all <- read.csv("Additional File 5 - ILVs.csv", header=TRUE, sep=",")

ILV_all[ILV_all=="<NA>"]=NA

# simulations revealed that power of NBDA to detect learning is highest at a cut-off of 11 sightings (see STAR methods)
# get list of IDs that have been seen at least 11 times
ILV <- subset(ILV_all, subset=ILV_all$Number_sightings>10)

IDs <- ILV$id_individual
length(IDs) # 310 individuals remaining

num <- which(colnames(SRI_all) %in% IDs) # extract positions of the 310 individuals in the networks
SRI <- SRI_all[num, num] # reduce association network to 310 individuals
dim(SRI) # ensure dimensions are 310x310
class(SRI) # ensure it is a matrix

# repeat for home range overlaps
num <- which(colnames(ecol_all) %in% IDs)
ecol <- ecol_all[num, num]
dim(ecol)
class(ecol)

# repeat for relatedness matrix
num <- which(colnames(relate_all) %in% IDs)
relate <- relate_all[num, num]
dim(ecol)
class(relate)


##### Order of acquisition as vector: Three-letter codes of individuals that have acquired shelling (in order of observation)
shellers <- c("WIM", "PAR", "GVY", "DEN", "JUL", "MSH", "HLM", "ZED", "BNN", "JON", "MOY", "ARN", "HLI", "DET", "GRE")

order <- NULL # create an object to store the vector of acquistion

for (i in 1:length(shellers)){ # for each sheller, extract the position in the networks and ILV data frame
  order[i] <- which(IDs==shellers[i])
}

order <- as.vector(order)
OAc <- order # order of acquisition is stored as a vector with the position of informed individuals in the network(s).

## prepare individual level variables
Sex <- ILV$Sex_1_0 # sex as 0.5 for males, -0.5 for females and 0 for unknown sex
Number_of_sightings <- ILV$Number_sightings-mean(ILV$Number_sightings)
Av_water_depth <- ILV$Av_water_depth-mean(ILV$Av_water_depth) # water depth averaged across each individual's sightings
Av_group_size <- ILV$Av_group_size-mean(ILV$Av_group_size) # an individual's average number of group members

## Here, we want to set E as the baseline level of the factor (i.e. all zeroes = Haplotype E).
# see STAR methods for details
HaplotypeH <- (ILV$Haplotype=="H")*1
HaplotypeH[is.na(HaplotypeH)]<-0
HaplotypeE <- (ILV$Haplotype=="E")*1
HaplotypeE[is.na(HaplotypeE)]<-0
HaplotypeD <- (ILV$Haplotype=="D")*1
HaplotypeD[is.na(HaplotypeD)]<-0
HaplotypeF <- (ILV$Haplotype=="F")*1
HaplotypeF[is.na(HaplotypeF)]<-0
HaplotypeNotED <- 1-HaplotypeE-HaplotypeD


###################### create NBDA object with CONSTANT ILVs and MULTIPLE s parameters (estimated individually for each network):

n.assMatrix <- 3 # number of matrices
assMatrix3 <- array(data = c(SRI,ecol,relate), dim=c(nrow(SRI), ncol(SRI), n.assMatrix)) # create an array with the three matrices
Sex <- matrix(data = Sex, nrow=length(IDs), byrow=F) # all ILVs need to go into a matrix
Number_of_sightings <- matrix(data = Number_of_sightings, nrow=length(IDs), byrow=F) 
Av_water_depth <- matrix(data = Av_water_depth, nrow=length(IDs), byrow=F)
Av_group_size <- matrix(data = Av_group_size, nrow=length(IDs), byrow=F)
HaplotypeH <- matrix(data = HaplotypeH, nrow=length(IDs), byrow=F)
HaplotypeE <- matrix(data = HaplotypeE, nrow=length(IDs), byrow=F)
HaplotypeD <- matrix(data = HaplotypeD, nrow=length(IDs), byrow=F)
HaplotypeF <- matrix(data = HaplotypeF, nrow=length(IDs), byrow=F)
HaplotypeNotED<- matrix(data = HaplotypeNotED, nrow=length(IDs), byrow=F)

# ILVs need to be combined into a character vector to create the NBDA data object
ILV <- c("Sex","Number_of_sightings","Av_water_depth","Av_group_size","HaplotypeH","HaplotypeNotED","HaplotypeD","HaplotypeF") 

# label for NBDA object
label <- "shelling1redhap"

# create the NBDA Data Object 
# since we are using the unconstrained version (where ILVs are allowed to separately influence social and asocial effects)
# we are using the unconstrained models, where ILVs can influence social and asocial learning independently
# we only define asoc_ilv and int_ilv, but do not define multi_ilv
# multiplicative models (multi_ilv) allow ILVs to influence both social and asocial learning
# but are assumed to do so to the same extent
nbdaDataSHELLING_H2asBaseline <- nbdaData(label=label, 
                                          assMatrix=assMatrix3, # array with association matrices
                                          asoc_ilv=ILV, # influencing asocial learning. refers to the character vector with the ILV matrices
                                          int_ilv=ILV, # influencing social learning. 
                                          multi_ilv="ILVabsent", # we are not fitting multiplicative models
                                          orderAcq=OAc, # individual three-letter codes in order of acquisition
                                          asocialTreatment="constant") # creates OADA object

# in the NBDADataObject only asoc_ilv and int_ilv are defined, multi_ilv are "absent"
nbdaDataSHELLING_H2asBaseline@asoc_ilv 
nbdaDataSHELLING_H2asBaseline@int_ilv
nbdaDataSHELLING_H2asBaseline@multi_ilv

print(nbdaDataSHELLING_H2asBaseline)


##########################################################################################################
#### prepare a matrix that specifies all possible combinations of networks and ILVs:
# each row determines one model. The first three position in each row correspond to the networks (assoc, ecol, relate)
# the following 8 positions to the asoc_ilv and the last 8 positions to the int_ilv
# entries of !=0 mean that the parameter is estimated, while 0 means it isn't


constraintsVectMatrix<-NULL

netVect<-c(1,0,0)
for(a in 0:1){
  for(b in 0:1){
    for(c in 0:1){
      for(d in 0:1){
          for(a2 in 0:1){
            for(b2 in 0:1){
              for(c2 in 0:1){
                for(d2 in 0:1){
                  #ASOCIAL AND SOCIAL EFFECTS PRESENT
                  #Full version of Haplotype
                  constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,1,1,1,1))
                  #F dropped from Haplotype
                  constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,1,1,1,0))
                  #F and H dropped from Haplotype
                  constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,1,1,0))

                  #ASOCIAL EFFECTS PRESENT
                  #Full version of Haplotype
                  constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,0,0,0,0))
                  #F dropped from Haplotype
                  constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,0,0,0,0))
                  #F and H dropped from Haplotype
                  constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,0,0,0))
                  
                  #SOCIAL EFFECTS PRESENT
                  #Full version of Haplotype
                  constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,1))
                  #F dropped from Haplotype
                  constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,0))
                  #F and H dropped from Haplotype
                  constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,1,1,0))
                  
                  #NO EFFECT
                  #HaplotypeAbsent
                  constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,0,0,0))
                }
              }
            }
          }
      }
    }
  }
}

netVect<-c(0,1,0)
for(a in 0:1){
  for(b in 0:1){
    for(c in 0:1){
      for(d in 0:1){
        for(a2 in 0:1){
          for(b2 in 0:1){
            for(c2 in 0:1){
              for(d2 in 0:1){
                #ASOCIAL AND SOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,1,1,1,1))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,1,1,1,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,1,1,0))
                
                #ASOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,0,0,0,0))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,0,0,0,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,0,0,0))
                
                #SOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,1))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,1,1,0))
                
                #NO EFFECT
                #HaplotypeAbsent
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,0,0,0))
              }
            }
          }
        }
      }
    }
  }
}

netVect<-c(0,0,1)
for(a in 0:1){
  for(b in 0:1){
    for(c in 0:1){
      for(d in 0:1){
        for(a2 in 0:1){
          for(b2 in 0:1){
            for(c2 in 0:1){
              for(d2 in 0:1){
                #ASOCIAL AND SOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,1,1,1,1))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,1,1,1,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,1,1,0))
                
                #ASOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,0,0,0,0))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,0,0,0,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,0,0,0))
                
                #SOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,1))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,1,1,0))
                
                #NO EFFECT
                #HaplotypeAbsent
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,0,0,0))
              }
            }
          }
        }
      }
    }
  }
}

netVect<-c(1,0,1)
for(a in 0:1){
  for(b in 0:1){
    for(c in 0:1){
      for(d in 0:1){
        for(a2 in 0:1){
          for(b2 in 0:1){
            for(c2 in 0:1){
              for(d2 in 0:1){
                #ASOCIAL AND SOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,1,1,1,1))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,1,1,1,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,1,1,0))
                
                #ASOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,0,0,0,0))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,0,0,0,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,0,0,0))
                
                #SOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,1))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,1,1,0))
                
                #NO EFFECT
                #HaplotypeAbsent
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,0,0,0))
              }
            }
          }
        }
      }
    }
  }
}

netVect<-c(0,1,1)
for(a in 0:1){
  for(b in 0:1){
    for(c in 0:1){
      for(d in 0:1){
        for(a2 in 0:1){
          for(b2 in 0:1){
            for(c2 in 0:1){
              for(d2 in 0:1){
                #ASOCIAL AND SOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,1,1,1,1))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,1,1,1,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,1,1,0))
                
                #ASOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,0,0,0,0))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,0,0,0,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,0,0,0))
                
                #SOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,1))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,1,1,0))
                
                #NO EFFECT
                #HaplotypeAbsent
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,0,0,0))
              }
            }
          }
        }
      }
    }
  }
}

netVect<-c(1,1,0)
for(a in 0:1){
  for(b in 0:1){
    for(c in 0:1){
      for(d in 0:1){
        for(a2 in 0:1){
          for(b2 in 0:1){
            for(c2 in 0:1){
              for(d2 in 0:1){
                #ASOCIAL AND SOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,1,1,1,1))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,1,1,1,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,1,1,0))
                
                #ASOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,0,0,0,0))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,0,0,0,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,0,0,0))
                
                #SOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,1))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,1,1,0))
                
                #NO EFFECT
                #HaplotypeAbsent
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,0,0,0))
              }
            }
          }
        }
      }
    }
  }
}

netVect<-c(1,1,1)
for(a in 0:1){
  for(b in 0:1){
    for(c in 0:1){
      for(d in 0:1){
        for(a2 in 0:1){
          for(b2 in 0:1){
            for(c2 in 0:1){
              for(d2 in 0:1){
                #ASOCIAL AND SOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,1,1,1,1))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,1,1,1,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,1,1,0))
                
                #ASOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,0,0,0,0))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,0,0,0,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,0,0,0))
                
                #SOCIAL EFFECTS PRESENT
                #Full version of Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,1))
                #F dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,1,1,1,0))
                #F and H dropped from Haplotype
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,1,1,0))
                
                #NO EFFECT
                #HaplotypeAbsent
                constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,0,0,0))
              }
            }
          }
        }
      }
    }
  }
}

netVect<-c(0,0,0)
for(a in 0:1){
  for(b in 0:1){
    for(c in 0:1){
      for(d in 0:1){
        #ASOCIAL EFFECTS PRESENT
        #Full version of Haplotype
        constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,0,0,0,0))
        #F dropped from Haplotype
        constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,0,0,0,0))
        #F and H dropped from Haplotype
        constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,0,0,0))

        #NO EFFECT
        #HaplotypeAbsent
        constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,0,0,0,a2,b2,c2,d2,0,0,0,0))
      }
    }
  }
}


# creating the cumSum for each row within the matrix: we receive a matrix where each parameter is independlty estimated
# by using same numbers (e.g. 1 for both the social and ecological network), their effects are constrained to be the same
constraintsVectMatrix<-t(apply(constraintsVectMatrix,1,cumsum))*constraintsVectMatrix


# running the AIC table function runs OADA across all defined models and calculates model support (AICc weights)
# a progress bar informs what % of models have run
tableSHELLING<-oadaAICtable(nbdadata=nbdaDataSHELLING_H2asBaseline, constraintsVectMatrix=constraintsVectMatrix,writeProgressFile = T)

# save the output object
save(tableSHELLING, file="shelling.AIC_table_indepInts_RedHapinOneCORRECTED.Rdata")
load("shelling.AIC_table_indepInts_RedHapinOneCORRECTED.Rdata")

# print summary results
print(tableSHELLING)
# if there are unfitted models, several columns will show NA (unfitted models are removed below)

# write raw AIC table
write.csv(as.data.frame(tableSHELLING@printTable), "shelling.AIC_table_indepInts_RedHapinOneCORRECTED.csv")


# Create a new object with a printTable that excludes unfitted models
newTableSHELLING<-tableSHELLING
newTableSHELLING@printTable<-tableSHELLING@printTable[!is.nan(tableSHELLING@printTable$aicc)&!is.na(tableSHELLING@printTable$aicc),]

# recalculate model support only including models that were fitted
newTableSHELLING@aicc<-tableSHELLING@aicc[!is.nan(tableSHELLING@aicc)&!is.na(tableSHELLING@aicc)]
newTableSHELLING@MLEs<-tableSHELLING@MLEs[!is.nan(tableSHELLING@aicc)&!is.na(tableSHELLING@aicc),]
newTableSHELLING@MLEilv<-tableSHELLING@MLEilv[!is.nan(tableSHELLING@aicc)&!is.na(tableSHELLING@aicc),]
newTableSHELLING@MLEint<-tableSHELLING@MLEint[!is.nan(tableSHELLING@aicc)&!is.na(tableSHELLING@aicc),]

newTableSHELLING@printTable<-newTableSHELLING@printTable[order(newTableSHELLING@printTable$aicc),]
newTableSHELLING@printTable$deltaAICc<-newTableSHELLING@printTable$aicc-newTableSHELLING@printTable$aicc[1]

newTableSHELLING@printTable$RelSupport<- exp(-0.5*newTableSHELLING@printTable$deltaAICc)
newTableSHELLING@printTable$AkaikeWeight<-newTableSHELLING@printTable$RelSupport/sum(newTableSHELLING@printTable$RelSupport)


newTableSHELLING@deltaAIC<-newTableSHELLING@aicc-min(newTableSHELLING@aicc)

# calculate model support and akaike weights for each model
newTableSHELLING@RelSupport<- exp(-0.5*newTableSHELLING@deltaAIC)
newTableSHELLING@AkaikeWeight<-newTableSHELLING@RelSupport/sum(newTableSHELLING@RelSupport)


# extract the number of unfitted models that were removed. 
dim(tableSHELLING@printTable)[1]-dim(newTableSHELLING@printTable)[1]
dim(tableSHELLING@printTable)[1]
## In our data set, 902 models could not be fitted put of 17'984
# they likely have too many parameters for the data set 

# extract support for each network combination (order of networks is association, ecology, relatedness)
networksSupport_shelling<-networksSupport(newTableSHELLING)
networksSupport_shelling
write.csv(networksSupport_shelling, file="networksSupport_shelling_int_sep.csv")
#75.4% support for transmission along the social network alone

# extract support for each variable
variable_support <- variableSupport(newTableSHELLING, includeAsocial = T)
variable_support
write.csv(variable_support, file="variable_support_shelling_int_sep.csv")

# extract weighted medians for each variable (model averaging)
MLE_med  <- modelAverageEstimates(newTableSHELLING,averageType = "median")
MLE_med
write.csv(MLE_med,file="modelWeightedMediansSHELLING.csv")

########################################################################################
# Here, we examine the haplotpe factor (to decide how many levels the haplotype factor should have)
#Full version of Haplotype
constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,1,a2,b2,c2,d2,1,1,1,1))
#F dropped from Haplotype
constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,1,1,1,0,a2,b2,c2,d2,1,1,1,0))
#F and H dropped from Haplotype
constraintsVectMatrix<-rbind(constraintsVectMatrix,c(netVect,a,b,c,d,0,1,1,0,a2,b2,c2,d2,0,1,1,0))

## The support for HaplotypeF gives the support for the full factor 0.44%
## The support for HaplotypeH - Haplotype F gives the support H dropped = 
100*(0.0845914-0.004425819)
#[1] 8.016558 8.01%
## The support for HaplotypeNotED or Haplotype D - support for HaplotypeF is the support for F and H dropped=
100*(0.9730177-0.0845914)
#[1] 88.84263 88.84%

## Strongly supports the use of the fully reduced haplotype factor 
## cut down to the models using this factor

## Create a new object with a printTable that excludes unfitted models and treats haplotype as a 3 level factor (hap E, H and other)
# This cuts the number of models down from <17000 to the <7000 that the results in the manuscript are based on
newTableSHELLING_HapRed<-newTableSHELLING
newTableSHELLING_HapRed@printTable<-newTableSHELLING_HapRed@printTable[(newTableSHELLING@printTable$CONS.ASOC.HaplotypeH==0)&(newTableSHELLING@printTable$CONS.ASOC.HaplotypeF==0)&(newTableSHELLING@printTable$CONS.SOCIAL.HaplotypeH==0)&(newTableSHELLING@printTable$CONS.SOCIAL.HaplotypeF==0),]


newTableSHELLING_HapRed@aicc<-newTableSHELLING_HapRed@aicc[(newTableSHELLING_HapRed@MLEilv[,5]==0)&(newTableSHELLING_HapRed@MLEilv[,8]==0)&(newTableSHELLING_HapRed@MLEint[,5]==0)&(newTableSHELLING_HapRed@MLEint[,8]==0)]
newTableSHELLING_HapRed@MLEs<-newTableSHELLING_HapRed@MLEs[(newTableSHELLING_HapRed@MLEilv[,5]==0)&(newTableSHELLING_HapRed@MLEilv[,8]==0)&(newTableSHELLING_HapRed@MLEint[,5]==0)&(newTableSHELLING_HapRed@MLEint[,8]==0),]
newTableSHELLING_HapRed@MLEilv<-newTableSHELLING_HapRed@MLEilv[(newTableSHELLING_HapRed@MLEilv[,5]==0)&(newTableSHELLING_HapRed@MLEilv[,8]==0)&(newTableSHELLING_HapRed@MLEint[,5]==0)&(newTableSHELLING_HapRed@MLEint[,8]==0),]
newTableSHELLING_HapRed@MLEint<-newTableSHELLING_HapRed@MLEint[(newTableSHELLING@MLEilv[,5]==0)&(newTableSHELLING@MLEilv[,8]==0)&(newTableSHELLING@MLEint[,5]==0)&(newTableSHELLING@MLEint[,8]==0),]

newTableSHELLING_HapRed@printTable$deltaAICc<-newTableSHELLING_HapRed@printTable$aicc-newTableSHELLING_HapRed@printTable$aicc[1]

newTableSHELLING_HapRed@printTable$RelSupport<- exp(-0.5*newTableSHELLING_HapRed@printTable$deltaAICc)
newTableSHELLING_HapRed@printTable$AkaikeWeight<-newTableSHELLING_HapRed@printTable$RelSupport/sum(newTableSHELLING_HapRed@printTable$RelSupport)

newTableSHELLING_HapRed@deltaAIC<-newTableSHELLING_HapRed@aicc-min(na.omit(newTableSHELLING_HapRed@aicc))

newTableSHELLING_HapRed@RelSupport<- exp(-0.5*newTableSHELLING_HapRed@deltaAIC)
newTableSHELLING_HapRed@AkaikeWeight<-newTableSHELLING_HapRed@RelSupport/sum(newTableSHELLING_HapRed@RelSupport)

networksSupport_HapRed<-networksSupport(newTableSHELLING_HapRed)
networksSupport_HapRed
write.csv(networksSupport_HapRed, file="networksSupport_shelling_int_sep_HapRed.csv")

#77.1% support for social network only again

variable_support_HapRed <- variableSupport(newTableSHELLING_HapRed, includeAsocial = T)
variable_support_HapRed
write.csv(variable_support, file="variable_support_shelling_int_sep_HapRed.csv")


MLE  <- modelAverageEstimates(newTableSHELLING_HapRed)
MLE  <- modelAverageEstimates(newTableSHELLING_HapRed, averageType="median")
write.csv(MLE, "MLE_shelling.csv")

#######################################################################################################################################
#Getting 95% confidence intervals using profile likelihood techniques
#This is vital for s parameters since CIs based on SEs will be highly misleading due to frequent assymetry in the profile likelihood 
#######################################################################################################################################
# extract the best model [16]
print(newTableSHELLING_HapRed)[1:10,]

# create object of best model
bestModelData<-constrainedNBDAdata(nbdadata=nbdaDataSHELLING_H2asBaseline,constraintsVect =constraintsVectMatrix[16,])
# run OADA on best model only
model.best.social<-oadaFit(bestModelData)
# check output parameters (1: social network, 2: haplotyeNOT.ED, 3: HaplotypeD, 4: AvGroupSize)
model.best.social@outputPar
#[1]  15.6793765 -21.1332499   3.8529092  -0.6906619
model.best.social@optimisation
model.best.social@aicc

# plot profile likelihoods 
# (which=1 plots for the social learning parameter for social transmission along the network)
# repeat for which=2-4 to extract profile likelihoods for all other parameters too
plotProfLik(which=1,model=model.best.social,range=c(0,150), resolution=20)
plotProfLik(which=1,model=model.best.social,range=c(0,20), resolution=50) # adjust range where plot crosses dashed line (lower)
profLikCI(which=1,model=model.best.social,upperRange=c(120,150),lowerRange=c(0,15)) 

#Lower CI   Upper CI 
#2.058049 144.933292  

# extract how many individuals have learned shelling through social learning rather than independent learning
prop.solve.social.byevent <- oadaPropSolveByST.byevent(nbdadata = model.best.social, model=model.best.social)
prop.solve.social <- oadaPropSolveByST(nbdadata = model.best.social, model=model.best.social)
# overall, 56.7% of shellers have learned socially

#To get the estimates for the lower bound we should really find the corresponding value of the other parameters to plug in when s1 is constrained to this value
bestModelDataS1LowerBound<-constrainedNBDAdata(nbdadata=nbdaDataSHELLING_H2asBaseline,constraintsVect =constraintsVectMatrix[16,],offset=c(2.058049,rep(0,18)))
bestModelS1LowerBound<-oadaFit(bestModelDataS1LowerBound,type="asocial")
bestModelS1LowerBound@outputPar
#Now plug into the prop solve function in one of these two ways:
prop.solve.social.lower <- oadaPropSolveByST(par= c(2.058049, bestModelS1LowerBound@outputPar),model=NULL, nbdadata = bestModelData)
prop.solve.social.lower <- oadaPropSolveByST(model=bestModelS1LowerBound, nbdadata = bestModelDataS1LowerBound)
# lower bound for % of individuals learned socially is at 41.2%

#Same for upper limit
bestModelDataS1UpperBound<-constrainedNBDAdata(nbdadata=nbdaDataSHELLING_H2asBaseline,constraintsVect =constraintsVectMatrix[16,],offset=c(144.933292,rep(0,18)))
bestModelS1UpperBound<-oadaFit(bestModelDataS1UpperBound,type="asocial")
bestModelS1UpperBound@outputPar
prop.solve.social.upper <- oadaPropSolveByST(par= c(144.933292, bestModelS1UpperBound@outputPar),model=NULL, nbdadata = bestModelData)
#or
prop.solve.social.upper <- oadaPropSolveByST(model=bestModelS1UpperBound, nbdadata = bestModelDataS1UpperBound)
# lower bound for % of individuals learned socially is at 73.8%

prop.solve.social
prop.solve.social.lower
prop.solve.social.upper
#57% (95% C.I.= 41-74) (horizontal) social transmission (excluding original innovator)

#for Hap not ED
## essentially the asocial learning rate for not ED is estimated at 0. We can get an upper limit on this:
plotProfLik(which=2,model=model.best.social,range=c(-10,0))
profLikCI(which=2,model=model.best.social,upperRange=c(-2,0))
#Lower CI  Upper CI 
#NA -1.132578 
# backtransformed
exp(1.132578)
#[1] 3.103647 Haplotype E is estimated to be at least 3x faster to learn asocially than ABIKFH

#for Hap D relative to haplotype E now (since hapE is treated as baseline)
plotProfLik(which=3,model=model.best.social,range=c(0,7),resolution=20)
profLikCI(which=3,model=model.best.social,upperRange=c(5,6),lowerRange=c(1,2))
#Lower CI Upper CI 
#1.726526 5.714828
## backtransformed
exp(c(model.best.social@outputPar[3],1.726526, 5.714828))
#[1]  47.129384   5.621092 303.332026
## hap D estimated 47x faster to learn (95% C.I.= 5.6-303) compared to hap E

#for effect of group size on s
rbind(model.best.social@varNames,model.best.social@outputPar)
plotProfLik(which=4,model=model.best.social,range=c(-2,0),resolution=20)
profLikCI(which=4,model=model.best.social,upperRange=c(-0.5,0),lowerRange=c(-1.5,-1))
#Lower CI   Upper CI 
#-1.3645667 -0.1640229 
## backtransformed
exp(c(-model.best.social@outputPar[4],1.3645667 ,0.1640229 ))
#[1] 1.995038 3.914027 1.178241
## for every 1 individual decrease in average group size, s increases by 2x (95% C.I.= 1.18-3.91)

## let us see if there is evidence for s1>0 when we include and thus control for the ecological network- find best model with both

# extract best model with transmission along the social and ecological network 
print(newTableSHELLING_HapRed)[newTableSHELLING_HapRed@printTable$netCombo=="1:2:0"|newTableSHELLING_HapRed@printTable$netCombo=="1:2:3",]
best12ModelData<-constrainedNBDAdata(nbdadata=nbdaDataSHELLING_H2asBaseline,constraintsVect =constraintsVectMatrix[12816,])
best12Model<-oadaFit(best12ModelData)
best12Model@outputPar
#[1]  15.6794237   0.0000000 -20.6526418   3.8528746  -0.6906612
#The same but with s2 (ecol) estimated at 0
best12Model@optimisation
best12Model@aicc

plotProfLik(which=1,model=best12Model,range=c(0,150), resolution=20)
profLikCI(which=1,model=best12Model,upperRange=c(120,150),lowerRange=c(0,15))
#Lower CI  Upper CI 
#2.05803 144.93329  
## Has no influence on the fit- just estimated at zero whatever so the bit below just gives the same answers as without the ecological network

## repeat the same for social+relatedness netwok
# extract best model with transmission along the social and relatedness network 
print(newTableSHELLING_HapRed)[newTableSHELLING_HapRed@printTable$netCombo=="1:0:2"|newTableSHELLING_HapRed@printTable$netCombo=="1:2:3",]
best12ModelData<-constrainedNBDAdata(nbdadata=nbdaDataSHELLING_H2asBaseline,constraintsVect =constraintsVectMatrix[7696,])
best12Model<-oadaFit(best12ModelData)
best12Model@outputPar
#[1]  15.6794237   0.0000000 -20.6526418   3.8528746  -0.6906612
#The same but with s3 (relatedness) estimated at 0
best12Model@optimisation
best12Model@aicc

plotProfLik(which=1,model=best12Model,range=c(0,150), resolution=20)
profLikCI(which=1,model=best12Model,upperRange=c(120,150),lowerRange=c(0,15))
#Lower CI  Upper CI 
#2.05803 144.93329  
## Also has no influence on the fit- just estimated at zero whatever so the bit below just gives the same answers as without the relatedness network

## to conclude, social learning along the association network is supported.
## the effect size remains the same even in presence of the ecolocial or relatedness network
