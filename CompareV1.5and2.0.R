cMat2 <- read.csv(
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/commMat2.0.csv',
          header = TRUE)
rNames <- cMat2[, 'X']
cMat2 <- as.matrix(cMat2[, -1])
cMat2 <- apply(cMat2, 2, as.numeric)
rownames(cMat2) <- rNames

cMat1 <- read.csv(
  'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/commMat.csv',
  header = TRUE)
rNames <- cMat1[, 'X']
cMat1 <- as.matrix(cMat1[, -1])
cMat1 <- apply(cMat1, 2, as.numeric)
rownames(cMat1) <- rNames

# doing this messes up the entries
cMat3 <- cMat2[c("Winds", "Surface Temperature", "Surface Salinity", "Bottom Temperature", "Bottom Salinity", "Habitat: Pelagic", "Recreational Groundfish Fishery", "Commercial Groundfish Fishery", 
                 "Commercial Pelagic Fishery", "Commercial Shellfish Fishery", "Source Water Proportions", "Habitat: Nearshore", "Stratification", "Habitat: Seafloor & Demersal", "Forage Fish", 
                 "Protected Species", "Primary Production", "Groundfish", "Fished Invertebrates", "Tidal Forcing", "Air Temperature", "Precipitation", "Mid Atlantic Groundfish", 
                 "Cultural Practices & Attachments", "Seafood", "Employment", "Profits", "Copepods & Micronekton", "Benthos", "Detritus & Bacteria", "Gelatinous Zooplankton"),
               c("Winds", "Surface.Temperature", "Surface.Salinity", "Bottom.Temperature", "Bottom.Salinity", "Habitat..Pelagic", "Recreational.Groundfish.Fishery", "Commercial.Groundfish.Fishery",
                        "Commercial.Pelagic.Fishery", "Commercial.Shellfish.Fishery", "Source.Water.Proportions", "Habitat..Nearshore", "Stratification", "Habitat..Seafloor...Demersal", "Forage.Fish", 
                        "Protected.Species", "Primary.Production", "Groundfish", "Fished.Invertebrates", "Tidal.Forcing", "Air.Temperature", "Precipitation", "Mid.Atlantic.Groundfish", 
                        "Cultural.Practices...Attachments", "Seafood", "Employment", "Profits", "Copepods...Micronekton", "Benthos", "Detritus...Bacteria", "Gelatinous.Zooplankton")]
all(rownames(cMat1)==rownames(cMat2))
all(colnames(cMat1)==colnames(cMat2))
all(cMat1==cMat2)
test1 <- cMat1!=cMat2

cMat1[c("Recreational Groundfish Fishery", "Commercial Groundfish Fishery", "Commercial Pelagic Fishery", "Commercial Shellfish Fishery"),
      c("Cultural.Practices...Attachments", "Seafood", "Employment", "Profits")] <- 0
det(-cMat1)
cMat1["Primary Production", c("Habitat..Nearshore", "Habitat..Seafloor...Demersal")] <- 0
det(-cMat1)
test2 <- cMat1==cMat2
all(test2)
det(-cMat2)

testOrder1 <- commMat[c("Winds", "Surface Temperature", "Surface Salinity", "Bottom Temperature", "Bottom Salinity", "Habitat: Pelagic", "Recreational Groundfish Fishery", "Commercial Groundfish Fishery", 
                 "Commercial Pelagic Fishery", "Commercial Shellfish Fishery", "Source Water Proportions", "Habitat: Nearshore", "Stratification", "Habitat: Seafloor & Demersal", "Forage Fish", 
                 "Protected Species", "Primary Production", "Groundfish", "Fished Invertebrates", "Tidal Forcing", "Air Temperature", "Precipitation", "Mid Atlantic Groundfish", 
                 "Cultural Practices & Attachments", "Seafood", "Employment", "Profits", "Copepods & Micronekton", "Benthos", "Detritus & Bacteria", "Gelatinous Zooplankton"),
               c("Winds", "Surface.Temperature", "Surface.Salinity", "Bottom.Temperature", "Bottom.Salinity", "Habitat..Pelagic", "Recreational.Groundfish.Fishery", "Commercial.Groundfish.Fishery",
                 "Commercial.Pelagic.Fishery", "Commercial.Shellfish.Fishery", "Source.Water.Proportions", "Habitat..Nearshore", "Stratification", "Habitat..Seafloor...Demersal", "Forage.Fish", 
                 "Protected.Species", "Primary.Production", "Groundfish", "Fished.Invertebrates", "Tidal.Forcing", "Air.Temperature", "Precipitation", "Mid.Atlantic.Groundfish", 
                 "Cultural.Practices...Attachments", "Seafood", "Employment", "Profits", "Copepods...Micronekton", "Benthos", "Detritus...Bacteria", "Gelatinous.Zooplankton")]
all(testOrder1==cMat1)
det(-testOrder1)
det(-commMat)
testOrder2 <- commMat[c("Winds", "Surface Temperature", "Surface Salinity", "Bottom Temperature", "Bottom Salinity", "Habitat: Pelagic", "Recreational Groundfish Fishery", "Commercial Groundfish Fishery", 
                        "Commercial Pelagic Fishery", "Commercial Shellfish Fishery", "Source Water Proportions", "Habitat: Nearshore", "Stratification", "Habitat: Seafloor & Demersal", "Forage Fish", 
                        "Protected Species", "Primary Production", "Groundfish", "Fished Invertebrates", "Tidal Forcing", "Air Temperature", "Precipitation", "Mid Atlantic Groundfish", 
                        "Cultural Practices & Attachments", "Seafood", "Employment", "Profits", "Copepods & Micronekton", "Benthos", "Detritus & Bacteria", "Gelatinous Zooplankton"),]
testOrder2 <- testOrder2[,c("Winds", "Surface.Temperature", "Surface.Salinity", "Bottom.Temperature", "Bottom.Salinity", "Habitat..Pelagic", "Recreational.Groundfish.Fishery", "Commercial.Groundfish.Fishery",
                        "Commercial.Pelagic.Fishery", "Commercial.Shellfish.Fishery", "Source.Water.Proportions", "Habitat..Nearshore", "Stratification", "Habitat..Seafloor...Demersal", "Forage.Fish", 
                        "Protected.Species", "Primary.Production", "Groundfish", "Fished.Invertebrates", "Tidal.Forcing", "Air.Temperature", "Precipitation", "Mid.Atlantic.Groundfish", 
                        "Cultural.Practices...Attachments", "Seafood", "Employment", "Profits", "Copepods...Micronekton", "Benthos", "Detritus...Bacteria", "Gelatinous.Zooplankton")]
all(testOrder2==cMat1)
det(-testOrder2)
det(-commMat)
all(testOrder2==testOrder1)
write.csv(testOrder1, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/testOrder.csv')

A <- c(1,-1,-1,0)
B <- c(0,-1,1,1)
C <- c(0,0,0,1)
D <- c(-1,1,0,0)

ABCD <- cbind(A,B,C,D)
det(-ABCD)

BACD <- cbind(B,A,C,D)
det(-BACD)

det(-rbind(BACD[2,],BACD[1,],BACD[3,],BACD[4,])) # same as first

switchAB <- diag(4)
switchAB[1,1] <- 0
switchAB[1,2] <- 1
switchAB[2,1] <- 1
switchAB[2,2] <- 0
switchAB%*%ABCD

cMat2 <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/commMat2.0.csv',
                  header = TRUE)
cMattest <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/testOrder.csv',
                  header = TRUE)
colnames(cMat2)
colnames(cMattest)

test3 <- merge(cMat2, cMattest, all=TRUE)
dim(test3)

rNames <- cMat2[, 'X']
cNames <- colnames(cMat2)
cMat2 <- as.matrix(cMat2[, !cNames%in%"X"])
dim(cMat2)
rownames(cMat2) <- cNames[-1]

indx <- sample(cNames[-1], length(cNames[-1]))
indx
#cMat2[indx,indx]
det(-cMat2)
det(-cMat2[indx,indx])

eigVals <- eigen(cMat2[indx,indx], symmetric=FALSE, only.values = TRUE)$values
all(Re(eigVals)<=0)
eigVals

# try with old matrix
cMat1 <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/commMat.csv',
  header = TRUE)

rNames <- cMat1[, 'X']
cNames <- colnames(cMat1)
cMat1 <- as.matrix(cMat1[, !cNames%in%"X"])
dim(cMat1)
rownames(cMat1) <- cNames[-1]

indx <- sample(cNames[-1], length(cNames[-1]))
indx
#cMat1[indx,indx]
det(-cMat1)
det(-cMat1[indx,indx])

eigVals <- eigen(cMat1[indx,indx], symmetric=FALSE, only.values = TRUE)$values
all(Re(eigVals)<=0)
eigVals
#-------------------------------------------------------
# Test re-arrangement of row-column pairs from Dambacher et al. 2003 Fig. 2

modelA <- matrix(c(-1,1,0,0,0,
                   -1,-1,1,0,0,
                   0,-1,-1,1,0,
                   0,0,-1,-1,1,
                   0,0,0,-1,-1),
                 nrow=5, ncol=5, byrow=TRUE)
rownames(modelA) <- LETTERS[1:5]
colnames(modelA) <- LETTERS[1:5]
det(-modelA)
modelAprime <- modelA[c("B","C","A","D","E"),c("B","C","A","D","E")]
det(-modelAprime)
modelAprime <- modelA[c("D","E","B","C","A"),c("D","E","B","C","A")]
det(-modelAprime)
modelAprime <- modelA[c("D","E","B","C","A"),c("D","E","B","C","A")]
det(-modelAprime)

indx <- sample(LETTERS[1:5], 5)
indx
modelA[indx,indx]
det(-modelA[indx,indx])

modelL <- matrix(c(-1,1,0,0,0,
                   -1,-1,0,1,1,
                   0,0,-1,1,0,
                   0,0,-1,-1,1,
                   1,-1,1,-1,-1),
                 nrow=5, ncol=5, byrow=TRUE)
rownames(modelL) <- LETTERS[1:5]
colnames(modelL) <- LETTERS[1:5]
det(-modelL)

indx <- sample(LETTERS[1:5], 5)
indx
modelL[indx,indx]
det(-modelL[indx,indx])

#----------------------------------------------------------
# Test calculations for a modified simple model that includes profits, employment, and seafood nodes, but a single commercial fishery
library(LoopAnalyst)
simpleMat <- read.csv('/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/simpCommMat2.0_wPrfEmplSeaf.csv', 
                      header = TRUE)
rNames <- simpleMat[, 'X']
simpleMat <- as.matrix(simpleMat[, -1])
simpleMat <- apply(simpleMat, 2, as.numeric)
rownames(simpleMat) <- rNames

# determine feedback and stability
feedback(simpleMat)
det(-simpleMat)
#(-1^(nrow(simpleMat)+1))*det(simpleMat)
# assess Lyapunov stability (Dambacher et al. 2003 Am Nat)
eigVals <- eigen(simpleMat, symmetric=FALSE, only.values = TRUE)$values
all(Re(eigVals)<=0)

# negative inverse matrix
#solve(-simpleMat)

# Get adjoint matrix
simpAdjMat <- make.adjoint(simpleMat, status = TRUE)

write.csv(simpAdjMat, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/simpAdjMat2.0_wPrfEmplSeaf.csv')

# Get absolute feedback matrix
simpTotMat <- make.T(simpleMat)

write.csv(simpTotMat, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/simpTotMat2.0_wPrfEmplSeaf.csv')

# get the weighted feedback matrix
simpWeight <- make.wfm(simpleMat, status = TRUE)

write.csv(simpWeight, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/simpWeight2.0_wPrfEmplSeaf.csv')
