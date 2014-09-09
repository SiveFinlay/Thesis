#09/09/2014

#Error checking for the skvent data
#7 different skulls, 3 pictures of each, 3 replicates of each picture

library(geomorph)

source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r")
#source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/DisparityFunctions_Variance_Range.r")
#source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/PvalueFunction_FromDistribution.r")
########################################################
#READ IN DATA;
########################################################
  setwd("C:/Users/sfinlay/Desktop/Thesis/Methods/analyses")

#1) Landmarks
#landmarks + curves file with the control lines removed
  land <- readland.tps(file="Skvent_7skulls_3pics_3reps_13landmarks+curve_edited.TPS")

#2) Sliders
#edited sliders file (top 2 rows removed and the words before slide after put in instead
  curves <- as.matrix(read.table("Skvent_7skulls_3pics_3reps_13landmarks+curve_edited_sliders_edited.NTS", header=TRUE))

#3) Taxonomy
#file that has the correct taxonomy for each of the images
  taxa <- read.csv ("Skvent_7skulls_3pics_3reps_13landmarks+curve_images+taxonomic.csv" , header=T)

#Combine the data
    mydata <- list(land=land, curves=curves, ID=taxa$ID, SpecID=taxa$SpecID, Order=taxa$Order,
                  Fam=taxa$Family, Genus=taxa$Genus, Species=taxa$Species, Binom=taxa$Binomial)

#####################
#PROCRUSTES
#####################
  mydataGPA <- gpagen(mydata$land, curves=mydata$curves, ProcD=TRUE,)
  
#List the coordinates with the taxonomic information
  Proc.co <- list(coords=mydataGPA$coords, csize=mydataGPA$Csize, ID=mydata$ID, SpecID=mydata$SpecID,
                  Order=mydata$Order, Fam=mydata$Fam, Genus=mydata$Genus, Species=mydata$Species, Binom=mydata$Binom)

################
#PCA
#################
#NB: I'm skipping the species averaging step from the main analysis because I want to look at the overall variation
  mydata.PCA <- plotTangentSpace(Proc.co$coords, axis1 = 1, axis2 = 2,warpgrids = TRUE, label = TRUE)

  #Gives 7 very distinct clusters of points
  
###################
#Procrustes distance
################
#Compare the Procrustes distances within specimens to the distances between different skulls

#Select the PC axes that account for 95% of the variation
  PC95axes <- selectPCaxes(mydata.PCA, 0.956, Proc.co$Binom)
  
#Euclidean distance matrix of all of the species
  Euc.dist <- as.matrix(dist(PC95axes, method="euclidean", diag=FALSE, upper=FALSE))
  
#Find the min, max and mean Procrustes distances within each set of pictures
  #NB: The function selects the minimum non 0 number
  #(there will always be 0s in the distance matrix because there's no difference between an image and itself)
 var.within <- Proc.dist.within(Euc.dist)

#Mean values for digitisation error
  min.within <- mean(var.within[,1])
  max.within <- mean(var.within[,2])
  mean.within <- mean(var.within[,3])

#Min, max and mean Procrustes distances between all of the pictures
  min.between <- min(Euc.dist[which(Euc.dist > 0)])         #Minimum non 0 difference
  max.between <- max(Euc.dist)
  mean.between <- mean(Euc.dist)
  
#Ratios of the variation within speicimens to variation between specimens
  var.ratio <- matrix(nrow=1, ncol=3)
    colnames(var.ratio) <- c("min", "max", "mean")
      var.ratio[1,1] <- (min.within/min.between)*100
      var.ratio[1,2] <- (max.within/max.between)*100
      var.ratio[1,3] <- (mean.within/mean.between)*100

# I think this means that 3.7% of the observed differences between species actually comes from morphometric error?
#Need to do the nested ANOVA test next