#31/07/2014

#Morphometric analyses of the full data sets
  #preparation for the convergence analyses
  
#Steps in analysis
  #1) Read in and clean up raw landmark data
  #2) Procrustes superimposition of all specimens
  #3) Find the average Procrustes shape coordinates for each species
  #4) PCA of the shape coordinates
  
#Output
  #1) PCA plot coloured by family (_allfam_PCA)
  #2) Average shape coordinates for each species (sps.mean)
  #3) Taxonomy for these shape coordinates (sp.fam)
  #4) Output of the PCA analysis (sps.meanPCA)
  
library(geomorph)
#General functions
  source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r")

########################################################
#READ IN DATA
  setwd("C:/Users/sfinlay/Desktop/Thesis/Convergence/data")

#SkDors data
    #1) Landmarks
        #landmarks + curves file with the control lines removed
        #land <- readland.tps(file="skdors/Skdors_16_12_13_10landmarks+4curves_edited.TPS")

    #2) Sliders
        #edited sliders file (top 2 rows removed and the words before slide after put in instead
        #curves <- as.matrix(read.table("skdors/Skdors_16_12_13_10landmarks+4curves_sliders_edited.NTS", header=TRUE))

    #3) Taxonomy
        #file that has the correct taxonomy for each of the images
        #taxa <- read.csv ("skdors/Skdors_16_12_13_10landmarks_images+specimens.csv" , header=T)

    #4) Specimens to remove
        #Null

#SkVent data
    #1) Landmarks
        #land <- readland.tps(file="skvent/SkVent_30_10_13_13landmarks+1curve_edited.TPS")
    #2) Sliders
        #curves <- as.matrix(read.table(file="skvent/SkVent_1skull_13landmarks+1curve_sliders_edited.tps", header=TRUE))     #this is a tps file and the others are nts but it doesn't make a difference
    #3) Taxonomy
        #taxa <- read.csv("skvent/SkVent_30_10_13_imagelist+specimens.csv" , header=TRUE)
    #4) Specimens to remove
        #rem <- read.csv("skvent/SkVent_remove_spec.csv", header=T)
        
#SkLat data
    #1) Landmarks
        #land <- readland.tps(file="sklat/SkLat_08_11_13_9landmarks_2curves_edited.TPS")
    #2) Sliders
        #curves <- as.matrix(read.table(file="sklat/SkLat_08_11_13_9landmarks_2curves_sliders_edited.NTS", header=TRUE))
    #3) Taxonomy
        #taxa <- read.csv("sklat/SkLat_08_11_13_Specimens+images.csv", header=TRUE)
    #4) Specimens to remove
        #rem <- read.csv("sklat/SkLat_remove_spec.csv", header=T)
        

#Mandibles data
    #1) Landmarks
        land <- readland.tps(file="mands/Mands_14_03_2014_7landmarks+4curves_edited.TPS")
    #2) Sliders
        curves <- as.matrix(read.table("mands/Mands_14_03_2014_7landmarks+4curves_sliders_edited.txt", header=TRUE))
    #3) Taxonomy
        taxa <- read.csv("mands/Mands_14_03_2014_Images+Specimens.csv", header=T)
    #4) Specimens to remove
        rem <- read.csv("mands/Mands_remove_spec.csv", header=T)
        
#######################################
#CLEAN UP THE DATA
#######################################

#Combine the taxonomic information with the array of landmark data
  combine <- list(land=land, curves=curves, ID=taxa$ID,SpecID=taxa$SpecID, Order=taxa$Order_05,
                  Fam=taxa$Family_05, Genus=taxa$Genus_05, Species=taxa$Species_05, Binom=taxa$Binomial_05)

# Remove the _sp. specimens
  sp <- which(combine$Species=="sp.")

  combine <- remove.from.list(combine, sp)
  combine <- droplevels.from.list(combine)

#********************************************
#Remove broken or damaged specimens
  #doesn't apply to the skdors data because rem is NULL

#find the ID numbers of specimens with missing data
  matching <- matching.id(rem$SpecID, combine$SpecID)
    combine <- remove.from.list(combine, matching)
    combine <- droplevels.from.list(combine)
#*********************************************
#Use all of the families
  mydata <-combine
  
  #Or remove the Notoryctidae because there's only one specimen
      #mydata <- remove.from.list(combine, which(combine$Fam=="Notoryctidae"))
      #mydata <- droplevels.from.list(mydata)
      
#######################################
#PROCRUSTES SUPERIMPOSTION
#######################################

#General Procrustes Alignment of all of the scaled coordinates
  mydataGPA <- gpagen(mydata$land, curves=mydata$curves, ProcD=TRUE,)


#List the coordinates with the taxonomic information
  Proc.co <- list(coords=mydataGPA$coords,csize=mydataGPA$Csize,ID=mydata$ID,SpecID=mydata$SpecID,
                  Order=mydata$Order, Fam=mydata$Fam, Genus=mydata$Genus, Species=mydata$Species, Binom=mydata$Binom)

#######################################
#SPECIES AVERAGING
#######################################

#group the arrays of coordinates according to species
  group.sps.coords <- species.coordinates(Proc.co$coords, Proc.co$Binom)

#average coordinate values for each species
  sps.mean <- mean.coords(group.sps.coords)

#######################################
#PRINCIPAL COMPONENTS ANALYSIS
#######################################

sps.meanPCA <- plotTangentSpace(sps.mean$meanshape, axis1 = 1, axis2 = 2,warpgrids = TRUE, label = TRUE)

#Re-create the PCA plot
  xaxis <- sps.meanPCA$x[,1]
  yaxis <- sps.meanPCA$x[,2]

#colour points by family
  #data frame of the unique family and binomical combinations
  sp.fam <- as.data.frame(unique(cbind(as.matrix(Proc.co$Fam), as.matrix(Proc.co$Binom))))
    colnames(sp.fam) <- c("Family","Binomial")


#palette of colours to go with the family
  #levels(sp.fam$Fam) gives the order of the families; alphabetical
  Fam.col <- c("magenta", "black", "forestgreen", "saddlebrown", "darkorchid", "dodgerblue","red1")
  #set this new vector as the colour palette for the PCA plots
  palette(Fam.col)

########################
#OUTPUT
#####################
setwd("C:/Users/sfinlay/Desktop/Thesis/Convergence/output")

#Save the PCA graph without the legend (just gets in the way)
  #pdf is better quality than jpeg or tiff


#SkDors
  #pdf(file="skdors/skdors_allfam_PCA.pdf")
#SkVent
  #pdf(file="skvent/skvent_allfam_PCA.pdf")
#SkLat
  #pdf(file="sklat/sklat_allfam_PCA.pdf")
#Mands
  #pdf(file="mands/mands_allfam_PCA.pdf")

  plot(xaxis, yaxis, xlab="", ylab="", las=1,
       col=sp.fam$Family, pch=16, bty="l", cex.lab=1, cex=1.2, xaxt="n", yaxt="n")
    #draw the min,max and 0 values on the x axis
      axis(side=1, at=c(round(min(xaxis),3), 0, round(max(xaxis),3)), las=1, cex=1.5)
    #same for the y axis
      axis(side=2, at=c(round(min(yaxis),3), 0, round(max(yaxis),3)), las=1, cex=1.5)
    #add dotted lines along 0,0
      abline(0,0, h=0, v=0, lty=2, lwd=1)

      legend('topright', legend = levels(sp.fam$Family), col = Fam.col,  border="black",
              cex = 0.8, pch = 16, bty="o", box.lty=1, ncol=1)
      #identify points on the graph
      #identify(xaxis,yaxis,labels=(sp.fam$Binom))

  #dev.off()

#Just the legend on the background of a blank graph (need to crop the picture later)
  #plot(xaxis, yaxis, xlab="", ylab="", las=1, col="white", bty="n", cex.lab=1, cex=1.2, xaxt="n", yaxt="n")
  #legend('center', legend = levels(sp.fam$Family), col = Fam.col,  border="black",
              #cex = 0.8, pch = 16, bty="o", box.lty=1, ncol=1)


#Data: average shape coordinates, accompanying taxonomy and sps.meanPCA results

#SkDors
  #dput(sps.mean, file="skdors/skdors_allfam_sps.mean.txt")
  #write.table(file="skdors/skdors_allfam_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #dput(sps.meanPCA, file="skdors/skdors_allfam_sps.meanPCA.txt")

#SkVent
  #dput(sps.mean, file="skvent/skvent_allfam_sps.mean.txt")
  #write.table(file="skvent/skvent_allfam_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #dput(sps.meanPCA, file="skvent/skvent_allfam_sps.meanPCA.txt")

#SkLat
  #dput(sps.mean, file="sklat/sklat_allfam_sps.mean.txt")
  #write.table(file="sklat/sklat_allfam_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #dput(sps.meanPCA, file="sklat/sklat_allfam_sps.meanPCA.txt")

#Mands
  #dput(sps.mean, file="mands/mands_allfam_sps.mean.txt")
  #write.table(file="mands/mands_allfam_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #dput(sps.meanPCA, file="mands/mands_allfam_sps.meanPCA.txt")

    
    