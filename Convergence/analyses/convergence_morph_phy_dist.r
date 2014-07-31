#31/07/2014

#Compare morphological and phylogenetic distances

#Pairwise comparisons among tenrecs and all other species; not interested in comparisons that don't involve tenrecs


#Use output data from the convergence_procrustes and convergence_phylogenies scripts
  setwd("C:/Users/sfinlay/Desktop/Thesis/Convergence/output")

library(ape)

source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r")


#Input data: taxonomy, PCA results, phylogenies

#Skdors
   sp.fam <- read.table(file="skdors/skdors_allfam_sps.mean_taxonomy.txt")
   sps.meanPCA <- dget("skdors/skdors_allfam_sps.meanPCA.txt")
   mytrees <- read.trees("phylogenies/skdors_allspecies_100trees.phy")


########################################
#MORPHOLOGICAL DISTANCE MATRIX
#######################################
#Morphological distance matrix calculated from PC scores (cf calculating disparity from PC scores)

#Select PC axes that account for 95% of the variation
    PC95axes <- selectPCaxes(sps.meanPCA, 0.956, sp.fam$Binomial)

#Euclidean distance matrix from these axes
  morph.dist <-dist (PC95axes, method="euclidean", diag=FALSE, upper=FALSE)
#turn it into an ordinary matrix
  morph.dist <- as.matrix(morph.dist)
  
#Order the distance matrix alphabetically by rowname
  morph.dist.order <- morph.dist[order(rownames(morph.dist)),]
  
#Order the taxonomic matrix alphabetically by binomial name
  sp.fam.order <- sp.fam[order(sp.fam$Binomial),]
  
#Select comparisons that involve tenrecs: tenrec vs. tenrec or tenrec vs. other species
  #i.e. remove comparisons that don't involve tenrecs

  tenrec.morph.comp <- morph.dist.order[which(sp.fam.order$Family=="Tenrecidae"),]
#-------------------------------------------------------------
#Phylogenetic distance matrices
  phylo.dist <- NULL
    for (i in 1:length(mytrees)){
      phylo.dist[[i]] <- cophenetic(mytrees[[i]])
    }

#Order all of the phylogenetic distance matrices alphabetically rowname (species)
  phylo.dist.order <- NULL
    for (i in 1: length(phylo.dist)){
      phylo.dist.order[[i]] <- phylo.dist[[i]][order(rownames(phylo.dist[[i]])),]
    }

#Select just the comparisons that involve tenrecs
  tenrec.phylo.comp <- NULL
    for (i in 1: length(phylo.dist.order)){
      tenrec.phylo.comp[[i]] <- phylo.dist.order[[i]][which(sp.fam.order$Family=="Tenrecidae"),]
    }
#----------------------------------------------

#Compare morphological and phylogenetic distances

#Doesn't make sense to just do a direct comparison like this
plot(tenrec.phylo.comp[[52]], tenrec.morph.comp)

#I need to do some more meaningful analyses
  #e.g. Wheatsheaf index or Stayton's weighted count metric
  #(Muschick's method requires re-visiting the issue of shape simulations)
  
  #Wheatsheaf index will probably be the best one: strength of convergence in focal groups
    #Groups:
      #Spiny:Echinops, Setifer, Hemicentetes?, Tenrec vs. hedgehogs, gymnures, solenodons
      #Mole-like: Oryzorcites vs. talpidae, golden moles
      #Shrew-like: Microgale, Geogale vs. soricidae
      #Aquatic?: Limnogale, potamogale vs. gymnures, notoryctes (probably don't have a good enough sample size)