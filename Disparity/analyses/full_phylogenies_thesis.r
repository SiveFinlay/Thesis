#04/03/15
#Phylogenies for the thesis
#Phylogenies of all golden moles and tenrecs
#Use them to show how my species sampling is distributed across the phylogeny

#Steps
#1) Read in data
#2) Prune trees partially (remove genera that aren't in the data)
#3) Add missing species to the trees

library(ape)
library(geiger)
library(phytools)

source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r")
source("http://www.phytools.org/add.species.to.genus/v0.1/add.species.to.genus.R")

#distribution of mammal supertrees with resolved polytomies from Kuhn et al 2011
  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity")
  supertree<-read.nexus("data/FritzTree.rs200k.100trees.tre")
#--------------------------------------------------------------
#Read in data
  #This example uses the skulls dorsal data but all analyses used included the same species

#1) Mean shape coordinates for each species
  sps.mean <- dget(file="output/shape_data/skdors/SkDors_tenrec+gmole_sps.mean.txt")
#2) Taxonomy for those coordinates
  sps.tax <- read.table(file="output/shape_data/skdors/SkDors_tenrec+gmole_sps.mean_taxonomy.txt")

##############################################
#Extract all species in the trees
  all.species <- supertree[1]$tip.label[[1]]

#Matrix of genus and species names of all of the tip labels
  sup.spe <- split.binom(all.species)

#Matrix of genus and species names for species in mydata
  morph.spe <- split.binom(as.character(sps.mean$Binom))

###############################################
#PRUNE TREES PARTIALLY
###############################################

#Genera in the morphological data that are currently in the phylogenies
  common <- intersect(sup.spe$Genus, morph.spe$Genus)
  
#Add golden mole genera that are not in the morphological data
  #Wilson and Reeder, mammal species of the world 2005 includes Neamblysomus
  #All other MSW genera are represented within the morphological data
  
  common <- c(common, "Neamblysomus")
  
  
#Select just those genera from the supertrees
  com.gen <- as.list(rep(NA, length(common)))
    for (i in 1: length(common)){
      com.gen[[i]] <- which(sup.spe$Genus == common[i])
    }

  com.id <- unlist(com.gen)

#list of id numbers that are not common
  lose.tip <- all.species[-com.id]

#Prune the trees down to a more manageable size

   sup.prun <- NULL
   for (i in 1:length(supertree)){
    sup.prun[[i]] <- remove.missing.species.phy(supertree[[i]], lose.tip)
    }
  class(sup.prun) <- "multiPhylo"


###############################################
#ADD MISSING SPECIES TO THE TREES
###############################################

#Microgale in the super trees and in the morphological data
  sup.mic <- sup.spe[which(sup.spe$Genus == "Microgale"),]
  morph.mic <- morph.spe[which(morph.spe$Genus == "Microgale"),]

#Microgale missing from the trees
  miss.mic <- morph.mic[-(which(morph.mic[,2] %in% sup.mic[,2])),]
  microgale <- paste(miss.mic$Genus, miss.mic$Species, sep="_")

#  "Microgale_jobihely"     "Microgale_soricoides"   "Microgale_fotsifotsy"
 #"Microgale_gymnorhyncha" "Microgale_monticola"    "Microgale_grandidieri"


#add the missing microgale species at random to the rest of the microgale genus
  sup.prun.micro <- sup.prun
    for(i in 1:length(sup.prun.micro)){
      for (m in 1:length(microgale)){
        sup.prun.micro[[i]] <- add.species.to.genus(sup.prun.micro[[i]],microgale[m],where="random")
      }
    }

#Find the missing golden mole species

#species in the pruned trees
  prun.spe <- split.binom(sup.prun.micro[[1]]$tip.label)

#golden moles in the morphological data
  morph.gm <- subset.matrix(sps.tax, sps.tax$Family, "Chrysochloridae")
  morph.gm.spe <- split.binom(as.character(morph.gm[,2]))


#golden moles missing from the trees
  miss.gm <- morph.gm.spe[-(which(morph.gm.spe[,2] %in% prun.spe[,2])),]
  miss.gm <- paste(miss.gm$Genus, miss.gm$Species, sep="_")

  #Eremitalpa granti and Cryptochloris wintoni

#golden mole Genera currently in the trees (the first 17 rows in sup.prun.micro[[1]]$tip.label)
  tree.gm <- setdiff(sup.prun.micro[[1]]$tip.label[1:17], miss.gm)

#Add the first missing golden mole species
  sup.prun.micro.gm <- add.species.to.MRCA(sup.prun.micro, tree.gm, miss.gm[1])

#Add the second missing golden mole species
  sup.prun.micro.gm2 <- add.species.to.MRCA(sup.prun.micro.gm, tree.gm, miss.gm[2])
  class(sup.prun.micro.gm2) <- "multiPhylo"

#Resolve polytomies at random
  mytrees <- calc.each.list(sup.prun.micro.gm2, multi2di)
  
##########################################
#Identify species that are only in the trees  and not in the morphological data

  TreeOnly <- tree.only(sup.prun.micro.gm2.res, sps.mean$Binom)
  
  #Golden moles in the tree only
  gmole.to <- TreeOnly[[1]]
  
#"Calcochloris_tytonis"       "Chlorotalpa_sclateri"       "Amblysomus_septentrionalis"
#"Amblysomus_marleyi"         "Neamblysomus_gunningi"      "Neamblysomus_julianae"
#"Chrysochloris_visagiei"

#Mark these tip labels differently to the rest of the labels in bold

#Plot one of the distribution of trees

#Select the tree only golden moles
ident.tip <- NULL

for (i in 1:length(gmole.to)){
  ident.tip[[i]] <- which(mytrees[[10]]$tip.label == gmole.to[i])
  }

#Create a font identity vector
font.label <- rep(1, length(mytrees[[10]]$tip.label))

#ident.tip species get a different font (bold italics)
font.label[ident.tip] <- 4

  

plot(mytrees[[10]], cex=0.75, font=font.label)

#Saved the tree as tc+gm_phylo


