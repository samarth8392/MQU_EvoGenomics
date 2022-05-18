###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 01/11/22                  Last Modified: 02/07/22 ###
###########################################################################
###########################################################################
###                   population_tree.R                                 ###
###########################################################################

#### PREREQUISITES #####
library(ggplot2)
library(plyr)
library(dplyr)
library(gridExtra)
library(ggtree)
library(phytools)
library(RColorBrewer)
library(phangorn)
library(corrplot)
library(rcartocolor)
# Load colors
#display_carto_all(colorblind_friendly = TRUE)

col_az <- carto_pal(12,"Vivid")[11]
col_wtx <- carto_pal(7,"Purp")[5]
col_etx <-carto_pal(7,"Peach")[5]
col_mx <- "#3598FF"
colss <-c(rep(col_az,28),rep(col_mx,4),rep(col_wtx,31),rep(col_etx,3))
col_azwtx <- "#90CCC2"
cols <- c(col_az,col_wtx,col_etx, col_mx,col_azwtx)
cols_l <- adjustcolor(cols, alpha.f = 0.5)

setwd("~/Documents/Thesis_Research/Final Results/Ch3/revise/tassel/")

poplist <- read.table("best66.samples.list",header = T)

temp_data <- read.tree("best66.tree.txt")
temp_data <- ape::root(temp_data, node = 74)
sample_idx_g <- match(temp_data$tip.label, poplist$SampleID)
loc.name_g <- poplist$Pop[sample_idx_g]
temp_data$tip.label <- as.character(loc.name_g)
groupInfo <- split(temp_data$tip.label, gsub("_\\w+", "", as.character(loc.name_g)))
tree_g2 <- groupOTU(temp_data, groupInfo, group_name = "Location")
tree_g2 <- ape::root(tree_g2, node = 68)
par(mar = c(0.1, 0.1,0.1,0.1))

p <- ggtree(tree_g2, layout="daylight",ladderize=T, aes(color=Location),right = T,branch.length="none", size=1.5) + 
  geom_tippoint(size=4)+#geom_tiplab(size=3, hjust = -0.3)+
  theme(legend.position = "none")+#geom_text(aes(label=node), hjust=-.3)+
  scale_color_manual(values=c("black",cols[1],cols[3],cols[4],cols[2]))+
  theme(plot.margin = unit(c(4,4,4,4), "mm"))
p#+geom_treescale(-2,2)

# From co-ancestry tree
library(GPArotation)
library(paran)
source("~/Documents/GitHub/Mathur_ch3/RScripts/FinestructureRcode/FinestructureLibrary.R")
source("~/Documents/GitHub/Mathur_ch3/RScripts/FinestructureRcode/FinestructureDendrogram.R")
setwd("~/Documents/Thesis_Research/Final Results/Ch3/fineS/run2/")
treefile<-"xml/best66.noSing.chr1_linked_tree.xml" ## finestructure tree file
treexml<-xmlTreeParse(treefile) ## read the tree as xml format
ttree<-extractTree(treexml) ## extract the tree into ape's phylo format

ttree <- ape::root(ttree, node = 74)

sample_idx_g <- match(ttree$tip.label, poplist$SampleID)
loc.name_g <- poplist$Pop[sample_idx_g]
ttree$tip.label <- as.character(loc.name_g)
groupInfo <- split(ttree$tip.label, gsub("_\\w+", "", as.character(loc.name_g)))
tree_g2 <- groupOTU(ttree, groupInfo, group_name = "Location")
tree_g2 <- ape::root(tree_g2, node = 67)
par(mar = c(0.1, 0.1,0.1,0.1))

p <- ggtree(tree_g2, layout="daylight",ladderize=T, aes(color=Location),right = T,branch.length="none", size=1.5) + 
  geom_tippoint(size=4)+#geom_tiplab(size=3, hjust = -0.3)+
  theme(legend.position = "none")+#geom_text(aes(label=node), hjust=-.3)+
  scale_color_manual(values=c("black",cols[1],cols[3],cols[4],cols[2]))+
  theme(plot.margin = unit(c(4,4,4,4), "mm"))
p+geom_treescale(-2,2)

