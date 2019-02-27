#!/bin/bash/R

########### INPUT FOR WORKING ############
otutable <- "/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/MF_and_OTU_edited/otu_table_text.txt"
MFFP <- "/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/MF_and_OTU_edited/MF_withalpha.txt"
dmFP <- "/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/beta_div/bray_curtis_dm.txt"


############## LOAD AND FILTER DATA #####################
library(MASS)
library(vegan)
set.seed(10234)
setwd("/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/")


# Load mappingfile and dm and otutable
MF <- read.delim(paste0(MFFP), header = TRUE, row.names = 1)
dm <- read.delim(paste0(dmFP), header=TRUE, row.names=1)

# Look at dm
any(dm < 0) # No negative values
all(rowSums(dm) == colSums(dm)) # symmetric; row sums = col sums

# Make smaller dm to test out with
smdm <- dm[1:5,1:5]

# Make PCOA 
??pcoa
