#!/bin/bash/R

############ Analysis for 5 species dataset ############

#### FILTERING BD NOTES:
# First, checks all numbers to see if less than indivTHRESH. If less than indivTHRESH, it is changed to '0'
# Then, it see is if at least 2 are NOT zero and the third is more than 50. 
# If the third is less than 50 AND the other two measurements are zero, they are all changed to zeros.

library(MASS) # for isoMDS
library(vegan) # for betadisper 
library(tidyverse) # for data manipulation
library(ggplot2)
library(gridExtra) # for grid arrange
library("rstanarm")
library(rstan) #extract


#################################
##### PATHWAYS #######
homedir = "/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset"
output = "Preliminary"
# otuPWD = "MF_and_OTU_edited/otu_table_text.txt"
otuPWD = "MF_and_OTU_edited/otu_table_r5000.txt"
mfPWD = "MF_and_OTU_edited/MF_withalpha_osseadded.txt"
inhibPWD <- "ANTIFUNGAL/inhibitory_metadata_MANUAL.txt"

indivTHRESH = 5 # BD individual threshold
thirdTHRESH = 50 # BD 3rd sample threshold
thresh = 1 # number of OTUs in the table for OTU to be kept
minOTUSample = 5
minOTUTab = 100


dmPWD = "beta_div/bray_curtis_dm.txt"
############## LOAD AND FILTER DATA #####################

setwd(homedir)
dir.create(output)
setwd(output)

# Load mappingfile and dm and otutable
mf <- read.delim(paste0("../",mfPWD), header = TRUE, as.is = TRUE)
# Below line is commented out bc it takes a long time to load and I am trying to save time
otu <- read.delim(paste0("../",otuPWD), header=TRUE, as.is = TRUE, skip=1)
# dm 
dm <- read.delim(paste0("../",dmPWD), header=1, row.names=1, as.is=TRUE)
# inhib
inhib <- read.delim(paste0("../",inhibPWD), header=FALSE, as.is=TRUE)

##### ADJUSTING DATA: Make same order, make 'time' variable ######
{
    
    #### Edit metadata/mapping file #####
    toKeepMF <-  c("X.SampleID"
                   ,"Bd_Run_1","Bd_Run_2","Bd_Average_Run_3" # infection levels; will transform to be average
                   ,"SPEC" # Species code
                   ,"TREATMENT_GROUP" # whether it was pre or post; should rename
                   ,"BD100KZSP_3122011" # whether or not individual was exposed to BD: will rename
                   ,"ANONYMIZED_NAME" # individual frog ID
                   , "PRE_POST_BD_NUM" # timepoint
                   ,"shannon_even_5000_alpha" #shannon div
                   ,"observed_otus_even_5000_alpha" #observed otus (richness)
    )
    newNames <- c("SampleID"
                  , "BD1","BD2","BD3"
                  , "species"
                  , "prepost"
                  , "BD_infected"
                  , "toadID"
                  , "timepoint"
                  , "shannon"
                  , "observed_otus"
                  )
    mf.tb <- mf %>%
        as_tibble() %>% # make into tibble
        dplyr::select(one_of(toKeepMF)) %>% # filter to only relevant variables
        rename_at(vars(toKeepMF), ~ newNames) %>%#rename variable names
        filter(prepost == "Pre" | prepost == "Pos") %>% # get rid of things in "tre", which is a different experiment
        separate(timepoint, into=c("exposure","time"), sep="_", convert=TRUE) %>% # create a time variable
        mutate( time = ifelse(exposure == "Pre", time, time + 5)) %>% # time variable "restarts" at BD exposure point, but I want it to be continuous
        filter(species != "None")
    
    ### Check to make sure there are no duplicates
    #+ message=FALSE, echo=FALSE, warnings=FALSE, include=FALSE
    any(table(dplyr::select(mf.tb, toadID, time)) > 1)
    # Oh no! There is a duplicate! Check which one it is.
    #+ message=FALSE, echo=FALSE, warnings=FALSE, include=FALSE
    table(dplyr::select(mf.tb, toadID, time))>1
    # Raca10 has 2 timepoint 3's
    # Buma11 has 2 time point 3's
    # I think this is just a mistake in the data entry. I've gone through the rest of the metadata and determined that 
    #   it is probably mis-labeled so:
    # I'm going to manually change:
    # I think mck115Pre.3Buma.11.452405 Buma11 timepoint 3 is actually Buma9, time 3
    # I think mck105Pre.3Raca.10.451190 Raca10 timepoint 3 is actually Raca9, time 3
    mf.tb[mf.tb$SampleID == "mck115Pre.3Buma.11.452405",c("toadID")] <- "Buma_9"
    mf.tb[mf.tb$SampleID == "mck105Pre.3Raca.10.451190", c("toadID")] <- "Raca_9"
    
    # Get a "mean" Bd infection number
    indivTHRESH <- 5
    thirdTHRESH <- 50
    # If 2 BD samples are 0 and the third is less than 50, then set to zero.
    # Also, if anything is less than 5, make it zero anyway.
    allBd <- mf.tb %>% 
        dplyr::select(BD1,BD2,BD3)
    for ( r in 1:nrow(allBd) ) {
        for ( c in 1:ncol(allBd) ) {
            if ( (allBd[r,c] < indivTHRESH) | is.na(allBd[r,c]) ) {
                allBd[r,c] <- 0
            }
        }
        if ( (sum(allBd[r,] == 0) == 2) & (max(allBd[r,],na.rm=TRUE) < thirdTHRESH) ) {
            allBd[r,] <- c(0,0,0)
        }
        
    }
    
    #rename original so it doesn't get lost
    mf.tb <- mf.tb %>%
        rename(BD1_raw=BD1, BD2_raw=BD2, BD3_raw=BD3) %>%
        mutate(BD3_raw = if_else(is.na(BD3_raw),0,BD3_raw))
    mf.tb[,c("BD1","BD2","BD3")] <- allBd
    # Now, make an "average" bd load
    mf.tb <- mf.tb %>%
        rowwise() %>%
        mutate(aveBD = mean(c(BD1,BD2,BD3)), aveBD = max(c(BD1,BD2,BD3))) %>%
        # Finally, make BD presence absence
        mutate(PABD = ifelse(aveBD>0, 1, 0))
    
    # copy mf so that there is an unfiltered version
    mf.tb.unfilt <- mf.tb
    # filter to ONLY include those things in otu table
    mf.tb <- mf.tb %>%
        filter(SampleID %in% colnames(otu))
    
    #+ include=FALSE
    mf.tb$shannon <- as.numeric(mf.tb$shannon)
    mf.tb$observed_otus <- as.numeric(mf.tb$observed_otus)
    mf.tb$logRich <- log(mf.tb$observed_otus)
    
    #Samples to keep in otu table
    OTU_names <- otu$X.OTU.ID
    otu.tb <- otu %>%
        as_tibble() %>%
        dplyr::select(one_of(c(mf.tb$SampleID))) %>%
        # replace(.<minOTUSample, 0) %>%
        # mutate(rowsums=rowSums(.)) %>%
        mutate(OTUID = OTU_names) #%>%
        # filter(rowsums > minOTUTab) %>%
        # dplyr::select(-rowsums)

    ##### Beta diversity turnover #####
    # Get distance to sample directly before for every sample
    # AKA: How similar was the current sample to the sample JUST before that time point for that individual toad?
    # Note; tried to do this with dplyr but it is VERY slow. Just used base R subsetting instead
    mf.tb$distance <- NA
    for ( i in 1:nrow(mf.tb) ) {
        currentSample <- mf.tb$SampleID[i]
        current <- mf.tb[i, c("toadID","time")]
        prevSample <- mf.tb$SampleID[mf.tb$toadID==as.character(current[1,1]) & mf.tb$time == as.numeric(current[1,2]-1)]
        if ( length(prevSample) > 0) {
            distTemp <- dm[currentSample,as.character(prevSample)]
            if ( length(distTemp) > 0) {
                mf.tb[i,"distance"] <- distTemp
                
            }
        }
        
        # print(paste0("done",i," out of ",nrow(mf.tb)))
    }
    
    
    ### Get proportion and count of inhibitory otus from inhibitory otu metadata
    inhib.tb <- inhib %>%
        as_tibble() %>%
        rename(Name=V1, inhib=V2, num=V3, seq=V4,sampleID=V5, taxa=V6)
    
    otu.tb.inhib <- otu.tb %>%
        mutate(inhib = (inhib.tb[match(otu.tb$OTUID, inhib.tb$seq),"inhib"])$inhib)
    
    # transpose it, and then make column for "percent inhib" and "count inhib"
    # total counts
    sampleCounts <- otu.tb.inhib %>%
        dplyr::select(-c("inhib","OTUID")) %>%
        colSums()
    # counts of total inhibitory sequences
    inhibCounts <- otu.tb.inhib %>%
        filter(inhib=="inhibitory") %>%
        dplyr::select(-c("inhib","OTUID")) %>%
        colSums()
    # richness of inhibitory sequences
    inhibRich <- otu.tb.inhib %>%
        filter(inhib=="inhibitory") %>%
        dplyr::select(-c("inhib","OTUID")) %>%
        replace(.>0, 1) %>%
        colSums()
    ## Now, add this information into the mf
    mf.tb.inhib <- mf.tb %>%
        mutate(n=sampleCounts[match(SampleID, names(sampleCounts))]
               ,inhibCounts =inhibCounts[match(SampleID, names(inhibCounts))]
               , inhibRich=inhibRich[match(SampleID, names(inhibRich))]) %>%
        mutate(percInhib = inhibCounts/n)
  
}
# mf.tb = mf[!(mf$TREATMENT_GROUP == "Tre"),]
# mf.tb <- mf %>%
#     select(TREATMENT_GROUP!="Tre")
# 
# ### Fix typos in dataset; make a list of all sample types (ie species names)
# # Get the species and control name and get rid of all numbers and underscores
# mf.tb[,"SampleType"] <- gsub('[0-9]|_','',mf.tb$SPECIES_NAME)
# # It also turns out there is one called 'Conrol' instead of "Control" so fix that too
# mf.tb[,"SampleType"] <- gsub("Conrol","Control",mf.tb[,"SampleType"])
# # Also, "cricket" vs "crickets"
# mf.tb[,"SampleType"] <- gsub("Crickets","Cricket",mf.tb[,"SampleType"])
# # Also, typo in sterile water
# mf.tb[,"SampleType"] <- gsub("Sterlile","Sterile",mf.tb[,"SampleType"])
# mf.tb <- mf.tb[rownames(mf.tb) %in% colnames(otu),]
# 
# ## Now, filter OTU Table to include only these samples
# otu.filt = otu[,rownames(mf.tb)]
# ## Get rid of low abundance OTUs; less than thresh
# otu.filt = otu.filt[rowSums(otu.filt) >= thresh,]
# 
# # Make a 'time' variable
# # Take the "PRE_POST_BD_NUM" variable in MF and change it:
# # If the "TREATMENT_GROUP" variable is "Pos", then you add 5 to the time variable
# # If the "TREATMENT_GROUP" variable is "Tre", then you add 16 to the time variable
# # The first 5 timepoints are "Pre" so you don't need to change this.
# # Now, you have a variable that is continuous across the entire treatment
# time <- as.numeric(gsub(".*_","",mf.tb$PRE_POST_BD_NUM))
# mf.tb <- cbind(mf.tb, time)
# # make a 'time2' variable for tracking across all time
# prepos <- mf.tb$TREATMENT_GROUP
# timetotal <- time
# for ( i in 1:length(prepos) ) {
#     if ( prepos[i] == "Pos" ) {
#         timetotal[i] <- as.numeric(timetotal[i]) + 5
#     } else if ( prepos[i] == "Tre" ) {
#         timetotal[i] <- as.numeric(timetotal[i]) + 16
#     }
# }
# mf.tb <- cbind(mf.tb, timetotal)
# 
# # Change "NA"s for BD run to '0'
# # BD runs 1 and 2 have '0's and no NAs, whereas BD run 3 has all NAs and seemingly no zeros.
# # I change them all to zero because I don't actually know which ones are 'true' NAs and which ones are not
# mf.tb$Bd_Average_Run_3[is.na(mf.tb$Bd_Average_Run_3)] <- 0
# mf.tb$Bd_Run_1[is.na(mf.tb$Bd_Run_1)] <- 0
# mf.tb$Bd_Run_2[is.na(mf.tb$Bd_Run_2)] <- 0
# }
# ########### Plot BD load in histogram ############
# {# Combine all BD values into a single table
# Bd1 <- mf.tb$Bd_Run_1
# Bd2 <- mf.tb$Bd_Run_2
# Bd3 <- mf.tb$Bd_Average_Run_3
# allBd <- cbind(Bd1,Bd2,Bd3) # Combine for later
# rownames(allBd) <- rownames(mf.tb)
# mf.tb[,c("Bd1","Bd2","Bd3")] <- allBd
# 
# # Now, chance all zeros to "NA" for this histogram because I can't differentiate
# # between a true zero and a "NA" zero
# Bd1[Bd1==0] <- NA
# Bd2[Bd2==0] <- NA
# Bd3[Bd3==0] <- NA


# Make a histogram of all these BD loads to see if there are 2 modes
# If there are 2 modes (or a strangely skewed single mode), that might suggest
# that there is some level of 'noise' that exists in data
# Then, we can use a cut-off to determine whether a BD read should actually be zero or not
    
pdf(file="histogram_bdload.pdf", height = 10, width =5)
    mf.tb %>%
        dplyr::select(BD1_raw,BD2_raw,BD3_raw) %>%
        # mutate(aveBD_raw = (BD1_raw+BD2_raw+BD3_raw)/3) %>%
        # mutate(logAveBD_raw = log(aveBD_raw)) %>%
        gather(key=Run, value=BD_load) %>%
        ggplot(mapping=aes(x=log(BD_load))) +
        geom_histogram(aes(y=..density..), bins=25) +
        facet_wrap(~Run,nrow=3)
dev.off()
    
pdf(file="histogram_bdload_afterAdjustments.pdf", height = 10, width =5)
mf.tb %>%
    dplyr::select(BD1,BD2,BD3) %>%
    # mutate(aveBD_raw = (BD1_raw+BD2_raw+BD3_raw)/3) %>%
    # mutate(logAveBD_raw = log(aveBD_raw)) %>%
    gather(key=Run, value=BD_load) %>%
    ggplot(mapping=aes(x=log(BD_load))) +
    geom_histogram(aes(y=..density..), bins=25) +
    facet_wrap(~Run,nrow=3)
dev.off()


######## Check to see how many toads/frogs were NEVER infected ########
infectionRate <- mf.tb %>%
    group_by(toadID, prepost, BD_infected) %>%
    summarise(maxBd1=max(BD1), maxBd2=max(BD2), maxBd3=max(BD3)) %>%
    ungroup() %>%
    mutate(prepost=as.character(prepost)) %>%
    arrange(as.character(prepost))
# infectionRate <- aggregate(allBd, by = list(CON = mf.tb$BD100KZSP_3122011, ID = mf.tb$ANONYMIZED_NAME, TREAT = mf.tb$TREATMENT_GROUP), FUN = max)
write.table(as.data.frame(infectionRate), file = "infectionRate.txt", sep = "\t")


####### Plot infection per amphibian over time #############
# I want a table and graph that shows each individual amphibian's infection load over time.

# Get list of species
speciesList <- levels(factor(mf.tb$species))

# Make a time list
timeList <- levels(factor(mf.tb$time))

# makes list of species over time; track infection occurance. Print out.
{
    ## NOTE: use mf.tb.unfilt because we want the 'zeros' from samples that werent diverse enough
    bdInfectTime <- list()
    mf.tb.unfilt[,c("maxInfect","meanInfect")] <- matrix(ncol = 2, nrow = nrow(mf.tb.unfilt))
    for ( sp in speciesList ) {
        # Make matrix
        n_s <- mf.tb.unfilt %>%
            filter(species==sp) %>%
            group_by(toadID) %>%
            summarise() %>%
            nrow()
        n_t <- mf.tb.unfilt %>%
            filter(species==sp) %>%
            pull(time) %>%
            max()
        tempMat <- matrix(ncol = n_t, nrow = n_s)
        colnames(tempMat) <- seq(1,n_t,by = 1)
        rownames(tempMat) <- seq(1,n_s, by = 1)
        tempMatAve <- matrix(ncol = n_t, nrow = n_s)
        colnames(tempMatAve) <- seq(1,n_t,by = 1)
        rownames(tempMatAve) <- seq(1,n_s, by = 1)
        for ( s in 1:n_s ) {
            for ( t in 1:n_t ) {
                bdload <- which((mf.tb.unfilt$time == t) & (mf.tb.unfilt$toadID == paste0(sp,"_",s)))
                if ( length(bdload) > 0 ) {
                    tempMat[paste0(s),paste0(t)] <- max((mf.tb.unfilt[bdload,c("BD1","BD2","BD3")]), na.rm=TRUE)
                    tempMatAve[paste0(s),paste0(t)] <- rowMeans(mf.tb.unfilt[bdload,c("BD1","BD2","BD3")], na.rm=TRUE)
                }
            }
        }
        bdInfectTime[[sp]] <- list()
        bdInfectTime[[sp]][["Max"]] <- tempMat
        bdInfectTime[[sp]][["Mean"]] <- tempMatAve
    }
    
    capture.output(bdInfectTime, file = "bdInfecTime.txt")
}


# Plot each individual taod to see how their infection changes over time
# Each individual is going to be a color of the rainbow
rainbowCol <- c("darkred","red","orange","yellow","green","darkgreen","blue","darkblue","purple","grey")
# Get the max log load of all samples so that we can compare across species
maxLoadLog <- log(max(unlist(bdInfectTime), na.rm=TRUE) +1)
# let list of controls vs not
which_con <- mf.tb[,c("toadID","BD_infected","time")] %>%
    group_by(toadID) %>%
    summarise(BD_infected=unique(BD_infected)) %>%
    separate(toadID, into=c("species","indiv"), remove = FALSE)
# Make a 6-panel plot for Max and Mean infection over time
{
pdf("InfectionOverTime_all_Max.pdf", height=5, width=7)
par(mfrow=c(2,3))
for ( sp in speciesList ) {
    tableTemp <- bdInfectTime[[sp]][["Max"]]
    minT <- 1
    maxT <- ncol(tableTemp)
    controls <- which_con %>%
        filter(species==sp) %>%
        mutate(indiv=as.numeric(indiv)) %>%
        arrange(indiv) %>%
        dplyr::select(BD_infected,indiv) %>%
        transmute(lty = if_else(BD_infected=="y",2,1), indiv=indiv) 
    plot(NULL, ylim=c(0,maxLoadLog), xlim=c(minT,maxT),xlab="Time", ylab="BD Load (ln)", xaxt="n", main = paste0(sp))
    axis(side=1, at = seq(minT,maxT), labels = seq(minT,maxT), cex.axis=0.8)
    for ( r in 1:nrow(tableTemp) ) {
        points(log(tableTemp[r,]+1), type = "l", col=rainbowCol[r], lty=pull(controls[controls$indiv==r,"lty"]))
    }
    abline(v = 6, lty = 2, col = "black")
}
# Last panel is the legend
plot(0,0,pch="",axes=FALSE,xlab="",ylab="")
legend("left"
       , legend= c(sapply(seq(1,9),function(x) paste0("Individual ",x)))
       , lty = 1
       , col = c(rainbowCol, NA)
       , bty = "n"
       , cex=0.75)
dev.off()

pdf("InfectionOverTime_all_Mean.pdf", height = 5, width=7)
par(mfrow=c(2,3))
for ( sp in speciesList ) {
    tableTemp <- bdInfectTime[[sp]][["Mean"]]
    minT <- 1
    maxT <- ncol(tableTemp)
    controls <- which_con %>%
        filter(species==sp) %>%
        mutate(indiv=as.numeric(indiv)) %>%
        arrange(indiv) %>%
        dplyr::select(BD_infected,indiv) %>%
        transmute(lty = if_else(BD_infected=="y",2,1), indiv=indiv) 
    plot(NULL, ylim=c(0,maxLoadLog), xlim=c(minT,maxT),xlab="Time", ylab="BD Load (ln)", xaxt="n", main = paste0(sp))
    axis(side=1, at = seq(minT,maxT), labels = seq(minT,maxT), cex.axis=0.8)
    for ( r in 1:nrow(tableTemp) ) {
        points(log(tableTemp[r,]+1), type = "l", col=rainbowCol[r], lty=pull(controls[controls$indiv==r,"lty"]))
    }
    abline(v = 6, lty = 2, col = "black")
}
# Last panel is the legend
plot(0,0,pch="",axes=FALSE,xlab="",ylab="")
legend("left"
       , legend= c(sapply(seq(1,9),function(x) paste0("Individual ",x)))
       , lty = 1
       , col = c(rainbowCol, NA)
       , bty = "n"
       , cex=0.75)
dev.off()


for ( sp in speciesList ) {
    pdf(paste0("InfectionOverTime_all_Mean_",sp,".pdf"), height = 5, width=9)
    par(mfrow=c(2,5))
    tableTemp <- bdInfectTime[[sp]][["Mean"]]
    minT <- 1
    maxT <- ncol(tableTemp)
    controls <- which_con %>%
        filter(species==sp) %>%
        mutate(indiv=as.numeric(indiv)) %>%
        arrange(indiv) %>%
        dplyr::select(BD_infected, indiv) %>%
        transmute(lty = if_else(BD_infected=="y","CON","EXP_BD"), indiv= indiv)
    for ( r in 1:nrow(tableTemp) ) {    
        plot(NULL, ylim=c(0,maxLoadLog), xlim=c(minT,maxT),xlab="Time", ylab="BD Load (ln)", xaxt="n", main = paste0("Indiv_",r, "_",controls[controls$indiv==r,"lty"]))
        axis(side=1, at = seq(minT,maxT), labels = seq(minT,maxT), cex.axis=0.8)
        points(log(tableTemp[r,]+1), type = "l", col=rainbowCol[r])
        abline(v = 6, lty = 2, col = "black")
    }
    dev.off()

}



}
########## PLOTTING DATA #########


# Filter out controls for later
# mf.filt.nocon <- mf.tb[-grep("(c|C)ontrol",rownames(mf.tb)),]
# otu.filt.nocon <- otu.filt[, rownames(mf.filt.nocon)]

keep <- which(rownames(dm) %in% colnames(otu.tb.inhib))
dm.filt <- dm[keep,keep]
otu.filt <- otu.tb.inhib[,rownames(dm.filt)]
mf.filt <- mf.tb[match(rownames(dm.filt), mf.tb$SampleID),]

# Set colors for eachs species
colSpecies <- cbind(colors = c("red"
                               , "blue"
                               ,"yellow"
                               ,"purple"
                               ,"green"
                               ,"burlywood4"
                               ,"black"
                               ,"darkgreen"
                               ,"lightblue")
                    , type = levels(factor(mf.tb[,"SampleType"]
                                           , levels = c("Bufo Boreas"
                                                        ,"Bufo Marinus"
                                                        ,"Osteopilus septentrionalis"
                                                        ,"Rana pipiens"
                                                        ,"Rana Catesbeiana"
                                                        ,"ControlCricket"
                                                        ,"ControlHoltfreter"
                                                        ,"ControlMold"
                                                        ,"ControlSterile Water")
                    )))

#### PLOTTING NMDS ####
{# MAKE NMDS
    # Make NMDS of distance matrix
    set.seed(1017937)
    nmds.all <- isoMDS(dist(dm.filt), k = 2)
    nmds <- nmds.all$points
    # dimensions for all plots
    maxlim <- max(c(abs(min(nmds)),max(nmds)))
    }


### Isolate by SPECIES, TIMEPOINT, TREATMENT
{
splitNMDS <- list()
splitDM <- list()
for ( sp in speciesList ) {
    splitNMDS[[sp]] <- list()
    splitDM[[sp]] <- list()
    # First, get overall NMDS for that species ONLY. This is so we can plot all at once
    splitNMDS[[sp]][["all"]] <- nmds[grep(sp, rownames(nmds)),]
    splitDM[[sp]][["all"]] <- dm.filt[grep(sp,rownames(dm.filt)),grep(sp,colnames(dm.filt))]
    # Now, separate out by time point so we can make a gif :)
    for ( t in timeList ) { # SPLIT BY TIME POINTS
        splitNMDS[[sp]][[t]] <- list()
        splitDM[[sp]][[t]] <- list()
        
        rowsTemp <- mf.filt$SampleID %in% rownames(splitNMDS[[sp]][["all"]])
        rowsTempDM <- mf.filt$SampleID %in% rownames(splitDM[[sp]][["all"]])
        timesTemp <- mf.filt$time %in% t
        for ( treat in c("y","n") ) { # WITHIN TIME POINT, SPLIT BY TREAT GROUP
            treatTemp <- mf.filt$BD_infected %in% treat
            # Find intersect
            tokeep <- rowsTemp & timesTemp & treatTemp
            tokeepDM <- rowsTempDM & timesTemp & treatTemp

            samplesToKeep <- pull(mf.filt[tokeep,"SampleID"])
            samplesToKeepDM <- pull(mf.filt[tokeepDM,"SampleID"])
            
            splitNMDS[[sp]][[t]][[treat]] <- nmds[samplesToKeep,]
            splitDM[[sp]][[t]][[treat]] <- dm.filt[samplesToKeepDM,samplesToKeepDM]
        }
    }
}
}

#### PLOT ALL INDIVIDUALS IN EACH SPECIES TOGETHER ######
pre <- c("1","2","3","4","5")
post <- c("6","7","8","9","10","11","12","13","14","15","16")
{
pdf("NMDS_byspecies_treat.pdf",7,5)
par(mfrow=c(2,3), mar = c(4.1,4.1,3.1,1))
for ( sp in speciesList ) {
    plot(NULL, xlim=c(-maxlim,maxlim),ylim=c(-maxlim,maxlim)
         , xlab = "NMDS1", ylab = "NMDS2"
         , main = sp)
    for ( time in timeList ) {
        for ( treat in c("y","n")) {
            if ( time %in% pre ) {
                if ( treat == "y" ) {
                    coltemp <- "lightblue"
                } else if (treat == "n" ) {
                    coltemp <- "blue"
                } else { # should never happen
                    coltemp <- "black"
                }
            } else if ( time %in% post ) {
                if ( treat == "y" ) {
                    coltemp <- "red"
                } else if (treat == "n" ) {
                    coltemp <- "pink"
                } else { # should never happen
                    coltemp <- "black"
                }
            } else { # should never happen
                coltemp <- "black"
            }
            points(splitNMDS[[sp]][[time]][[treat]], bg = coltemp, col= "black", pch = 21)
        }
    }
}
plot(NULL, xlim=c(-maxlim,maxlim),ylim=c(-maxlim,maxlim),bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
legend("center", legend=c("Pre-BD: Control", "Pre-BD: Treatment", "Post-BD: Control", "Post-BD: Treatment"), pch = 21, pt.bg = c("lightblue","blue","pink","red"), col = c("black","black","black","black"))
dev.off()}


##### PLOT EACH TIMEPOINT SEPARATELY ######
dir.create("SEP_by_TIME")
for ( time in timeList ) {
    pdf(paste0("SEP_by_TIME/NMDS_bytime_speciestreat",time,".pdf"),7,5)
    par(mfrow=c(2,3), mar = c(4.1,4.1,3.1,1))
    for ( sp in speciesList ) {
        plot(NULL, xlim=c(-maxlim,maxlim),ylim=c(-maxlim,maxlim)
             , xlab = "NMDS1", ylab = "NMDS2"
             , main = sp)
        for ( treat in c("y","n")) {
            if ( time %in% pre ) {
                if ( treat == "y" ) {
                    coltemp <- "lightblue"
                } else if (treat == "n" ) {
                    coltemp <- "blue"
                } else { # should never happen
                    coltemp <- "black"
                }
            } else if ( time %in% post ) {
                if ( treat == "y" ) {
                    coltemp <- "red"
                } else if (treat == "n" ) {
                    coltemp <- "pink"
                } else { # should never happen
                    coltemp <- "black"
                }
            } else { # should never happen
                coltemp <- "black"
            }
            sampleTemp <- rownames(splitNMDS[[sp]][[time]][[treat]])
            
            outlinecol <- sapply(rowSums(mf.filt[match(sampleTemp,mf.filt$SampleID),c("BD1","BD2","BD3")]) > 0, function(x) {
                if (x) {"yellow"}
                else {"black"}
            })
            points(splitNMDS[[sp]][[time]][[treat]], bg = coltemp, col= outlinecol, pch = 21)
        }
    }
    plot(NULL, xlim=c(-maxlim,maxlim),ylim=c(-maxlim,maxlim),bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
    legend("center", legend=c("Pre-BD: Control", "Pre-BD: Treatment", "Post-BD: Control", "Post-BD: Treatment"), pch = 21, pt.bg = c("lightblue","blue","pink","red"), col = c("black","black","black","black"))
    dev.off()
}

##### CALCULATE NMDS SEPARATELY FOR EACH SPECIES AND PLOT ######
{
    nmds.sep.list <- list()
    pdf("NMDS_byspecies_treat_sepNMDS.pdf",7,5)
par(mfrow=c(2,3), mar = c(4.1,4.1,3.1,1))
for ( sp in speciesList ) {
    nmds.sep.list[[sp]] <- isoMDS(dist(splitDM[[sp]][["all"]]), k=2)
    maxlim <- max(abs(nmds.sep.list[[sp]]$points))
    plot(NULL, xlim=c(-maxlim,maxlim),ylim=c(-maxlim,maxlim)
         , xlab = "NMDS1", ylab = "NMDS2"
         , main = sp)
    # Make colors
    samplesTemp <- rownames(nmds.sep.list[[sp]]$points)
    colorsTemp <- vector(length=length(samplesTemp))
    vec.pre <- mf.filt$prepost[match(samplesTemp,mf.filt$SampleID)]
    prePos <- vec.pre %in% "Pre" 
    # posPos <- grep("pos", samplesTemp)
    vec.bd <- mf.filt$BD_infected[match(samplesTemp,mf.filt$SampleID)]
    bdPos <- vec.bd %in% "y" 
    # nobdpos <- grep("n",vec.bd)
    
    # Now assign colors
    colorsTemp[prePos & bdPos] <- "blue"
    colorsTemp[prePos & !bdPos] <- "lightblue"
    colorsTemp[!prePos & bdPos] <- "red"
    colorsTemp[!prePos & !bdPos] <- "pink"
    colorsTemp[is.na(colorsTemp)] <- "black" # Should never happen
    
    points(nmds.sep.list[[sp]]$points, bg = colorsTemp, col= "black", pch = 21)
}
plot(NULL, xlim=c(-maxlim,maxlim),ylim=c(-maxlim,maxlim),bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
legend("center", legend=c("Pre-BD: Control", "Pre-BD: Treatment", "Post-BD: Control", "Post-BD: Treatment"), pch = 21, pt.bg = c("lightblue","blue","pink","red"), col = c("black","black","black","black"))
dev.off()
}

##### PLOT BETA DISPERSION OVER TIME FOR EACH SPECIES #####
{
    # Set up empty list
    bdisp <- list()
    for ( sp in speciesList ) {
            bdisp[[sp]] <- list() # need two; one for control, one for treat
            # split dm for control and treat
            samplesTemp <- rownames(splitDM[[sp]][["all"]])
            vec.treat <- mf.filt$BD_infected[match(samplesTemp, mf.filt$SampleID)]
            treatdmTemp <- splitDM[[sp]][["all"]][vec.treat == "y", vec.treat == "y"]
            condmTemp <- splitDM[[sp]][["all"]][vec.treat == "n", vec.treat == "n"]
            alldmTemp <- splitDM[[sp]][["all"]]
            # now, calculate betadisper over time
            samplesTemptreat <- rownames(treatdmTemp)
            samplesTempcon <- rownames(condmTemp)
            samplesTempall <- rownames(alldmTemp)
            
            vec.time.treat <- mf.filt$time[match(samplesTemptreat, mf.filt$SampleID)]
            vec.time.con <- mf.filt$time[match(samplesTempcon, mf.filt$SampleID)]
            vec.time.all <- mf.filt$time[match(samplesTempall, mf.filt$SampleID)]
            
            bdisp[[sp]][["treat"]] <- betadisper(dist(treatdmTemp), group = as.factor(vec.time.treat), type="centroid")
            bdisp[[sp]][["con"]] <- betadisper(dist(condmTemp), group = as.factor(vec.time.con), type="centroid")
            bdisp[[sp]][["all"]] <- betadisper(dist(alldmTemp), group = as.factor(vec.time.all), type="centroid")
            
        
    }
    
    ##### Now, plot all betadispersions
    pdf("betadisp_bytime.pdf", height = 5, width = 7)
    par(mfrow = c(2,3))
    for ( sp in speciesList ) {
        if ( !is.null(bdisp[[sp]][["treat"]]) ) {
            means_sd_treat <- aggregate(bdisp[[sp]][["treat"]]$distances, by=list(bdisp[[sp]][["treat"]]$group), FUN=function(x)(c(mean=mean(x),sd=sd(x))))
            means_sd_con <- aggregate(bdisp[[sp]][["con"]]$distances, by=list(bdisp[[sp]][["con"]]$group), FUN=function(x)(c(mean=mean(x),sd=sd(x))))
            
            maxval <- max(c(means_sd_treat[,2][,1] + means_sd_treat[,2][,2],means_sd_con[,2][,1] + means_sd_con[,2][,2]), na.rm = TRUE)
            plot(NULL, xlim=range(as.numeric(timeList)), ylim = c(0,maxval)
                 , ylab = "Distance to centroid"
                 , xlab = "Time Point"
                 , main=sp)
            abline(v=6, lty=2, col="grey")
            #treatment lines
            lines(means_sd_treat[,2][,1], col="red")
            arrows(1:nrow(means_sd_treat), means_sd_treat[,2][,1]-means_sd_treat[,2][,2], 1:nrow(means_sd_treat), means_sd_treat[,2][,1]+means_sd_treat[,2][,2], length=0.05, angle=90, code=3, col="red")
            #control lines
            lines(means_sd_con[,2][,1], col="blue")
            arrows(1:nrow(means_sd_con), means_sd_con[,2][,1]-means_sd_con[,2][,2], 1:nrow(means_sd_con), means_sd_con[,2][,1]+means_sd_con[,2][,2], length=0.05, angle=90, code=3, col="blue")
        }
    }
    plot(NULL, xlim=c(0,1), ylim=c(0,1), xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    legend("center", legend=c("Control","Bd-Treated"), lty=c(1,1), col=c("blue","red"))
    dev.off()
    
    ##### PLOT BETA DISPERSION OVER TIME FOR EACH INDIVIDUAL; INFECTION HIGHLIGHTED #####
    
    ##### Now, plot all betadispersions
    pdf("betadisp_bytime_individual.pdf", height=5, width=7)
    par(mfrow = c(2,3))
    for ( sp in speciesList ) {
        # get vector of individuals and list of unique individuals
        allDisttrecon <- c(bdisp[[sp]][["treat"]]$distances,bdisp[[sp]][["con"]]$distances)
        allGrouptrecon <- c(bdisp[[sp]][["treat"]]$group,bdisp[[sp]][["con"]]$group)
        # Get maxval
        maxval <- max(allDisttrecon, na.rm=TRUE)
        
        # now we can plot all individual lines
        plot(NULL, xlim=range(as.numeric(timeList)), ylim = c(0,maxval)
             , ylab = "Distance to centroid"
             , xlab = "Time Point"
             , main=sp)
        abline(v=6, lty=2, col="grey")
        vec.indiv <- as.vector(mf.filt$toadID[match(c(names(allDisttrecon)), mf.filt$SampleID)])
        unique.indiv <- unique(vec.indiv)
        for ( indiv in unique.indiv ) { # iterate through each individual and aggregate by time
            # Find out if individual is treat or control
            trecon <- as.character(mf.filt$BD_infected[match(indiv, mf.filt$toadID)])
            if (trecon == "y") {
                colTemp <- "red"
            } else if (trecon == "n") {
                colTemp <- "blue"
            } else { colTemp <- "black"} # should never happen
            # Find out vector of infected for this individual
            indivDistTemp <- allDisttrecon[(vec.indiv %in% indiv)]
            infectCol <- rep(colTemp, length(indivDistTemp))
            
            infectedPos <- which(rowSums(mf.filt[mf.filt$SampleID %in% names(indivDistTemp),c("BD1","BD2","BD3")])>1)
            if (length(infectedPos) > 0) {
                infectCol[infectedPos] <- "yellow"
            }
            indivDistTempagg <- aggregate(indivDistTemp, by=list(allGrouptrecon[(vec.indiv %in% indiv)]), FUN = mean)
            infectCol <- aggregate(infectCol, by=list(allGrouptrecon[(vec.indiv %in% indiv)]), FUN = function(x) {x[1]})
            
            lines(indivDistTempagg$x, col=colTemp, lty = 1)
            points(indivDistTempagg$x, col=as.character(infectCol$x), pch = 19)
        }
    }
    plot(NULL, xlim=c(0,1), ylim=c(0,1), xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    legend("center", legend=c("Control","Bd-Treated", "BD detected"), lty=c(1,1,NA), col=c("blue","red","yellow"), pch=c(19,19,19))
    dev.off()
    
    for ( sp in speciesList ) {
    pdf(paste0("betadisp_bytime_individual_",sp,".pdf"), height=5, width=9)
    par(mfrow = c(2,5))
        # get vector of individuals and list of unique individuals
        allDisttrecon <- c(bdisp[[sp]][["all"]]$distances)
        allGrouptrecon <- c(bdisp[[sp]][["all"]]$group)
        # Get maxval
        maxval <- max(allDisttrecon, na.rm=TRUE)
        
        # now we can plot all individual lines
        vec.indiv <- as.vector(mf.filt$toadID[match(c(names(allDisttrecon)), mf.filt$SampleID)])
        unique.indiv <- unique(vec.indiv)
        for ( indiv in unique.indiv ) { # iterate through each individual and aggregate by time
            plot(NULL, xlim=range(as.numeric(timeList)), ylim = c(0,maxval)
                 , ylab = "Distance to centroid"
                 , xlab = "Time Point"
                 , main=sp)
            abline(v=6, lty=2, col="grey")
            # Find out if individual is treat or control
            trecon <- as.character(mf.filt$BD_infected[match(indiv, mf.filt$toadID)])
            if (trecon == "y") {
                colTemp <- "red"
            } else if (trecon == "n") {
                colTemp <- "blue"
            } else { colTemp <- "black"} # should never happen
            # Find out vector of infected for this individual
            indivDistTemp <- allDisttrecon[(vec.indiv %in% indiv)]
            infectCol <- rep(colTemp, length(indivDistTemp))
            infectedPos <- which(rowSums(mf.filt[mf.filt$SampleID %in% names(indivDistTemp),c("BD1","BD2","BD3")])>1)
            if (length(infectedPos) > 0) {
                infectCol[infectedPos] <- "yellow"
            }
            indivDistTempagg <- aggregate(indivDistTemp, by=list(allGrouptrecon[(vec.indiv %in% indiv)]), FUN = mean)
            infectCol <- aggregate(infectCol, by=list(allGrouptrecon[(vec.indiv %in% indiv)]), FUN = function(x) {x[1]})
            
            lines(indivDistTempagg$x, col=colTemp, lty = 1)
            points(indivDistTempagg$x, col=as.character(infectCol$x), pch = 19)
        }
        dev.off()
    }
}

##### PLOT ALPHA DIVERSITY FOR EACH TREATMENT #####
{
    #### WORKING: PLOT INDIIVDUAL ALPHA DIVERSITY; THEN CORRELATE ALPHA/BETA INDIVDIUALLY
    # Set up empty list
    alphadiv <- list()
    for ( sp in speciesList ) {
        alphadiv[[sp]] <- list() # need two; one for control, one for treat
        for ( treat in c("y","n") ) {
            tokeepSP <- mf.filt$species == sp
            tokeepTREAT <- mf.filt$BD_infected == treat
            samplesTemp <- mf.filt$SampleID[tokeepSP & tokeepTREAT]

            alphadiv[[sp]][[treat]] <- mf.filt[match(samplesTemp,mf.filt$SampleID),c("shannon","observed_otus","time", "SampleID")]
        }
    }
    ##### Now, plot all alpha divs
    # have to do each metric
    metrics <- c("shannon","observed_otus")
    for ( m in metrics ) {
        mName <- gsub("_even.*","",m)
        pdf(paste0("alphadiv_bytime_",m,".pdf"), height=5, width=7)
        par(mfrow = c(2,3))
        for ( sp in speciesList ) {
            a_con <- alphadiv[[sp]][["n"]]
            a_tre <- alphadiv[[sp]][["y"]]

            # Find maximum value to plot
            means_sd_con <- aggregate(as.numeric(a_con[[m]]), by=list(a_con[["time"]]), FUN=function(x){c(mean=mean(x),sd=sd(x))})
            means_sd_treat <- aggregate(as.numeric(a_tre[[m]]), by=list(a_tre[["time"]]), FUN=function(x){c(mean=mean(x),sd=sd(x))})
            maxval <- max(c(means_sd_treat[,2][,1] + means_sd_treat[,2][,2],means_sd_con[,2][,1] + means_sd_con[,2][,2]), na.rm = TRUE)
            plot(NULL, xlim=range(as.numeric(timeList)), ylim = c(0,maxval)
                 , ylab = paste0("Alpha diversity (",mName,")")
                 , xlab = "Time Point"
                 , main=sp)
            abline(v=6, lty=2, col="grey")
                #treatment lines
                lines(means_sd_treat[,2][,1], col="red")
                arrows(1:nrow(means_sd_treat), means_sd_treat[,2][,1]-means_sd_treat[,2][,2], 1:nrow(means_sd_treat), means_sd_treat[,2][,1]+means_sd_treat[,2][,2], length=0.05, angle=90, code=3, col="red")
                #control lines
                lines(means_sd_con[,2][,1], col="blue")
                arrows(1:nrow(means_sd_con), means_sd_con[,2][,1]-means_sd_con[,2][,2], 1:nrow(means_sd_con), means_sd_con[,2][,1]+means_sd_con[,2][,2], length=0.05, angle=90, code=3, col="blue")
                abline(v=6, lty=2, col="grey")

        }
        plot(NULL, xlim=c(0,1), ylim=c(0,1), xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
        legend("center", legend=c("Control","Bd-Treated"), lty=c(1,1), col=c("blue","red"))
        dev.off()
    }


    ##### PLOT ALPHA DIV OVER TIME FOR EACH INDIVIDUAL; INFECTION HIGHLIGHTED #####

    ##### Now, plot all alphadivs
    for ( m in metrics ) {
        mName <- gsub("_even.*","",m)
        pdf(paste0("alphadiv_bytime_individual_",mName,".pdf"), height=5, width=7)
        par(mfrow = c(2,3))
        for ( sp in speciesList ) {
            # get vector of individuals and list of unique individuals
            a_con <- alphadiv[[sp]][["n"]][c(m,"SampleID")]
            a_tre <- alphadiv[[sp]][["y"]][c(m,"SampleID")]

            all <- t(rbind((a_con[c(m,"SampleID")]),(a_tre[c(m,"SampleID")])))
            allgroups <- c(alphadiv[[sp]][["n"]]$time,alphadiv[[sp]][["y"]]$time)
            maxval <- max(as.numeric(all), na.rm=TRUE)

            # now we can plot all individual lines
            plot(NULL, xlim=range(as.numeric(timeList)), ylim = c(0,maxval)
                 , ylab = paste0("Alpha Diversity (",mName,")")
                 , xlab = "Time Point"
                 , main=sp)
            abline(v=6, lty=2, col="grey")

            vec.indiv <- as.vector(mf.filt$toadID[match(all["SampleID",], mf.filt$SampleID)])
            unique.indiv <- unique(vec.indiv)
            for ( indiv in unique.indiv ) { # iterate through each individual and aggregate by time
                # Find out if individual is treat or control
                trecon <- as.character(mf.filt$BD_infected[match(indiv, mf.filt$toadID)])
                if (trecon == "y") {
                    colTemp <- "red"
                } else if (trecon == "n") {
                    colTemp <- "blue"
                } else { colTemp <- "black"} # should never happen
                # Find out vector of infected for this individual
                indivDistTemp <- as.numeric(all[m,(vec.indiv %in% indiv)])
                names(indivDistTemp) <- all["SampleID",(vec.indiv %in% indiv)]
                infectCol <- rep(colTemp, length(indivDistTemp))
                infectedPos <- which(rowSums(mf.filt[match(names(indivDistTemp),mf.filt$SampleID),c("BD1","BD2","BD3")])>1)
                if (length(infectedPos) > 0) {
                    infectCol[infectedPos] <- "yellow"
                }
                indivDistTempagg <- aggregate(indivDistTemp, by=list(allgroups[(vec.indiv %in% indiv)]), FUN = mean)
                infectCol <- aggregate(infectCol, by=list(allgroups[(vec.indiv %in% indiv)]), FUN = function(x) {x[1]})

                lines(indivDistTempagg$x, col=colTemp, lty = 1)
                points(indivDistTempagg$x, col=as.character(infectCol$x), pch = 19)
            }
        }
        plot(NULL, xlim=c(0,1), ylim=c(0,1), xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
        legend("center", legend=c("Control","Bd-Treated", "BD detected"), lty=c(1,1,NA), col=c("blue","red","yellow"), pch=c(19,19,19))
        dev.off()

    }
}

# ###### Is alpha and beta diversity lower/higher right before, during, or after infection? ######
# 
# {
#     allInfectionCases <- list()
#     for ( sp in speciesList ) {
#         # get all names of samples from this species
#         namesTreatTemp <- rownames(alphadiv[[sp]]$y)
#         namesConTemp <- rownames(alphadiv[[sp]]$n)
#         # Find names that are infected in each
#         namesInfectedTreat <- namesTreatTemp[rowSums(allBd[namesTreatTemp,])>0]
#         namesInfectedCon <- namesConTemp[rowSums(allBd[namesConTemp,])>0]
#         
#         # Now, go through names:
#         allInfectedNames <- c(namesInfectedTreat,namesInfectedCon )
#         infectionCases <- matrix(nrow = length(allInfectedNames), ncol=14, dimnames = list(allInfectedNames,c("first","last","distBefore","distAfter","alphaBefore","alphaDuring","alphaAfter","firstinfect","only","load"
#                                                                                                               # ,"centrDistBefore","centrDistDuring","centrDistAfter"
#                                                                                                               , "nonAlphamean","nonAlphasd","nonDistmean", "nonDistsd"
#                                                                                                               # ,"noncentrDistmean","noncentrDistsd"
#                                                                                                               )))
#         for ( i in c(namesInfectedTreat,namesInfectedCon ) ){
#             # First, I need to check if this name is a consecutive sample that have remained infected.
#             
#             #IDING THE INDIVIDUAL
#             # What individual are we looking at?
#             indiv <- as.character(mf.filt.nocon.dm[i,"ANONYMIZED_NAME"])
#             # What timepoint was it sampled at?
#             indiv_time <- mf.filt.nocon.dm[i,"timetotal"]
#             # What was mean infection load?
#             load <- mean(allBd[i,])
#             # What are all the other individuals' names?
#             all_indiv <- rownames(mf.filt.nocon.dm)[mf.filt.nocon.dm$ANONYMIZED_NAME==indiv]
#             # Were they infected at all?
#             all_indiv_infect <- rowSums(allBd[all_indiv,])>0
#             # WHich samples were also infected?
#             all_infected_names <- all_indiv[all_indiv_infect]
#             # Which time points were also infected then?
#             all_indiv_infect_time <- mf.filt.nocon.dm[all_infected_names,"timetotal"]
#             
#             
#             # WHAT ARE ITS NEIGHBOURS' STATUS'?
#             # Is this the FIRST sample in this individual that we detected bd?
#             firstInfect <- NA
#             if ( indiv_time == min(all_indiv_infect) ) {
#                 firstInfect <- TRUE
#             } else {
#                 firstInfect <- FALSE
#             }
#             # What are its names right before and after timepoint and was it infected?
#             neighbourpriorName <- rownames(mf.filt.nocon.dm)[(mf.filt.nocon.dm$ANONYMIZED_NAME==indiv) & (mf.filt.nocon.dm$timetotal==c(indiv_time-1))]
#             neighbourpostName <- rownames(mf.filt.nocon.dm)[(mf.filt.nocon.dm$ANONYMIZED_NAME==indiv) & (mf.filt.nocon.dm$timetotal==c(indiv_time+1))]
#             neighbours <- match(c(indiv_time - 1, indiv_time+1),all_indiv_infect_time)
#             # check if it's the first in its clump
#             first <- NA
#             if ( is.na(neighbours[1]) ) {
#                 first <- TRUE
#             } else { first <- FALSE }
#             # check if it's the last in its clump
#             last <- NA
#             if ( is.na(neighbours[2])) {
#                 last <- TRUE
#             } else { last <- FALSE }
#             # check if it's only in its clump
#             only <- NA
#             if ( (first == TRUE) & (last == TRUE) ) {
#                 only <- TRUE
#             } else { only <- FALSE }
#             
#             # WHAT IS ITS TRAVEL DISTANCE FROM LAST POINT
#             distBefore <- NA
#             if ( length(neighbourpriorName) > 0 ) {
#                 distBefore <- splitDM[[sp]][['all']][i,neighbourpriorName]
#             }
#             
#             # WHAT IS ITS TRAVEL DISTANCE IN NEXT STEP
#             distAfter <- NA
#             if ( length(neighbourpostName) > 0 ) {
#                 distAfter <- splitDM[[sp]][['all']][i,neighbourpostName]
#             }
#             # WHAT IS ITS ALPHA DIV PRIOR INFECTION?
#             alphaBefore <- NA
#             if ( length(neighbourpriorName) > 0 ) {
#                 alphaBefore <- mf.filt.nocon.dm[neighbourpriorName,"shannon_even_5000_alpha"]
#             }
#             # WHAT IS ITS ALPHA DIV DURING INFECTION?
#             alphaDuring <- NA
#             alphaDuring <- mf.filt.nocon.dm[i,"shannon_even_5000_alpha"]
#             
#             # WHAT IS ITS ALPHA DIV POST INFECTION?
#             alphaAfter <- NA
#             if ( length(neighbourpostName) > 0 ) {
#                 alphaAfter <- mf.filt.nocon.dm[neighbourpostName,"shannon_even_5000_alpha"]
#             }
#             # # WHAT IS ITS DIST TO CENTROID PRIOR INFECTION?
#             # centrDistBefore <- NA
#             # if ( length(neighbourpriorName) > 0 ) {
#             #     splitDM[[sp]][['all']][]
#             #     centrDistBefore <- mf.filt.nocon.dm[neighbourpriorName,"shannon_even_5000_alpha"]
#             # }
#             # # WHAT IS ITS DIST TO CENTROID DURING INFECTION?
#             # centrDistDuring <- NA
#             # centrDistDuring <- mf.filt.nocon.dm[i,"shannon_even_5000_alpha"]
#             # 
#             # # WHAT IS ITS DIST TO CENTROID POST INFECTION?
#             # centrDistAfter <- NA
#             # if ( length(neighbourpostName) > 0 ) {
#             #     centrDistAfter <- mf.filt.nocon.dm[neighbourpostName,"shannon_even_5000_alpha"]
#             # }
#             # WHAT WAS AVERAGE DISTANCE TRAVELLED AND ALPHA DIV IN NON_INFECTED SAMPLES
#             # First, find out what time points are neight infected nor adjacent to infected individuals
#             time_noninfect_opposite <- unique(c(all_indiv_infect_time + 1, all_indiv_infect_time - 1, all_indiv_infect_time))
#             names_noninfect <- rownames(mf.filt.nocon.dm)[(mf.filt.nocon.dm$ANONYMIZED_NAME == indiv) & is.na(match(mf.filt.nocon.dm$timetotal,time_noninfect_opposite))]
#             nonAlphamean <- mean(as.numeric(mf.filt.nocon.dm[names_noninfect,"shannon_even_5000_alpha"]))
#             nonAlphasd <- sd(as.numeric(mf.filt.nocon.dm[names_noninfect,"shannon_even_5000_alpha"]))
#             nonDistmean <- mean(splitDM[[sp]][['all']][names_noninfect,names_noninfect][lower.tri(splitDM[[sp]][['all']][names_noninfect,names_noninfect])])
#             nonDistsd <- sd(splitDM[[sp]][['all']][names_noninfect,names_noninfect][lower.tri(splitDM[[sp]][['all']][names_noninfect,names_noninfect])])
#             
#             infectionCases[i,c("first"
#                                ,"last"
#                                ,"distBefore"
#                                ,"distAfter"
#                                ,"alphaBefore"
#                                ,"alphaDuring"
#                                ,"alphaAfter"
#                                ,"firstinfect"
#                                ,"only"
#                                ,"load"
#                                # ,"centrDistBefore"
#                                # ,"centrDistDuring"
#                                # ,"centrDistAfter"
#                                , "nonAlphamean","nonAlphasd", "nonDistmean","nonDistsd"
#                                # ,"noncentrDistmean","noncentrDistsd"
#                                )] <- c(first,last,distBefore,distAfter,alphaBefore,alphaDuring,alphaAfter,firstInfect,only,load
#                                        # ,centrDistBefore,centrDistDuring,centrDistAfter
#                                        , nonAlphamean, nonAlphasd, nonDistmean, nonDistsd
#                                        # , noncentrDistmean, noncentrDistsd
#                                        )
#             
#         }
#         
#         allInfectionCases[[sp]] <- infectionCases
#     }
#     
# }
# 
# 
# ###### TESTING
# x <- allInfectionCases$Bubo[,"distBefore"]
# y <- allInfectionCases$Bubo[,"load"]
# con <- allInfectionCases$Bubo
# 
# x <- "distBefore"
# y <- "load"
# # con <- "nonDistmean"
# 
# allDataTest <- matrix(ncol=2)
# for ( sp in speciesList ) {
#     sp <- "Bubo"
#     conTemp <- allInfectionCases[[sp]][,con]
#     xTemp <- allInfectionCases[[sp]][,x]
#     yTemp <- allInfectionCases[[sp]][,y]
#     diff <- (yTemp-conTemp)
#     
#     allDataTestcon <- rbind(allDataTest,cbind(x = xTemp,y = diff))
#     allDataTest <- rbind(allDataTest,cbind(x = xTemp,y = yTemp))
#     
# }
# 
# quartz()
# plot( allDataTestcon )
# 
# quartz()
# plot( allDataTest )




########## GLM ###########


#' The 5 species dataset was collected by Dr. Valerie McKenzie. It looks at how amphibian microbiomes
#' and BD infection influence each other. In the dataset, there are 5 species and 10 individual amphibians with
#' each species (50 total individual amphibians). Skin microbiome samples for each amphibian were sampled across
#' 16 time points. Between time points 5 and 6, half of the individuals from each species were exposed to BD.
#' \
#' The goal of this dataset are to determine:\
#' (1) What factors of the microbiome influence the infection risk and intensity of BD on amphibains across different species?\
#' (2) What factors of the microbiome change after exposure to BD?\
#' \
#' 
#' 
#'
#' #### STEP ZERO: What does the distribution of BD load look like? ####
#' 
#' 
#' Before we build a model, let's look at what our data looks like. Specifically, we want to see whether
#' the distribution of our response variable is normal, lognomral, etc \
#' \
#' \
#' First, let's filter out all the samples that didn't get infected. Also, we will fit some 
#' distributions to the log-transformed BD loads to see how normal it looks.
mf.tb.infected <- mf.tb.inhib %>%
    filter(!(aveBD == 0)) %>%
    mutate(logBD = log(aveBD))
#make x sequence
x <- seq(min(mf.tb.infected$logBD)-sd(mf.tb.infected$logBD),max(mf.tb.infected$logBD)+sd(mf.tb.infected$logBD), length.out = 100)
#fit gamma distribution
fit.gam <- fitdistr(mf.tb.infected$logBD, "Gamma")
y.gam <- dgamma(x,shape=fit.gam$estimate[1], rate=fit.gam$estimate[2])
#fit normal
fit.norm <- fitdistr(mf.tb.infected$logBD, "Normal")
y.norm <- dnorm(x,mean=fit.norm$estimate[1], sd=fit.norm$estimate[2])
#fit poisson
fit.pois <- fitdistr(round(mf.tb.infected$logBD), "Poisson")
y.pois <- dpois(round(x), lambda = fit.pois$estimate[1])
#' Then, we can plot what the distribution of BD looks like with no transformation versus with
#'  a log transformation.
#+ fig.height=8, fig.width=5
pdf("aveBD.pdf")
grid.arrange(ggplot(mf.tb.infected) + 
                 geom_histogram(mapping=aes(x=aveBD, y=..density..), bins=20)
             , mf.tb.inhib %>%
                 filter(!(aveBD == 0)) %>%
                 mutate(logBD = log(aveBD+0.1)) %>%
                 ggplot(aes(x=logBD)) + 
                 geom_histogram(aes(y=..density..), bins=20)  +
                 geom_density(col="blue") +
                 geom_line(data=as.data.frame(x=x, y=y.norm), aes(x,y.norm),col="green")+
                 geom_line(data=as.data.frame(x=x,y=y.gam), aes(x,y.gam), col="red")+
                 geom_line(data=as.data.frame(x=x, y=y.pois), aes(x, y.pois), col="purple")
             , nrow=2)
dev.off()
#' Above, the blue is the kernel density line. Green shows normal. Red shows gamma. After log-transoformation
#' it looks significantly more normal. I think I should note that a poisson distribution might
#' actually be the best fit for this data, but since it has to be counts and my averages are not whole
#' numbers, that I would have to "round" my averages in order to fit a poisson distribution. \
#' \
#' Now, let's see how the histograms look if we split up by species:

#+ echo=FALSE, fig.width=10, fig.height=10
# Fit a few different distributions to the data
gamma.df <- NULL
for ( i in levels(factor(mf.tb.infected$species))) {
    datatemp <- mf.tb.infected %>% 
        filter(species==i) %>%
        dplyr::select(logBD) %>%
        pull()
    x1 <- seq(min(datatemp)-sd(datatemp), max(datatemp)+sd(datatemp), length.out=100)
    fit.temp <- fitdistr(datatemp, "Gamma")
    y1 <- dgamma(x1, fit.temp$estimate[1], fit.temp$estimate[2])
    gamma.df <- as_tibble(rbind(gamma.df, data.frame(x1=x1,y1=y1,species=i)))
}
norm.df <- NULL
for ( i in levels(factor(mf.tb.infected$species))) {
    datatemp <- mf.tb.infected %>% 
        filter(species==i) %>%
        dplyr::select(logBD) %>%
        pull()
    x1 <- seq(min(datatemp)-sd(datatemp), max(datatemp)+sd(datatemp), length.out=100)
    fit.temp <- fitdistr(datatemp, "Normal")
    y1 <- dnorm(x1, fit.temp$estimate[1], fit.temp$estimate[2])
    norm.df <- as_tibble(rbind(norm.df, data.frame(x1=x1,y1=y1,species=i)))
}
pois.df <- NULL
for ( i in levels(factor(mf.tb.infected$species))) {
    datatemp <- mf.tb.infected %>% 
        filter(species==i) %>%
        dplyr::select(logBD) %>%
        mutate(logBD=round(logBD))%>%
        pull()
    x1 <- seq(min(datatemp)-sd(datatemp), max(datatemp)+sd(datatemp), length.out=100)
    fit.temp <- fitdistr(datatemp, "Poisson")
    y1 <- dnorm(x1, fit.temp$estimate[1])
    pois.df <- as_tibble(rbind(pois.df, data.frame(x1=x1,y1=y1,species=i)))
}
mf.tb.inhib %>%
    filter(!(aveBD == 0)) %>%
    mutate(logBD = log(aveBD)) %>%
    ggplot(aes(x=logBD)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_density(col="blue") +
    geom_line(data=gamma.df, aes(x=x1,y=y1), col="red") +
    geom_line(data=norm.df, aes(x=x1, y=y1), col="green") +
    geom_line(data=pois.df, aes(x=x1, y=y1), col="purple") +
    facet_wrap(facets = ~species)

#' Given the plots above, I think a normal distribution is probably sufficient. (Here, purple is poisson) \

#' ### Draw experiment scales and nestedness


#+ fig.width=5, fig.height= 10
# spread data so we can get a visualization of what the grouping looks like
mf.tb.inhib %>%
    dplyr::select(toadID, BD_infected, time, PABD) %>%
    separate(toadID, c("species","Individual")) %>%
    mutate(Individual = factor(Individual, levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))) %>%
    ggplot(aes(x=time, y=Individual)) +
    theme_bw() +
    geom_tile(aes(fill=BD_infected), width=0.8, height=0.8) +
    facet_wrap(~species, nrow=5) +
    geom_vline(aes(xintercept=5.5), color="orange")

mf.tb.inhib %>%
    dplyr::select(toadID, BD_infected, time, PABD, aveBD) %>%
    separate(toadID, c("species","Individual")) %>%
    mutate(Individual = factor(Individual, levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))) %>%
    ggplot(aes(x=time, y=Individual)) +
    theme_bw() +
    geom_tile(aes(fill=log(aveBD)), width=0.8, height=0.8) +
    facet_wrap(~species, nrow=5) +
    geom_vline(aes(xintercept=5.5), color="orange")

#' In the plot above, we see all individuals who were part of the "control" group (not infected with BD) and all individuals who were part of the "treatment" group (infected with BD)
#' #' Also, the orange verticle line indicates where the treatment individuals were first exposed to BD.
#' 

#' ### STEP ONE: Figure out what the "baseline" for each species looks like ####

#+ echo=FALSE, include=FALSE, message=FALSE
# violin
gg_logRich <- mf.tb.inhib %>% 
    filter(prepost=="Pre", !is.na(logRich)) %>%
    separate(toadID, c("species","Identity")) %>%
    mutate(Identity = factor(Identity, levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))) %>%
    ggplot(height=5, width=10) +
    ggtitle("OTU richness")+
    geom_violin(aes(y=logRich, fill=species, x=Identity), show.legend = FALSE) +
    facet_wrap(~species, nrow=1)
gg_shannon <- mf.tb.inhib %>% 
    filter(prepost=="Pre", !is.na(shannon)) %>%
    separate(toadID, c("species","Identity")) %>%
    mutate(Identity = factor(Identity, levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))) %>%
    ggplot(height=5, width=10) +
    ggtitle("OTU diversity")+
    geom_violin(aes(y=shannon, fill=species, x=Identity), show.legend = FALSE) +
    facet_wrap(~species, nrow=1)
gg_distance <- mf.tb.inhib %>% 
    filter(prepost=="Pre", !is.na(distance)) %>%
    separate(toadID, c("species","Identity")) %>%
    mutate(Identity = factor(Identity, levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))) %>%
    ggplot(height=5, width=10) +
    ggtitle("Beta diversity")+
    geom_violin(aes(y=distance, fill=species, x=Identity), show.legend = FALSE) +
    facet_wrap(~species, nrow=1)
gg_inhibRich <- mf.tb.inhib %>% 
    filter(prepost=="Pre", !is.na(inhibRich)) %>%
    separate(toadID, c("species","Identity")) %>%
    mutate(Identity = factor(Identity, levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))) %>%
    ggplot(height=5, width=10) +
    ggtitle("Inhibitory Richness")+
    geom_violin(aes(y=inhibRich, fill=species, x=Identity), show.legend = FALSE) +
    facet_wrap(~species, nrow=1)
gg_percInhib <- mf.tb.inhib %>% 
    filter(prepost=="Pre", !is.na(percInhib)) %>%
    separate(toadID, c("species","Identity")) %>%
    mutate(Identity = factor(Identity, levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))) %>%
    ggplot(height=5, width=10) +
    ggtitle("Inhibitory percent composition of community")+
    geom_violin(aes(y=percInhib, fill=species, x=Identity), show.legend = FALSE) +
    facet_wrap(~species, nrow=1)
#+ fig.height=15, fig.width=10
# grid.arrange(gg_logRich, gg_shannon, gg_distance, gg_inhibRich, gg_percInhib, nrow=5)
grid.arrange(gg_logRich,gg_inhibRich, gg_percInhib, nrow=3)
#' In the plot above, we are plotting ONLY pre-infection toads. ALL toads are shown where diversity measurements
#' were available, although some toads will have dropped out due to low sequencing depth.
#' 
#' ### STEP TWO: Figure out how each of our explanatory variables vary with species ####
# First, group by toadID and pre/post
#+ warnings=FALSE, results="hide", message=FALSE
mf.indiv <- mf.tb.inhib %>%
    group_by(toadID, prepost) %>%
    summarise(rich = mean(logRich, na.rm=TRUE), div = mean(shannon, na.rm=TRUE), dist = mean(distance, na.rm=TRUE)
              , inhibRich = mean(inhibRich, na.rm=TRUE), percInhib = mean(percInhib, na.rm=TRUE)
              , maxBD = max(aveBD, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate(PABD = ifelse(maxBD>0,1,0))  %>%
    gather(variable,value,-c(toadID:prepost)) %>%
    unite(temp, prepost, variable) %>%
    spread(temp, value)

#+ warnings=FALSE, results="hide", message=FALSE
# add in osse 5 because we are missing that
# justOsse <- mf.tb.inhib %>%
    # dplyr::select(toadID, species, BD_infected, toadID, shannon, logRich, distance, inhibRich, percInhib, time) %>%
    # filter(time == 4, species=="Osse")

mf.indiv5 <- mf.tb.inhib %>%
    dplyr::select(toadID, species, BD_infected, toadID, shannon, logRich, distance, inhibRich, percInhib, time) %>%
    arrange(toadID) %>%
    filter(time == 5) %>%
    # rbind(justOsse) %>%
    full_join(mf.indiv) %>%
    separate(toadID, into=c("species","indiv"), remove = FALSE) # some time 5 were missing
# Need to fix BD infected column since that wasn't carried over
mf.indiv5$BD_infected  <- pull(mf.tb.inhib[match(mf.indiv5$toadID, mf.tb.inhib$toadID),"BD_infected"] )

p2 <- mf.indiv5 %>%
    filter(BD_infected=="y") %>%
    ggplot(aes(y=log(Pos_maxBD+0.1)))

#----
# Between species
pdf("infection_btwn_sp.pdf")
p2 + geom_point(aes(x=species, color=species, group=toadID), position=position_jitter(width=0.1, height=0), cex=5)
dev.off()

#'
#' Here, we are ONLY plotting whether or not an intentionally infected individuals. We are NOT showing individuals
#' who were infected but not intentionally, because we can no longer control infection probability.

# So it looks like Bubo is particularily susceptible to BD infection once it takes hold...
#----

# as a factor of community diversity 
# p2 + geom_point(aes(x=Pre_div, color=species), cex=5) +
#     geom_smooth(aes(x=Pre_div), method="lm")
pdf("shannon_bd.pdf")
p2 + geom_point(aes(x=Pre_div, color=species), cex=5) +
    labs(x="Shannon Diversity", y="Log max BD +0.1")
dev.off()
#' above plot shows log max bd + 0.1 infection of all EXPOSED individuals post-exposure
# as a factor of richness
# p2 + geom_point(aes(x=Pre_rich, color=species), cex=5) + 
#     geom_smooth(aes(x=Pre_rich), method="lm")
pdf("otu_bd.pdf")
p2 + geom_point(aes(x=Pre_rich, color=species), cex=5) + 
    labs(x="OTU richness", "Log max BD +0.1")
dev.off()
#----
# as a factor of inhibitory richness or percent of microbiome?
# p2 + geom_point(aes(x=Pre_inhibRich, color=species),position=position_jitter(width=0.1, height=0.1), cex=5) +
#     geom_smooth(aes(x=Pre_inhibRich), method="lm")
pdf("inhibrich_bd.pdf")
p2 + geom_point(aes(x=Pre_inhibRich, color=species),position=position_jitter(width=0.1, height=0.1), cex=5) +
    labs(x="Richness of BD-inhibitory bacteria", y="Log max BD +0.1")
dev.off()

pdf("percrich_bdf.pdf")
p2 + geom_point(aes(x=Pre_percInhib, color=species),position=position_jitter(width=0.1, height=0.1), cex=5) +
    labs(x="Percent BD-inhibitory bacteria (of community)", y="Log max BD +0.1")
dev.off()
# p2 + geom_point(aes(x=Pre_percInhib, color=species),position=position_jitter(width=0.1, height=0.1), cex=5) +
#     geom_smooth(aes(x=Pre_percInhib), method="lm")

# There seems to be a correlation between the richness and proportion of microbiome that is inhibitory, but not actually microiome richness

# What about distance travelled in "multidimensional space" over time?
pdf("dist_bd.pdf")
p2 + geom_point(aes(x=Pre_dist, color=species), cex=5) +
    geom_smooth(aes(x=Pre_dist), method="lm")
dev.off()

##### STEP TWO (B): Look at all individuals who were exposed, infected is categorical and not continuous #####

# jitter <- position_jitter(width=0, height=1)
p2b <- mf.indiv5 %>%
    filter(BD_infected == "y") %>%
    ggplot(aes(y=Pos_PABD))
# diversity
p2b + geom_point(aes(x=Pre_div, color=species), cex=5) +
    geom_smooth(aes(x=Pre_div), method="glm", method.args=list(family="binomial"))
# OTU richness
p2b + geom_point(aes(x=Pre_rich,color=species), cex=5) +
    geom_smooth(aes(x=Pre_rich), method="glm", method.args=list(family="binomial"))
# distance travelled by individual
p2b + geom_point(aes(x=Pre_dist, color=species), cex=5) + 
    geom_smooth(aes(x=Pre_dist), method="glm", method.args=list(family="binomial"))
#richness of inhibitory OTUs
p2b + geom_point(aes(x=Pre_inhibRich, color=species), position=position_jitter(width=0, height=0.05), cex=5) + 
    geom_smooth(aes(x=Pre_inhibRich), method="glm", method.args=list(family="binomial"))
# percent of community of inhibitory OTUs
p2b + geom_point(aes(x=Pre_percInhib, color=species), cex=5) + 
    geom_smooth(aes(x=Pre_percInhib), method="glm", method.args=list(family="binomial"))

#' #### SUMMARY ####
#' It seems that microbiome richness and diversity is not correlated with infection rates or infection intensity. However, 
#' the richness and proportion of BD-inhibitory bacteria seem to be correlated with both infection intensity and rates.\


#' #### Testing with a GLM ####
#' Since BD load is logarithmic, we'll be using a GLM. 
#' Species will be a fixed effect
#' Individual would be a random effect, but there's only one individual per group...\

mf.exposed <- mf.tb.inhib %>%
    filter(prepost == "Pos", BD_infected == "y") %>%
    mutate(aveBD = round(aveBD))

mf.exposed.t5 <- mf.indiv5 %>%
    filter(BD_infected == "y") %>%
    mutate(maxBD = round(Pos_maxBD))

#' Also, here is the first equation for lmer, where j is toadID, S is species, I is inhibitory bacterial richness and/or abundance.\
#' ToadID (random effect) is nested within the species (fixed effect), I think? But that doesn't change the way
#' you input the formula into R. The model equations should remain constant, because glm
#' automatically treats each toadID as unique anyway since they ARE all unique identifiers.\
#' Bi ~ N( uj[i], sigma_i)\
#' uj[i] = beta_s*S + beta_i*I\

#' #### Does the richness or proportion of BD-inhibitory bacteria affect infection PROBABILITY? ####
#' \

# Standardize variables so they are on the same scale
mf.exposed.t5 <- mf.exposed.t5 %>%
    mutate(st_Pre_inhibRich = (Pre_inhibRich-mean(Pre_inhibRich))/sd(Pre_inhibRich)+abs(min((Pre_inhibRich-mean(Pre_inhibRich))/sd(Pre_inhibRich))), st_Pre_percInhib= (Pre_percInhib-mean(Pre_percInhib))/sd(Pre_percInhib)+abs(min((Pre_percInhib-mean(Pre_percInhib))/sd(Pre_percInhib))))

#+ warnings=FALSE, results="hide", message=FALSE, cache=TRUE
glm_toad_inhibRich <- stan_glm(Pos_PABD ~ -1 + species + st_Pre_inhibRich, data=mf.exposed.t5, family=binomial(link="logit"))
#+ warnings=FALSE, results="hide", message=FALSE, cache=TRUE
glm_toad_percInhib <- stan_glm(Pos_PABD ~ -1 + species + st_Pre_percInhib, data=mf.exposed.t5, family=binomial(link="logit"))
#+ warnings=FALSE, results="hide", message=FALSE, cache=TRUE
glm_toad_inhib <- stan_glm(Pos_PABD ~ -1 + species + st_Pre_inhibRich:st_Pre_percInhib, data=mf.exposed.t5, family=binomial(link="logit"))

samps1 <- extract(glm_toad_inhibRich$stanfit)
samps2 <- extract(glm_toad_percInhib$stanfit)
samps3 <- extract(glm_toad_inhib$stanfit)
df1 <- data.frame(samp = as.vector(samps1$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000)
                                ,rep("Osse",4000)
                                ,rep("Raca",4000),rep("Rapi",4000),rep("XInhib",4000))
                  , model=rep("InhibRich",4000*6)
) 
df2 <- data.frame(samp = as.vector(samps2$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("XInhib",4000))
                  , model=rep("percInhib",4000*6)
) 
df3 <- data.frame(samp = as.vector(samps3$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("XInhib",4000))
                  , model=rep("Inhib",4000*6)
) 
overZero<- data.frame(x=c(0,0,0),y=c(500,500,500), label=NA, model=c("InhibRich","percInhib","Inhib"))
# Calculate number of parameters that were larger than 0
overZero[1,3] <- sum(df1 %>%
                         filter(parameter=="XInhib") %>%
                         dplyr::select(samp) %>%
                         pull()
                     >0)
overZero[2,3] <-sum(df2 %>%
                        filter(parameter=="XInhib") %>%
                        dplyr::select(samp) %>%
                        pull()
                    >0)
overZero[3,3] <-sum(df3 %>%
                        filter(parameter=="XInhib") %>%
                        dplyr::select(samp) %>%
                        pull()
                    >0)
overZero[,3] <- overZero[,3]/4000
#+ fig.width=10, fig.height=5
pdf("glm_inhib_on_infection.pdf")
rbind(df1,df2,df3) %>%
    filter(parameter=="XInhib") %>%
    ggplot() +
    geom_histogram(aes(x=samp), bins=20) +
    geom_vline(mapping=aes(xintercept=0), col="red")+
    geom_text(data=overZero, aes(x=x,y=y,label=label), col="blue") +
    facet_wrap(~model, nrow=1, scales="free")
dev.off()
#' It looks like the multiplicative ( nope; inhib rich) model is actually the best model? The samples suggest that there are the fewest
#' parameter estimates that are over 0 than the other two models.
#' Let's look at the model fit for parameter number three:

# Let's look at the model fit
samps_jointinhibmodel <- samps1$beta %>%
    as_tibble() %>%
    rename("Bubo"=V1,"Buma"=V2,"Osse"=V3,"Raca"=V4,"Rapi"=V5,"Inhib"=V6)

# get predictor x vector
inhib_jointeffect <- mf.exposed.t5 %>%
    # transmute(joint_eff = Pre_inhibRich*Pre_percInhib) %>%
    transmute(joint_eff = st_Pre_inhibRich) %>%
    pull()
jointeffect_x <- seq(0, max(inhib_jointeffect)+sd(inhib_jointeffect), length.out = 100)

# run for all species
#Bubo
pred_Bubo <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_jointinhibmodel$Bubo + jointeffect_x[i]*samps_jointinhibmodel$Inhib
    pred_Bubo[,i] <- exp(temp)/(1+exp(temp))
}
summary_Bubo <- data.frame(mean=colMeans(pred_Bubo), upr=apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), lwr=-apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), x=jointeffect_x, species=rep("Bubo",100))
#Buma
pred_Buma <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_jointinhibmodel$Buma + jointeffect_x[i]*samps_jointinhibmodel$Inhib
    pred_Buma[,i] <- exp(temp)/(1+exp(temp))
}
summary_Buma <- data.frame(mean=colMeans(pred_Buma), upr=apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), lwr=-apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), x=jointeffect_x, species=rep("Buma",100))
#Osse
pred_Osse <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_jointinhibmodel$Osse + jointeffect_x[i]*samps_jointinhibmodel$Inhib
    pred_Osse[,i] <- exp(temp)/(1+exp(temp))
}
summary_Osse <- data.frame(mean=colMeans(pred_Osse), upr=apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), lwr=-apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), x=jointeffect_x, species=rep("Osse",100))

#Raca
pred_Raca <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_jointinhibmodel$Raca + jointeffect_x[i]*samps_jointinhibmodel$Inhib
    pred_Raca[,i] <- exp(temp)/(1+exp(temp))
}
summary_Raca <- data.frame(mean=colMeans(pred_Raca), upr=apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), lwr=-apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), x=jointeffect_x, species=rep("Raca",100))

#Rapi
pred_Rapi <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_jointinhibmodel$Rapi + jointeffect_x[i]*samps_jointinhibmodel$Inhib
    pred_Rapi[,i] <- exp(temp)/(1+exp(temp))
}
summary_Rapi <- data.frame(mean=colMeans(pred_Rapi), upr=apply(pred_Rapi,MARGIN=2,sd) + colMeans(pred_Rapi), lwr=-apply(pred_Rapi,MARGIN=2,sd) + colMeans(pred_Rapi), x=jointeffect_x, species=rep("Rapi",100))
# make into one
summary_all <- rbind(summary_Bubo, summary_Buma, summary_Osse, summary_Raca, summary_Rapi)
#+ fig.width=12, fig.height=5
quartz()
ggplot(data=mf.exposed.t5) +
    geom_point(mapping=aes(x=st_Pre_inhibRich, y=Pos_PABD), position=position_jitter(width=0, height=0.05)) +
    geom_line(data=summary_all, aes(x=x, y=mean, color=species)) +
    geom_ribbon(data=summary_all, aes(x=x, ymin=lwr, ymax=upr, bg=species), alpha=0.2) +
    facet_wrap(~species, nrow=1) +
    labs(x="Richness of inhibitory bacteria", y="Presence/Absence BD")

#' It looks like there are certain species who are more likely to get infected than others, and there is a lower chance of infection
#' if your microbiome is RICHER (but not necessarily more proportionally inhibitory bacteria).
#' The effect of a slightly richer microbiome is very subtle, but it is fairly precise\


#' ##### Does richness or proportion of BD-inhibitory bacteria affect INTENSITY of infection? #####

mf.infected.t5 <- mf.indiv5 %>%
    filter(BD_infected == "y", Pos_PABD == 1) %>%
    mutate(logmaxBD = log(Pos_maxBD +0.1), 
           st_Pre_inhibRich = (Pre_inhibRich-mean(Pre_inhibRich))/sd(Pre_inhibRich)+abs(min((Pre_inhibRich-mean(Pre_inhibRich))/sd(Pre_inhibRich))), st_Pre_percInhib= (Pre_percInhib-mean(Pre_percInhib))/sd(Pre_percInhib)+abs(min((Pre_percInhib-mean(Pre_percInhib))/sd(Pre_percInhib))))

#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glm_toad_inhibRich_intensity <- stan_glm(logmaxBD ~ -1 + species + st_Pre_inhibRich, data=mf.infected.t5)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glm_toad_percInhib_intensity <- stan_glm(logmaxBD ~ -1 + species + st_Pre_percInhib, data=mf.infected.t5)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glm_toad_inhib_intensity <- stan_glm(logmaxBD ~ -1 + species + st_Pre_percInhib:st_Pre_inhibRich, data=mf.infected.t5)

samps1 <- extract(glm_toad_inhibRich_intensity$stanfit)
samps2 <- extract(glm_toad_percInhib_intensity$stanfit)
samps3 <- extract(glm_toad_inhib_intensity$stanfit)
df1 <- data.frame(samp = as.vector(samps1$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("XInhib",4000))
                  , model=rep("InhibRich",4000*5)
) 
df2 <- data.frame(samp = as.vector(samps2$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("XInhib",4000))
                  , model=rep("percInhib",4000*5)
) 
df3 <- data.frame(samp = as.vector(samps3$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("XInhib",4000))
                  , model=rep("Inhib",4000*5)
) 
overZero<- data.frame(x=c(0,0,0),y=c(500,500,500), label=NA, model=c("InhibRich","percInhib","Inhib"))
# Calculate number of parameters that were larger than 0
overZero[1,3] <- sum(df1 %>%
                         filter(parameter=="XInhib") %>%
                         dplyr::select(samp) %>%
                         pull()
                     >0)
overZero[2,3] <-sum(df2 %>%
                        filter(parameter=="XInhib") %>%
                        dplyr::select(samp) %>%
                        pull()
                    >0)
overZero[3,3] <-sum(df3 %>%
                        filter(parameter=="XInhib") %>%
                        dplyr::select(samp) %>%
                        pull()
                    >0)
overZero[,3] <- overZero[,3]/4000
#+ fig.height=10, fig.width=10
pdf("glm_inhib_intensityBD.pdf")
rbind(df1,df2,df3) %>%
    filter(parameter=="XInhib") %>%
    as.data.frame() %>%
    ggplot(aes(x=samp)) +
    geom_histogram(aes(y=..count..),bins=50) +
    facet_grid(~model, scale="free")+
    geom_vline(aes(xintercept=0), col="red") +
    geom_text(data=overZero, aes(x=x,y=y,label=label), col="blue")
dev.off()

# Let's look at the model fit
samps_inhibRichmodel <- samps3$beta %>%
    as_tibble() %>%
    rename("Bubo"=V1,"Buma"=V2,"Osse"=V3,"Raca"=V4,"Inhib"=V5)

# get predictor x vector
jointeffect_x <- seq(min(mf.infected.t5$st_Pre_inhibRich)*min(mf.infected.t5$st_Pre_percInhib)
                         , max(mf.infected.t5$st_Pre_inhibRich*mf.infected.t5$st_Pre_percInhib), length.out = 100)

# # run for all species
# #Bubo
# pred_Bubo <- matrix(ncol=100, nrow=4000)
# for ( i in 1:length(jointeffect_x)) {
#     temp <- samps_inhibRichmodel$Bubo + jointeffect_x[i]*samps_inhibRichmodel$Inhib
#     pred_Bubo[,i] <- temp
# }
# summary_Bubo <- data.frame(mean=colMeans(pred_Bubo), upr=apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), lwr=-apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), x=jointeffect_x, species=rep("Bubo",100))
# #Buma
# pred_Buma <- matrix(ncol=100, nrow=4000)
# for ( i in 1:length(jointeffect_x)) {
#     temp <- samps_inhibRichmodel$Buma + jointeffect_x[i]*samps_inhibRichmodel$Inhib
#     pred_Buma[,i] <- temp
# }
# summary_Buma <- data.frame(mean=colMeans(pred_Buma), upr=apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), lwr=-apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), x=jointeffect_x, species=rep("Buma",100))
# #Osse
# pred_Osse <- matrix(ncol=100, nrow=4000)
# for ( i in 1:length(jointeffect_x)) {
#     temp <- samps_inhibRichmodel$Osse + jointeffect_x[i]*samps_inhibRichmodel$Inhib
#     pred_Osse[,i] <- temp
# }
# summary_Osse <- data.frame(mean=colMeans(pred_Osse), upr=apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), lwr=-apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), x=jointeffect_x, species=rep("Osse",100))
# 
# #Raca
# pred_Raca <- matrix(ncol=100, nrow=4000)
# for ( i in 1:length(jointeffect_x)) {
#     temp <- samps_inhibRichmodel$Raca + jointeffect_x[i]*samps_inhibRichmodel$Inhib
#     pred_Raca[,i] <- temp
# }
# summary_Raca <- data.frame(mean=colMeans(pred_Raca), upr=apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), lwr=-apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), x=jointeffect_x, species=rep("Raca",100))
# 
# # make into one
# summary_all <- rbind(summary_Bubo, summary_Buma, summary_Osse, summary_Raca)
# #+ fig.width=12, fig.height=5
# quartz()
# ggplot(data=mf.infected.t5) +
#     geom_point(mapping=aes(x=st_Pre_inhibRich*st_Pre_percInhib, y=logmaxBD), position=position_jitter(width=0, height=0.05)) +
#     geom_line(data=summary_all, aes(x=x, y=mean, color=species)) +
#     geom_ribbon(data=summary_all, aes(x=x, ymin=lwr, ymax=upr, bg=species), alpha=0.2) +
#     facet_wrap(~species, nrow=1) +
#     labs(x="Interaction btwn rich*perc of inhibitory bacteria", y="Log max BD")
# 

#' The relationship between infection intensity (once infected) and inhibitory richness or proportion is actually very weak, at best. Percent
#' inhibitory bacteria is a better predictor than inhibitory richness (which actually suggests infection intensity is HIGHER when inhibitory bacteria are high)
#' but neither are very good predictors. Thus, the intensity of infection seems to be determined almost exclusively by the species
#' of the amphibian (or perhaps by other, unknown predictors? \

#' \
#' I wonder what is a predictor of inhibitory species richness then? I suspect species will be an important indicator of species richness. However, I wonder if 
#' overall diversity or richness may also be a predictor for inhibitory species richness. We also know that proportion of BD inhibitory species tend to decrease
#' in microbiomes over time, so I'll include that as a continuous fixed factor.
#'\
#' First, let's check if inhibitory richness has a normal distirbution
mf.notinfect <- mf.tb.inhib %>%
    filter(BD_infected=="n", PABD == 0) 
# Look at distribution
ggplot(mf.notinfect, aes(x=inhibRich)) +
    geom_histogram(bins=20) +
    facet_grid(~species)

# Looks pretty normal
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
lmer_controls_inhibRich_logRich <- stan_glmer(inhibRich ~ -1 + species + observed_otus + (1|time) + (1|toadID), data=mf.notinfect, family=poisson)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
lmer_controls_inhibRich_div <- stan_glmer(inhibRich ~ -1 + species + shannon + (1|time) + (1|toadID), data=mf.notinfect, family=poisson)


samps1 <- extract(lmer_controls_inhibRich_logRich$stanfit)
samps2 <- extract(lmer_controls_inhibRich_div$stanfit)
df1 <- data.frame(samp = as.vector(samps1$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("Xdiv",4000))
                  , model=rep("rich",4000*6)
) 
byToad1 <- ranef(lmer_controls_inhibRich_logRich)$toadID

df2 <- data.frame(samp = as.vector(samps2$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("Xdiv",4000))
                  , model=rep("div",4000*6)
) 
byToad2 <- ranef(lmer_controls_inhibRich_div)$toadID
byTime2 <- ranef(lmer_controls_inhibRich_div)$time

#+ fig.height=8, fig.width=10
rbind(df1,df2)%>%
    filter(parameter=="Xdiv") %>%
    as.data.frame() %>%
    ggplot(aes(x=samp)) +
    geom_histogram(aes(y=..count..),bins=50) +
    geom_vline(aes(xintercept=0), col="red")+
    facet_grid(~model, scale="free")
# 
# # Let's look at the model fit of the best model
# samps_shannonmodel <- samps2$beta %>%
#     as_tibble() %>%
#     rename("Bubo"=V1,"Buma"=V2,"Osse"=V3,"Raca"=V4,"Rapi"=V5,"shannon"=V6)
# 
# # get predictor x vector
# inhib_effect <- mf.notinfect %>%
#     dplyr::select(shannon) %>%
#     pull()
# jointeffect_x <- seq(max(0,min(inhib_effect,na.rm=TRUE)), max(inhib_effect, na.rm=TRUE), length.out = 100)
# 
# # run for all species
# #Bubo
# pred_Bubo <- matrix(ncol=100, nrow=4000)
# for ( i in 1:length(jointeffect_x)) {
#     temp <- samps_shannonmodel$Bubo + jointeffect_x[i]*samps_shannonmodel$shannon
#     pred_Bubo[,i] <- (temp)
# }
# summary_Bubo <- data.frame(mean=colMeans(pred_Bubo), upr=apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), lwr=-apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), x=jointeffect_x, species=rep("Bubo",100))
# #Buma
# pred_Buma <- matrix(ncol=100, nrow=4000)
# for ( i in 1:length(jointeffect_x)) {
#     temp <- samps_shannonmodel$Buma + jointeffect_x[i]*samps_shannonmodel$shannon
#     pred_Buma[,i] <-(temp)
# }
# summary_Buma <- data.frame(mean=colMeans(pred_Buma), upr=apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), lwr=-apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), x=jointeffect_x, species=rep("Buma",100))
# #Osse
# pred_Osse <- matrix(ncol=100, nrow=4000)
# for ( i in 1:length(jointeffect_x)) {
#     temp <- samps_shannonmodel$Osse + jointeffect_x[i]*samps_shannonmodel$shannon
#     pred_Osse[,i] <- (temp)
# }
# summary_Osse <- data.frame(mean=colMeans(pred_Osse), upr=apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), lwr=-apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), x=jointeffect_x, species=rep("Osse",100))
# #Raca
# pred_Raca <- matrix(ncol=100, nrow=4000)
# for ( i in 1:length(jointeffect_x)) {
#     temp <- samps_shannonmodel$Raca + jointeffect_x[i]*samps_shannonmodel$shannon
#     pred_Raca[,i] <- (temp)
# }
# summary_Raca <- data.frame(mean=colMeans(pred_Raca), upr=apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), lwr=-apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), x=jointeffect_x, species=rep("Raca",100))
# #Rapi
# pred_Rapi <- matrix(ncol=100, nrow=4000)
# for ( i in 1:length(jointeffect_x)) {
#     temp <- samps_shannonmodel$Rapi + jointeffect_x[i]*samps_shannonmodel$shannon
#     pred_Rapi[,i] <- (temp)
# }
# summary_Rapi <- data.frame(mean=colMeans(pred_Rapi), upr=apply(pred_Rapi,MARGIN=2,sd) + colMeans(pred_Rapi), lwr=-apply(pred_Rapi,MARGIN=2,sd) + colMeans(pred_Rapi), x=jointeffect_x, species=rep("Rapi",100))
# 
# # make into one
# summary_all <- rbind(summary_Bubo, summary_Buma, summary_Osse, summary_Raca, summary_Rapi)
# # Adjust the shannon values to account for time and individual
# # mean_Timeeff <- as.vector(colMeans(samps_shannonmodel[,"Time"]))
# mf.notinfect <- mf.notinfect %>%
#     mutate(adjinhibRich = inhibRich - byToad2[toadID,1] - byTime2[time,1])
# #+ fig.width=12, fig.height=5
# # filter out everything except first time point
# ggplot(data=mf.notinfect) +
#     geom_point(mapping=aes(x=shannon, y=adjinhibRich)) +
#     geom_line(data=summary_all, aes(x=x, y=(mean), color=species)) +
#     geom_ribbon(data=summary_all, aes(x=x, ymin=(lwr), ymax=(upr), bg=species), alpha=0.2) +
#     facet_wrap(~species, nrow=1)
# 


#' Overall richness doesn't influence richness of inibitory microbes as much as "diversity"-- 
#' so balanced microbiota that are more even tend to have more inhibitory bacteria. (In addition, time seems to affect diversity slightly as well). In general, however,
#' there is actually not that strong of an effect of either shannon or OTU richness on inhibitory bacteria richness... or at least, it doesn't appear that way.




#' Okay, so finally, I would like to see what kinds of correlations we find between BD infection probability/intensity and microbiome characteristics:
#' eg, do we find that overall microbiome diversity decreases when BD infection load is high?
#' First, let's check normality for all response variables
distance_gg <- mf.exposed %>%
    filter(BD_infected=="y", prepost == "Pos") %>%
    ggplot(aes(x=distance)) +
    geom_histogram(bins=15) +
    facet_grid(~species) +
    ggtitle("Beta diversity")
shannon_gg <- mf.exposed %>%
    filter(BD_infected=="y", prepost == "Pos") %>%
    ggplot(aes(x=shannon)) +
    geom_histogram(bins=15) +
    facet_grid(~species)+
    ggtitle("Shannon diversity")
rich_gg <- mf.exposed %>%
    filter(BD_infected=="y", prepost == "Pos") %>%
    ggplot(aes(x=logRich)) +
    geom_histogram(bins=15) +
    facet_grid(~species)+
    ggtitle("Log richness")
inhibRich_gg <- mf.exposed %>%
    filter(BD_infected=="y", prepost == "Pos") %>%
    ggplot(aes(x=inhibRich)) +
    geom_histogram(bins=15) +
    facet_grid(~species)+
    ggtitle("Inhibitory richness")
percInhib_gg <- mf.exposed %>%
    filter(BD_infected=="y", prepost == "Pos") %>%
    ggplot(aes(x=log(percInhib))) +
    geom_histogram(bins=15) +
    facet_grid(~species)+
    ggtitle("Log percent inhibitory")
#+ fig.height=15, fig.width=10, warnings=FALSE, message=FALSE
grid.arrange(distance_gg, shannon_gg, rich_gg, inhibRich_gg, percInhib_gg, nrow=5)

#' Overall richness and percent inhibitory both need to be lognormal; all others look normal enough I think.

#### Look at how number of inhibitory and percent composition of inhibitory bacteria change between treatments ####

mf.tb.inhib %>%
    mutate(prepost_full = ifelse(prepost=="Pos","Post-exposure","Pre-exposure"), BD_inf_full = ifelse(BD_infected == "y","Exposed","Not Exposed")) %>%
    mutate(prepost_full = factor(prepost_full, levels=c("Pre-exposure","Post-exposure"))) %>%
    ggplot() +
    geom_boxplot(aes(x=prepost_full, y=inhibRich, col=species)) +
    geom_point(aes(x=prepost_full, y=inhibRich, col=species)) +
    facet_grid(species~BD_inf_full) +
    labs(y="Richness of inhibitory bacteria", x="Pre- or Post- exposure")

mf.tb.inhib %>%
    mutate(prepost_full = ifelse(prepost=="Pos","Post-exposure","Pre-exposure"), BD_inf_full = ifelse(BD_infected == "y","Exposed","Not Exposed")) %>%
    mutate(prepost_full = factor(prepost_full, levels=c("Pre-exposure","Post-exposure"))) %>%
    ggplot() +
    geom_boxplot(aes(x=prepost_full, y=percInhib, col=species)) +
    geom_point(aes(x=prepost_full, y=percInhib, col=species)) +
    facet_grid(species~BD_inf_full) +
    labs(y="Percent of OTUs that are BD inhibitory", x="Pre- or Post- exposure")

### make exposed, but not infected ****


###### CORRELATION between diversity and inhibitory bacteria diversity ####
# inhibRich vs shannon
mf.tb.inhib %>%
    mutate(prepost_full = ifelse(prepost=="Pos","Post-exposure","Pre-exposure"), BD_inf_full = ifelse(BD_infected == "y","Exposed","Not Exposed")) %>%
    mutate(prepost_full = factor(prepost_full, levels=c("Pre-exposure","Post-exposure"))) %>%
    ggplot() +
    geom_point(aes(x=shannon, y=inhibRich, col=species)) +
    facet_grid(BD_inf_full ~ prepost_full) +
    labs(x="Shannon diversity", y="Richness of inhibitory bacteria") + 
    geom_smooth(aes(x=shannon, y=inhibRich, col=species), method = "lm")

# inhibRich vs oturich
mf.tb.inhib %>%
    mutate(prepost_full = ifelse(prepost=="Pos","Post-exposure","Pre-exposure"), BD_inf_full = ifelse(BD_infected == "y","Exposed","Not Exposed")) %>%
    mutate(prepost_full = factor(prepost_full, levels=c("Pre-exposure","Post-exposure"))) %>%
    ggplot() +
    geom_point(aes(x=logRich, y=inhibRich, col=species)) +
    facet_grid(BD_inf_full ~ prepost_full) +
    labs(x="Log OTU richness", y="Richness of inhibitory bacteria") + 
    geom_smooth(aes(x=logRich, y=inhibRich, col=species), method = "lm")

# inhibRich vs shannon
mf.tb.inhib %>%
    mutate(prepost_full = ifelse(prepost=="Pos","Post-exposure","Pre-exposure"), BD_inf_full = ifelse(BD_infected == "y","Exposed","Not Exposed")) %>%
    mutate(prepost_full = factor(prepost_full, levels=c("Pre-exposure","Post-exposure"))) %>%
    ggplot() +
    geom_point(aes(x=shannon, y=percInhib, col=species)) +
    facet_grid(BD_inf_full ~ prepost_full) +
    labs(x="Shannon diversity", y="Percent of OTUs that are inhibitory toward BD") + 
    geom_smooth(aes(x=shannon, y=percInhib, col=species), method = "lm")

# inhibRich vs oturich
mf.tb.inhib %>%
    mutate(prepost_full = ifelse(prepost=="Pos","Post-exposure","Pre-exposure"), BD_inf_full = ifelse(BD_infected == "y","Exposed","Not Exposed")) %>%
    mutate(prepost_full = factor(prepost_full, levels=c("Pre-exposure","Post-exposure"))) %>%
    ggplot() +
    geom_point(aes(x=logRich, y=percInhib, col=species)) +
    facet_grid(BD_inf_full ~ prepost_full) +
    labs(x="Log OTU richness", y="Percent of OTUs that are inhibitory toward BD") + 
    geom_smooth(aes(x=logRich, y=percInhib, col=species), method = "lm")


# inhibRich vs percRich
mf.tb.inhib %>%
    mutate(prepost_full = ifelse(prepost=="Pos","Post-exposure","Pre-exposure"), BD_inf_full = ifelse(BD_infected == "y","Exposed","Not Exposed")) %>%
    mutate(prepost_full = factor(prepost_full, levels=c("Pre-exposure","Post-exposure"))) %>%
    ggplot() +
    geom_point(aes(x=inhibRich, y=percInhib, col=species)) +
    facet_grid(BD_inf_full ~ prepost_full) +
    labs(x="Richness of BD-inhibitory bacteria", y="Percent of OTUs that are inhibitory toward BD") + 
    geom_smooth(aes(x=inhibRich, y=percInhib, col=species), method = "lm")


## Looking at the ID of inhibitory bacteria
# Make into relative abundance
otu.tb.inhib %>%
    dplyr::select(-c(OTUID,inhib)) %>%
colSums()
# confirming that it's rarefied to 5000

otu.tb.inhibOnly <- otu.tb.inhib %>%
    filter(inhib=="inhibitory")
otuIDs.inhib <- otu.tb.inhibOnly[,c("OTUID","inhib")]


otu.temp.inhibOnly <- otu.tb.inhibOnly %>%
    dplyr::select(-c(OTUID,inhib)) 
otu.temp.inhibOnly <- otu.temp.inhibOnly/5000

inhib_names <- inhib.tb[match(otuIDs.inhib$OTUID, inhib.tb$seq),"taxa"]
View(cbind(rowSums(otu.temp.inhibOnly), inhib_names))
