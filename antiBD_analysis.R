#!/bin/bash


######## Looking at anti-BD composition of community ##########

otufp='/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/CaptivCol_FR_SS_deblur_t90.final.withtax.txt'
mffp='/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/5Species_best_run_JGK_Modified.txt'
assign_taxafp='/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/ANTIFUNGAL/PICK_CLOSED_OTUS_UCLUST/uclust_ref_picked_otus/node-rep-from-5sp_otus__edited.txt'
bdmffp='/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/ANTIFUNGAL/Amphibian-skin_bacteria_metadata.txt'
thresh = 5 # Change this or lower to '0'

##### LOAD DATA #######
setwd("/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset/ANTIFUNGAL")
otu = read.delim(otufp, header = TRUE, row.names = 1)
mf = read.delim(mffp, header=TRUE, row.names = 1, stringsAsFactors = FALSE)
assign_taxa = read.delim(assign_taxafp,stringsAsFactors = FALSE, row.names=2,na.strings = c("","NA","n/a","na"), header=FALSE)
bdmf = read.delim(bdmffp, row.names = 1, stringsAsFactors = FALSE, header = TRUE)

dir.create("antiBD_R")
setwd("antiBD_R")

##### ADJUSTING DATA: Make same order, make 'time' variable ######
{## First, get rid of "tre"
    mf.sorted = mf[!(mf$TREATMENT_GROUP == "Tre"),]
    
    ### Fix typos in dataset; make a list of all sample types (ie species names)
    # Get the species and control name and get rid of all numbers and underscores
    mf.sorted[,"SampleType"] <- gsub('[0-9]|_','',mf.sorted$SPECIES_NAME)
    # It also turns out there is one called 'Conrol' instead of "Control" so fix that too
    mf.sorted[,"SampleType"] <- gsub("Conrol","Control",mf.sorted[,"SampleType"])
    # Also, "cricket" vs "crickets"
    mf.sorted[,"SampleType"] <- gsub("Crickets","Cricket",mf.sorted[,"SampleType"])
    # Also, typo in sterile water
    mf.sorted[,"SampleType"] <- gsub("Sterlile","Sterile",mf.sorted[,"SampleType"])
    
    
    ## Now, filter OTU Table to include only these samples
    otu.filt = otu[,rownames(mf.sorted)]
    ## Get rid of low abundance OTUs; less than thresh
    otu.filt = otu.filt[rowSums(otu.filt) >= thresh,]
    
    # Make a 'time' variable
    # Take the "PRE_POST_BD_NUM" variable in MF and change it:
    # If the "TREATMENT_GROUP" variable is "Pos", then you add 5 to the time variable
    # If the "TREATMENT_GROUP" variable is "Tre", then you add 16 to the time variable
    # The first 5 timepoints are "Pre" so you don't need to change this.
    # Now, you have a variable that is continuous across the entire treatment
    time <- as.numeric(gsub(".*_","",mf.sorted$PRE_POST_BD_NUM))
    mf.sorted <- cbind(mf.sorted, time)
    # make a 'time2' variable for tracking across all time
    prepos <- mf.sorted$TREATMENT_GROUP
    timetotal <- time
    for ( i in 1:length(prepos) ) {
        if ( prepos[i] == "Pos" ) {
            timetotal[i] <- as.numeric(timetotal[i]) + 5
        } else if ( prepos[i] == "Tre" ) {
            timetotal[i] <- as.numeric(timetotal[i]) + 16
        }
    }
    mf.sorted <- cbind(mf.sorted, timetotal)
    
    # Change "NA"s for BD run to '0'
    # BD runs 1 and 2 have '0's and no NAs, whereas BD run 3 has all NAs and seemingly no zeros.
    # I change them all to zero because I don't actually know which ones are 'true' NAs and which ones are not
    mf.sorted$Bd_Average_Run_3[is.na(mf.sorted$Bd_Average_Run_3)] <- 0
    mf.sorted$Bd_Run_1[is.na(mf.sorted$Bd_Run_1)] <- 0
    mf.sorted$Bd_Run_2[is.na(mf.sorted$Bd_Run_2)] <- 0
}
########### Plot BD load in histogram ############
# Combine all BD values into a single table
Bd1 <- mf.sorted$Bd_Run_1
Bd2 <- mf.sorted$Bd_Run_2
Bd3 <- mf.sorted$Bd_Average_Run_3
allBd <- cbind(Bd1,Bd2,Bd3) # Combine for later
rownames(allBd) <- rownames(mf.sorted)

# Now, chance all zeros to "NA" for this histogram because I can't differentiate
# between a true zero and a "NA" zero
Bd1[Bd1==0] <- NA
Bd2[Bd2==0] <- NA
Bd3[Bd3==0] <- NA

#### Aggregate antifungal stuff ####

{
    assign_taxa_winhib <- cbind(assign_taxa, "inhib"=bdmf[match(assign_taxa[,1], rownames(bdmf)),"Bd_inhibition"])
}


####### split OTU table into different species, individuals ####


# Get list of species
speciesList <- levels(factor(mf.sorted$SPEC))
speciesList <- speciesList[speciesList != "None"]

# Make a time list
timeList <- levels(factor(mf.sorted$timetotal))


# makes list of species over time; save filtered OTU table. save proportion of inhibitory microbes
{
    # In each species, you track individuals (n=6) over time (n=11). value is max or mean BD load for that individual at that time.
    proportion_inhibit <- list()
    split_otu <- list()
    # mf.sorted[,c("maxInfect","meanInfect")] <- matrix(ncol = 2, nrow = nrow(mf.sorted))
    for ( sp in speciesList ) {
        proportion_inhibit[[sp]] <- list()
        otu.temp <- otu.filt[,rownames(mf.sorted)[mf.sorted$SPEC==sp]]
        mf.temp <- mf.sorted[colnames(otu.temp),]
        proportion_inhibit[[sp]][["allotu"]] <- otu.temp
        proportion_inhibit[[sp]][["allmf"]] <- mf.temp
        for ( s in 1:10 ) {
            # sort in order
            mf.temp.temp <- mf.temp[mf.temp$ANONYMIZED_NAME == paste0(sp,"_",s),]
            mf.temp.temp <- mf.temp.temp[with(mf.temp.temp, order(timetotal)),]
            # now get OTU table
            otu.temp.temp <- otu.temp[,rownames(mf.temp.temp)]
            
            proportion_inhibit[[sp]][[paste0(s)]]  <- list()
            proportion_inhibit[[sp]][[paste0(s)]][["otu"]] <- otu.temp.temp
            proportion_inhibit[[sp]][[paste0(s)]][["mf"]] <- mf.temp.temp
        }
    }
    
}

#### Plot the amount of inhibitory bacteria #####

{
    amountInhib <- list()
    for ( sp in speciesList ) {
        amountInhib[[sp]] <- list()
        temp.mat.inhib <- matrix(ncol=16, nrow=10, dimnames = list(as.character(1:10),timeList))
        temp.mat.inhib.prop <- matrix(ncol=16, nrow=10, dimnames = list(as.character(1:10),timeList))
        
        for ( s in 1:10 ) {
            otu.temp <- proportion_inhibit[[sp]][[paste0(s)]][["otu"]]
            mf.temp <- proportion_inhibit[[sp]][[paste0(s)]][["mf"]]
            for ( c in 1:16 ) {
                colTemp <- rownames(mf.temp)[mf.temp[,"timetotal"] == c]
                tempComp <- as.data.frame(otu.temp[,colTemp],row.names = rownames(otu.temp))
                tempComp$inhib <- assign_taxa_winhib$inhib[match(rownames(tempComp), rownames(assign_taxa_winhib))]
                percinhib <- sum(tempComp$`otu.temp[, colTemp]`[which(tempComp$inhib == "inhibitory")])/sum(tempComp$`otu.temp[, colTemp]`)
                percinhibPROP <- sum(tempComp$`otu.temp[, colTemp]`[which(tempComp$inhib == "inhibitory")])/sum(tempComp$`otu.temp[, colTemp]`[which(!is.na(tempComp$inhib))])
                
                temp.mat.inhib[s,c] <- percinhib
                temp.mat.inhib.prop[s,c] <- percinhibPROP
                
            }
        }
        amountInhib[[sp]][["percinhib"]] <- temp.mat.inhib
        amountInhib[[sp]][["percinhibPROP"]] <- temp.mat.inhib.prop
    }
    
    pdf("perc_abx_highlight_bd.pdf",7,5)
    par(mfrow=c(2,3))
    for (sp in speciesList) {
        currenttab <-amountInhib[[sp]][["percinhib"]]
        plot(NULL, ylim=c(0,1), xlim=c(1,16), xlab="Time Point", ylab="Percent inhibitory bact", main=paste0(sp))
        for ( r in 1:nrow(currenttab)) {
            dotCol <- rep("black",16)
            sample_names <- rownames(mf.sorted)[mf.sorted$ANONYMIZED_NAME == paste0(sp,"_",r)]
            toYellow <- rowSums(allBd[sample_names,])>0
            dotCol[toYellow[order(mf.sorted[sample_names,"timetotal"])]] <- "yellow"
            lines(currenttab[r,])
            points(currenttab[r,], pch=19, col=dotCol)
        }
    }
    dev.off()
    
    pdf("perc_abx_prop_highlight_bd.pdf",7,5)
    par(mfrow=c(2,3))
    for (sp in speciesList) {
        currenttab <-amountInhib[[sp]][["percinhibPROP"]]
        plot(NULL, ylim=c(0,1), xlim=c(1,16), xlab="Time Point", ylab="Percent inhibitory bact", main=paste0(sp))
        for ( r in 1:nrow(currenttab)) {
            dotCol <- rep("black",16)
            sample_names <- rownames(mf.sorted)[mf.sorted$ANONYMIZED_NAME == paste0(sp,"_",r)]
            toYellow <- rowSums(allBd[sample_names,])>0
            dotCol[toYellow[order(mf.sorted[sample_names,"timetotal"])]] <- "yellow"
            lines(currenttab[r,])
            points(currenttab[r,], pch=19, col=dotCol)
        }
    }
    dev.off()
    
    for (sp in speciesList) {
        pdf(paste0("perc_abx_highlight_bd_",sp,".pdf"),9,5)
        par(mfrow=c(2,5))
            currenttab <-amountInhib[[sp]][["percinhib"]]
            for ( r in 1:nrow(currenttab)) {
                plot(NULL, ylim=c(0,1), xlim=c(1,16), xlab="Time Point", ylab="Percent inhibitory bact", main=paste0("Indiv_",r))
                dotCol <- rep("black",16)
                sample_names <- rownames(mf.sorted)[mf.sorted$ANONYMIZED_NAME == paste0(sp,"_",r)]
                toYellow <- rowSums(allBd[sample_names,])>0
                dotCol[toYellow[order(mf.sorted[sample_names,"timetotal"])]] <- "yellow"
                lines(currenttab[r,])
                points(currenttab[r,], pch=19, col=dotCol)
            }
        dev.off() 
    }
    
    
}


#### See if BD infection (OR time right before infection) correlates with any bacteria #######

bd.means <- rowMeans(allBd)
corr.otu.bd.TF <- list()
count <- 1
sigSpear <- c()
sigBinom <- c()
for ( o in rownames(otu.filt) ) {
    corr.otu.bd[[o]] <- list()
    abund.temp <- unlist(otu.filt[o,names(bd.means)])
    presabs.temp <- as.numeric(bd.means >0)
    spearTF <- cor.test(log(abund.temp+1), log(bd.means+1), method="spearman")$p.value < 0.05
    binomTF <- summary(glm(presabs.temp ~ abund.temp, family="binomial"))$coefficients[2,4] < 0.05
    
    corr.otu.bd.TF[[o]][["spear"]] <- spearTF
    corr.otu.bd.TF[[o]][["binom"]]  <- binomTF
    
    if ( spearTF ) {
        sigSpear <- c(sigSpear, o)
    }
    
    if ( binomTF ) {
        sigBinom <- c(sigBinom, o)
    }
    if ( count%%100 == 0) {
        print(paste0("done ",count ," out of ", nrow(otu.filt)))
    }
    count <- count + 1
}

# Okay; now that we have all OTUs that are correlated with BD, how many are abx?

sigBinom_inhib <- rownames(assign_taxa_winhib)[which(assign_taxa_winhib[sigBinom,"inhib"] == "inhibitory")]
sigSpear_inhib <- rownames(assign_taxa_winhib)[which(assign_taxa_winhib[sigSpear,"inhib"] == "inhibitory")]

length(sigBinom)
length(sigSpear)


# Crop down the list to include results only that are sign
otu.filt.sigbinom <- otu.filt[sigBinom,]
otu.filt.sigspear <- otu.filt[sigSpear,]

quartz()
par(mfrow=c(1,2))
hist(log(rowSums(otu.filt.sigbinom)))
hist(log(rowSums(otu.filt.sigspear)))


quartz()
abund.temp <- unlist(otu.filt[sigSpear_inhib[3],names(bd.means)])
presabs.temp <- as.numeric(bd.means >0)
plot(log(bd.means+1) ~ log(abund.temp+1))

quartz()
abund.temp <- unlist(otu.filt[sigBinom[4],names(bd.means)])
presabs.temp <- as.numeric(bd.means >0)
plot(log(bd.means+1) ~ log(abund.temp))
