#!/bin/bash
####### Turnover of microbiomes ########

##### PATHWAYS #######
homedir = "/Users/melissachen/Documents/PhD/Amphib_5_sp_dataset"
output = "Turnover"
otuPWD = "MF_and_OTU_edited/otu_table_text.txt"
mfPWD = "MF_and_OTU_edited/MF_withalpha.txt"
thresh = 1 # how many otu observations needed to keep OTU in dataset
indivTHRESH = 5 # BD individual threshold
thirdTHRESH = 50 # BD 3rd sample threshold

dmPWD = "beta_div/bray_curtis_dm.txt"

setwd(homedir)
dir.create(output)
setwd(output)

# Load mappingfile and dm and otutable
mf <- read.delim(paste0("../",mfPWD), header = TRUE, row.names = 1, na.strings = c("na","NA","N/A"))

# Below line is commented out bc it takes a long time to load and I am trying to save time
otu <- read.delim(paste0("../",otuPWD), header=TRUE, row.names = 1, skip=1)

# dm 
dm <- read.delim(paste0("../",dmPWD), header=TRUE, row.names=1)

##### ADJUSTING DATA: Make same order, make 'time' variable ######
## First, get rid of "tre"
mf.filt = mf[!(mf$TREATMENT_GROUP == "Tre"),]

### Fix typos in dataset; make a list of all sample types (ie species names)
# Get the species and control name and get rid of all numbers and underscores
mf.filt[,"SampleType"] <- gsub('[0-9]|_','',mf.filt$SPECIES_NAME)
# It also turns out there is one called 'Conrol' instead of "Control" so fix that too
mf.filt[,"SampleType"] <- gsub("Conrol","Control",mf.filt[,"SampleType"])
# Also, "cricket" vs "crickets"
mf.filt[,"SampleType"] <- gsub("Crickets","Cricket",mf.filt[,"SampleType"])
# Also, typo in sterile water
mf.filt[,"SampleType"] <- gsub("Sterlile","Sterile",mf.filt[,"SampleType"])


## Now, filter OTU Table to include only these samples
otu.filt = otu[,rownames(mf.filt)]
## Get rid of low abundance OTUs; less than thresh
otu.filt = otu.filt[rowSums(otu.filt) >= thresh,]

# Make a 'time' variable
# Take the "PRE_POST_BD_NUM" variable in MF and change it:
# If the "TREATMENT_GROUP" variable is "Pos", then you add 5 to the time variable
# If the "TREATMENT_GROUP" variable is "Tre", then you add 16 to the time variable
# The first 5 timepoints are "Pre" so you don't need to change this.
# Now, you have a variable that is continuous across the entire treatment
time <- as.numeric(gsub(".*_","",mf.filt$PRE_POST_BD_NUM))
mf.filt <- cbind(mf.filt, time)
# make a 'time2' variable for tracking across all time
prepos <- mf.filt$TREATMENT_GROUP
timetotal <- time
for ( i in 1:length(prepos) ) {
    if ( prepos[i] == "Pos" ) {
        timetotal[i] <- as.numeric(timetotal[i]) + 5
    } else if ( prepos[i] == "Tre" ) {
        timetotal[i] <- as.numeric(timetotal[i]) + 16
    }
}
mf.filt <- cbind(mf.filt, timetotal)

# Change "NA"s for BD run to '0'
# BD runs 1 and 2 have '0's and no NAs, whereas BD run 3 has all NAs and seemingly no zeros.
# I change them all to zero because I don't actually know which ones are 'true' NAs and which ones are not
mf.filt$Bd_Average_Run_3[is.na(mf.filt$Bd_Average_Run_3)] <- 0
mf.filt$Bd_Run_1[is.na(mf.filt$Bd_Run_1)] <- 0
mf.filt$Bd_Run_2[is.na(mf.filt$Bd_Run_2)] <- 0

########### Plot BD load in histogram ############
# Combine all BD values into a single table
Bd1 <- mf.filt$Bd_Run_1
Bd2 <- mf.filt$Bd_Run_2
Bd3 <- mf.filt$Bd_Average_Run_3
allBd <- cbind(Bd1,Bd2,Bd3) # Combine for later
rownames(allBd) <- rownames(mf.filt)

####### TURNOVER START ##########
# Make list of all amph individuals
indiv_vec = sort(unique(mf.filt$ANONYMIZED_NAME))
# get rid of "none"
indiv_vec = indiv_vec[-grep("None",indiv_vec)]
# x = indiv_vec[1]
# Now, make a list where the OTU tables are split by individuals.
otu.split.by.indiv = lapply(indiv_vec, function(x) {
    sampNames = rownames(mf.filt)[which(mf.filt$ANONYMIZED_NAME == x)]
    dates = mf.filt[sampNames,"COLLECTION_DATE"]
    dates = gsub("07","2007",dates)
    # dates_numeric = lapply(strsplit(dates,"/"),function(x) {as.numeric(x)})
    dates_formatted = as.Date(dates, format="%m/%d/%Y")
    newOTU = otu.filt[,sampNames[order(dates_formatted)]]
    # Change header names so they are time points instead
    newHeaders = mf.filt$timetotal[match(colnames(newOTU), rownames(mf.filt))]
    colnames(newOTU) = newHeaders
    return(newOTU)
})
names(otu.split.by.indiv) = indiv_vec

# Finally, loop through each individual's progression over time and count how many OTUs are gained or lost.
# Note that this works on a PRESENCE/ABSENCE basis!
# # First, get all time points
# allDates = sort(as.Date(gsub("/07","/2007",unique(mf.filt$COLLECTION_DATE)), format = "%m/%d/%Y"))

turnover = list()
turnover_byOTU = list()
for ( n in 1:(length(indiv_vec)) ) {
    sp = as.character(indiv_vec[n]) # get species name
    currentotutable = otu.split.by.indiv[[as.character(sp)]] # extract current OTU table to use
    # Get temptime list
    # temptime = mf.filt$timetotal[match(colnames(currentotutable),rownames(mf.filt))]
    
    # make empty table to store gained/lost info
    table_changes = matrix(ncol=ncol(currentotutable), nrow = 4, dimnames = list(c("gained","lost","gained_perc","lost_perc"),c(colnames(currentotutable))))
    for ( t in 1:(ncol(currentotutable)-1) ) {
        # t = 1 # for testing
        currentName = colnames(currentotutable)[t]
        current_otus = rownames(currentotutable)[!(currentotutable[,t] == 0)]
        old_total = length(current_otus)
        next_otus = rownames(currentotutable)[!(currentotutable[,t+1] == 0)]
        new_total = length(next_otus)
        
        gained = sum(is.na(match(current_otus, next_otus)))
        lost = sum(is.na(match(next_otus, current_otus)))
        
        table_changes["gained",currentName] = gained
        table_changes["lost",currentName] = lost
        table_changes["gained_perc",currentName] = gained/new_total
        table_changes["lost_perc",currentName] = lost/new_total
    }
    turnover[[sp]] = rbind(table_changes,as.numeric(colnames(currentotutable)))
    
    table_changes_byotu = matrix(ncol=ncol(currentotutable), nrow = 2, dimnames = list(c("constant","constant_perc"),c(colnames(currentotutable))))
    for ( t in 1:(ncol(currentotutable)-1) ) {
        # t = 1 # for testing
        currentName = colnames(currentotutable)[t]
        keept1 = which(!(currentotutable[,t] == 0))
        keept2 = which(!(currentotutable[,t+1] == 0))
        
        otu1 = currentotutable[,t]
        otu2 = currentotutable[,t+1]
        consistent = apply(cbind(otu1,otu2),MARGIN = 1,min) # taking the min between the two, so should be "consistent" OTUs
        
        table_changes_byotu["constant",currentName] = sum(consistent)
        table_changes_byotu["constant_perc",currentName] = sum(consistent)/sum(otu1)
    }
    turnover_byOTU[[sp]] = rbind(table_changes_byotu, as.numeric(colnames(currentotutable)))
    print(paste0("done ",n," out of ", length(indiv_vec)))
}

# Get averages for each species separated into infected vs not infected, with standard error
# FIND CONTROL SPECIES
CON = unique(mf.filt$ANONYMIZED_NAME[mf.filt$BD100KZSP_3122011 == "n"])
TREAT = unique(mf.filt$ANONYMIZED_NAME[mf.filt$BD100KZSP_3122011 == "y"])

allSp = unique(mf.filt$SPEC)
allSp = allSp[-grep("None",allSp)]
aggSpTurnover_lost_perc = list()
aggSpTurnover_gain_perc = list()
aggSpTurnover_constant_perc = list()
for ( i in allSp ) {
    # i = allSp[1]
    conTemp = CON[grep(i, CON)]
    treatTemp = TREAT[grep(i,TREAT)]
    aggSpTurnover_lost_perc[[i]] = list()
    
    for ( type in c("con","treat") ) {
        # type = "con"
        if (type == "con") {
            posSp = match(conTemp,names(turnover))
            tempTable_lost = matrix(ncol=16, nrow = length(posSp), dimnames=list(conTemp,seq(1,16)))
            tempTable_gain = matrix(ncol=16, nrow = length(posSp), dimnames=list(conTemp,seq(1,16)))
            tempTable_constant = matrix(ncol=16, nrow = length(posSp), dimnames=list(conTemp,seq(1,16)))
        } else if (type == "treat" ) {
            posSp = match(treatTemp,names(turnover))
            tempTable_lost = matrix(ncol=16, nrow = length(posSp), dimnames=list(treatTemp,seq(1,16)))
            tempTable_gain = matrix(ncol=16, nrow = length(posSp), dimnames=list(treatTemp,seq(1,16)))
            tempTable_constant = matrix(ncol=16, nrow = length(posSp), dimnames=list(treatTemp,seq(1,16)))
        }
        # posSp = grep(i, names(turnover))
        
        for ( j in posSp ) {
            # j = 34
            for ( k in colnames(turnover[[j]])) {
                # k = 2
                tempTable_lost[names(turnover[j]),paste0(k)] = turnover[[j]]["lost_perc",paste0(k)]
                tempTable_gain[names(turnover[j]),paste0(k)] = turnover[[j]]["gained_perc",paste0(k)]
                tempTable_constant[names(turnover[j]),paste0(k)] = turnover_byOTU[[j]]["constant_perc",paste0(k)]
            }
        }
        aggSpTurnover_lost_perc[[i]][[type]] = tempTable_lost
        aggSpTurnover_gain_perc[[i]][[type]] = tempTable_gain
        aggSpTurnover_constant_perc[[i]][[type]] = tempTable_constant
    }
    
    
}

## Now, let's try plotting

###### PLOTTING TURNOVER #######

### RETAINED PER OTU
pdf("percOTUs_retained.pdf")
par(mfrow=c(3,2),mar=c(4.1,4.1,4.1,1.1))
for ( n in names(aggSpTurnover_constant_perc) ) {
    i = aggSpTurnover_constant_perc[[n]]
    plot(NULL, xlim=c(1,16), ylim=c(0,1), xlab="Time point", ylab="Percent of retained reads", main=n)
    
    # Plot controls first
    # first ,summarize data
    meansTemp = colMeans(i$treat,na.rm = TRUE)
    sdTemp = apply(i$treat,MARGIN = 2,sd)
    x = as.numeric(colnames(i$treat))
    # now, draw dots, lines, etc
    lines(meansTemp ~ as.numeric(colnames(i$treat)), col="red")
    points(meansTemp ~ as.numeric(colnames(i$treat)), pch=19, col="red")
    arrows(x, meansTemp-sdTemp, x, meansTemp+sdTemp, length=0.05, angle=90, code=3, col="red")
    
    # Now, controls
    meansTemp = colMeans(i$con,na.rm = TRUE)
    sdTemp = apply(i$con,MARGIN = 2,sd)
    x = as.numeric(colnames(i$con))
    # now, draw dots, lines, etc
    lines(meansTemp ~ as.numeric(colnames(i$con)), col="blue")
    points(meansTemp ~ as.numeric(colnames(i$con)), pch=19, col="blue")
    arrows(x, meansTemp-sdTemp, x, meansTemp+sdTemp, length=0.05, angle=90, code=3, col="blue")
}
plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt = "n", bty="n")
legend("center", legend=c("Controls","BD-infected"), pch = c(19,19), lty = c(1,1), col=c("Blue","Red"), bty="n")
dev.off()

#### LOST OTUS
pdf("percOTUs_lost.pdf")
par(mfrow=c(3,2),mar=c(4.1,4.1,4.1,1.1))
for ( n in names(aggSpTurnover_lost_perc) ) {
    i = aggSpTurnover_lost_perc[[n]]
    plot(NULL, xlim=c(1,16), ylim=c(0,1), xlab="Time point", ylab="Percent of OTUs lost", main=n)
    
    # Plot controls first
    # first ,summarize data
    meansTemp = colMeans(i$treat,na.rm = TRUE)
    sdTemp = apply(i$treat,MARGIN = 2,sd)
    x = as.numeric(colnames(i$treat))
    # now, draw dots, lines, etc
    lines(meansTemp ~ as.numeric(colnames(i$treat)), col="red")
    points(meansTemp ~ as.numeric(colnames(i$treat)), pch=19, col="red")
    arrows(x, meansTemp-sdTemp, x, meansTemp+sdTemp, length=0.05, angle=90, code=3, col="red")
    
    # Now, controls
    meansTemp = colMeans(i$con,na.rm = TRUE)
    sdTemp = apply(i$con,MARGIN = 2,sd)
    x = as.numeric(colnames(i$con))
    # now, draw dots, lines, etc
    lines(meansTemp ~ as.numeric(colnames(i$con)), col="blue")
    points(meansTemp ~ as.numeric(colnames(i$con)), pch=19, col="blue")
    arrows(x, meansTemp-sdTemp, x, meansTemp+sdTemp, length=0.05, angle=90, code=3, col="blue")
}
plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt = "n", bty="n")
legend("center", legend=c("Controls","BD-infected"), pch = c(19,19), lty = c(1,1), col=c("Blue","Red"), bty="n")
dev.off()

#### GAINED OTUS
pdf("percOTUs_gained.pdf")
par(mfrow=c(3,2),mar=c(4.1,4.1,4.1,1.1))
for ( n in names(aggSpTurnover_gain_perc) ) {
    i = aggSpTurnover_gain_perc[[n]]
    plot(NULL, xlim=c(1,16), ylim=c(0,5), xlab="Time point", ylab="Percent of OTUs gained", main=n)
    
    # Plot controls first
    # first ,summarize data
    meansTemp = colMeans(i$treat,na.rm = TRUE)
    sdTemp = apply(i$treat,MARGIN = 2,sd)
    x = as.numeric(colnames(i$treat))
    # now, draw dots, lines, etc
    lines(meansTemp ~ as.numeric(colnames(i$treat)), col="red")
    points(meansTemp ~ as.numeric(colnames(i$treat)), pch=19, col="red")
    arrows(x, meansTemp-sdTemp, x, meansTemp+sdTemp, length=0.05, angle=90, code=3, col="red")
    
    # Now, controls
    meansTemp = colMeans(i$con,na.rm = TRUE)
    sdTemp = apply(i$con,MARGIN = 2,sd)
    x = as.numeric(colnames(i$con))
    # now, draw dots, lines, etc
    lines(meansTemp ~ as.numeric(colnames(i$con)), col="blue")
    points(meansTemp ~ as.numeric(colnames(i$con)), pch=19, col="blue")
    arrows(x, meansTemp-sdTemp, x, meansTemp+sdTemp, length=0.05, angle=90, code=3, col="blue")
}
plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt = "n", bty="n")
legend("center", legend=c("Controls","BD-infected"), pch = c(19,19), lty = c(1,1), col=c("Blue","Red"), bty="n")
dev.off()

#### histogram version
### consistent
pdf("percOTUs_retained_histogram.pdf")
par(mfrow=c(3,2),mar=c(4.1,4.1,4.1,1.1))
for ( n in names(aggSpTurnover_constant_perc) ) {
    # i = aggSpTurnover_constant_perc[[n]]
    hist(c(aggSpTurnover_constant_perc[[n]]$con,aggSpTurnover_constant_perc[[n]]$treat[,c(1,2,3,4,5)]), col=rgb(0,0,1,0.3), xlim=c(0,1)
         ,xlab="Percent reads retained"
         , main=n)
    hist(aggSpTurnover_constant_perc[[n]]$treat[,c(6,7,8,9,10,11,12,13,14,15,16)], col=rgb(1,0,0,0.3), xlim=c(0,1),add=T)
    
}
plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt = "n", bty="n")
legend("center", legend=c("Controls","BD-infected"), pch = c(19,19), lty = c(1,1), col=c("Blue","Red"), bty="n")
dev.off()

### lost
pdf("percOTUs_lost_histogram.pdf")
par(mfrow=c(3,2),mar=c(4.1,4.1,4.1,1.1))
for ( n in names(aggSpTurnover_lost_perc) ) {
    # i = aggSpTurnover_lost_perc[[n]]
    hist(c(aggSpTurnover_lost_perc[[n]]$con,aggSpTurnover_lost_perc[[n]]$treat[,c(1,2,3,4,5)]), col=rgb(0,0,1,0.3), xlim=c(0,1), freq=FALSE
         ,xlab="Percent reads retained"
         , main=n
         ,ylim=c(0,4))
    hist(aggSpTurnover_lost_perc[[n]]$treat[,c(6,7,8,9,10,11,12,13,14,15,16)], col=rgb(1,0,0,0.3), xlim=c(0,1),freq=FALSE,add=T)
    
}
plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt = "n", bty="n")
legend("center", legend=c("Controls","BD-infected"), pch = c(19,19), lty = c(1,1), col=c("Blue","Red"), bty="n")
dev.off()

### gained
pdf("percOTUs_gained_histogram.pdf")
par(mfrow=c(3,2),mar=c(4.1,4.1,4.1,1.1))
for ( n in names(aggSpTurnover_gain_perc) ) {
    # i = aggSpTurnover_gain_perc[[n]]
    hist(c(aggSpTurnover_gain_perc[[n]]$con,aggSpTurnover_gain_perc[[n]]$treat[,c(1,2,3,4,5)]), col=rgb(0,0,1,0.3), xlim=c(0,1), freq = FALSE
         ,xlab="Percent reads retained"
         , main=n)
    hist(aggSpTurnover_gain_perc[[n]]$treat[,c(6,7,8,9,10,11,12,13,14,15,16)], col=rgb(1,0,0,0.3), xlim=c(0,1),freq=FALSE,add=T)
    
}
plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt = "n", bty="n")
legend("center", legend=c("Controls","BD-infected"), pch = c(19,19), lty = c(1,1), col=c("Blue","Red"), bty="n")
dev.off()
# 
# ####### ALPHA DIVERSITY-- all samples ########
# ## Do all samples first
# # Examine how alpha diversity correlates with infection
# allBd.sorted <- allBd[rownames(mf.filt),]
# allBd.sorted.max <- apply(allBd.sorted, MARGIN=1, max)
# allBd.sorted.presabs <- allBd.sorted.max
# allBd.sorted.presabs[allBd.sorted.presabs > 0] <- 1
# 
# # CHANGE THIS to change alpha metrics
# alphadiv_shannon <- mf.filt$shannon_even_5000_alpha
# names(alphadiv_shannon) <- rownames(mf.filt)
# alphadiv_obs <- mf.filt$observed_otus_even_5000_alpha
# names(alphadiv_obs) <- rownames(mf.filt)
# 
# alpha_bd_shannon <- cbind(alpha = alphadiv_shannon, bd =allBd.sorted.presabs)
# alpha_bd_obs <- cbind(alpha = alphadiv_obs, bd =allBd.sorted.presabs)
# 
# # use a general linear model for binomial data
# glm.alpha.presabs.shannon <- glm(bd ~ alpha, family = binomial(link="logit"), data = as.data.frame(alpha_bd_shannon))
# capture.output(summary(glm.alpha.presabs.shannon), file="glm_binomial_shannon.txt")
# 
# glm.alpha.presabs.obs <- glm(bd ~ alpha, family = binomial(link="logit"), data = as.data.frame(alpha_bd_obs))
# capture.output(summary(glm.alpha.presabs.obs), file="glm_binomial_obs.txt")
# 
# # Plot results
# # Get prediction (probability of infection)- SHANNON
# rangeVal <- range(alphadiv_shannon, na.rm = TRUE)
# xeval <- seq(rangeVal[1], rangeVal[2],length.out = 100)
# yeval <- predict(glm.alpha.presabs.shannon, list(alpha = xeval), type="response")
# 
# pdf("glm_bd_shannon_all.pdf")
# plot(alpha_bd_shannon, xlab="Shannon Diversity", ylab="Probability infection")
# lines(yeval ~ xeval)
# dev.off()
# 
# # Get prediction (probability of infection)- OBSERVED OTUs
# rangeVal <- range(alphadiv_obs, na.rm = TRUE)
# xeval <- seq(rangeVal[1], rangeVal[2],length.out = 100)
# yeval <- predict(glm.alpha.presabs.obs, list(alpha = xeval), type="response")
# 
# pdf("glm_bd_obsOTUs_all.pdf")
# plot(alpha_bd_obs, xlab="Observed OTUs", ylab="Probability infection")
# lines(yeval ~ xeval)
# dev.off()

##### ALPHA DIVERSITY-- filter out pre #####
## Do post samples only
mf.postreatOnly <- mf.filt[(mf.filt$TREATMENT_GROUP == "Pos") & (mf.filt$BD100KZSP_3122011 == "y"),]
# Examine how alpha diversity correlates with infection
allBd.sorted <- allBd[rownames(mf.postreatOnly),]
allBd.sorted.max <- apply(allBd.sorted, MARGIN=1, max)
allBd.sorted.presabs <- allBd.sorted.max
allBd.sorted.presabs[allBd.sorted.presabs > 0] <- 1

# CHANGE THIS to change alpha metrics
alphadiv_shannon <- mf.postreatOnly$shannon_even_5000_normalized_alpha
names(alphadiv_shannon) <- rownames(mf.postreatOnly)
alphadiv_obs <- mf.postreatOnly$observed_otus_even_5000_normalized_alpha
names(alphadiv_obs) <- rownames(mf.postreatOnly)

alpha_bd_shannon <- cbind(alpha = alphadiv_shannon, bd =allBd.sorted.presabs)
alpha_bd_obs <- cbind(alpha = alphadiv_obs, bd =allBd.sorted.presabs)

# use a general linear model for binomial data
glm.alpha.presabs.shannon <- glm(bd ~ alpha, family = binomial(link="logit"), data = as.data.frame(alpha_bd_shannon))
capture.output(summary(glm.alpha.presabs.shannon), file="glm_binomial_shannon.txt")

glm.alpha.presabs.obs <- glm(bd ~ alpha, family = binomial(link="logit"), data = as.data.frame(alpha_bd_obs))
capture.output(summary(glm.alpha.presabs.obs), file="glm_binomial_obs.txt")

# Plot results
# Get prediction (probability of infection)- SHANNON
rangeVal <- range(alphadiv_shannon, na.rm = TRUE)
xeval <- seq(rangeVal[1], rangeVal[2],length.out = 100)
yeval <- predict(glm.alpha.presabs.shannon, list(alpha = xeval), type="response")

pdf("glm_bd_shannon_posOnly.pdf")
plot(alpha_bd_shannon, xlab="Shannon Diversity", ylab="Probability infection")
lines(yeval ~ xeval)
dev.off()

# Get prediction (probability of infection)- OBSERVED OTUs
rangeVal <- range(alphadiv_obs, na.rm = TRUE)
xeval <- seq(rangeVal[1], rangeVal[2],length.out = 100)
yeval <- predict(glm.alpha.presabs.obs, list(alpha = xeval), type="response")

pdf("glm_bd_obsOTUs_posOnly.pdf")
plot(alpha_bd_obs, xlab="Observed OTUs", ylab="Probability infection")
lines(yeval ~ xeval)
dev.off()



##### TRY TO CORRELATE ALPHA DIVERSITY AGAINST STABILITY NEXT ####
# Find all turnovers
turnover_alpha_agg <- cbind(alphadiv_obs,alphadiv_shannon,turnover=NA)
gained_alpha_agg <- cbind(alphadiv_obs,alphadiv_shannon,turnover=NA)
lost_alpha_agg <- cbind(alphadiv_obs,alphadiv_shannon,turnover=NA)
for ( i in rownames(turnover_alpha_agg) ) {
    datatemp <- mf.postreatOnly[i,c("SPEC","timetotal","ANONYMIZED_NAME")]
    turnover_alpha_agg[i,3] <- aggSpTurnover_constant_perc[[as.character(datatemp[1,1])]][["treat"]][as.character(datatemp[1,3]),as.character(datatemp[1,2])] 
    gained_alpha_agg[i,3] <- aggSpTurnover_gain_perc[[as.character(datatemp[1,1])]][["treat"]][as.character(datatemp[1,3]),as.character(datatemp[1,2])] 
    lost_alpha_agg[i,3] <- aggSpTurnover_lost_perc[[as.character(datatemp[1,1])]][["treat"]][as.character(datatemp[1,3]),as.character(datatemp[1,2])] 
    
}

#RETAINED
pdf("turnover_alpha_correlated_obs.pdf")
plot(turnover_alpha_agg[,c(1,3)], xlab="Observed Species", ylab="Retained Reads")
dev.off()
pdf("turnover_alpha_correlated_shannon.pdf")
plot(turnover_alpha_agg[,c(2,3)], xlab="Shannon Diversity", ylab="Retained Reads")
dev.off()
# 
# #gained
# quartz()
# plot(gained_alpha_agg[,c(1,3)], xlab="Observed Species", ylab="Gained OTUs")
# quartz()
# plot(gained_alpha_agg[,c(2,3)], xlab="Shannon Diversity", ylab="Gained OTUs")
# 
# #lost
# quartz()
# plot(lost_alpha_agg[,c(1,3)], xlab="Observed Species", ylab="Lost OTUs")
# quartz()
# plot(lost_alpha_agg[,c(2,3)], xlab="Shannon Diversity", ylab="Lost OTUs")

#### 

#### Are there any OTUs that were missing RIGHT before infection?

#### Is turnover along time linear, or is it "random" around a centroid?
