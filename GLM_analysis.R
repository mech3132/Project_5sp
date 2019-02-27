#' ---
#' title: Community microbiome and BD infection rates
#' author: Melissa Chen
#' output: github_document
#' ---
#' 

#+ message=FALSE, echo=FALSE, warnings=FALSE
###### INPUT VARIABLES AND PATHWAYS ######
library(tidyverse)
library(ggplot2)
library(vegan)
library(rstan)
library(rstanarm)
library(MASS)
library(gridExtra)
library(rethinking)
#+ message=FALSE, echo=FALSE, warnings=FALSE, cache=TRUE
# Mapping file with alpha div
mfPWD <- "../MF_and_OTU_edited/MF_withalpha_osseadded.txt"
# otu table
otuPWD <- "../MF_and_OTU_edited/otu_table_r5000.txt"
# bray curtis distance matrix
dmPWD <- "../beta_div/bray_curtis_dm.txt"
# inhibitory metadata
inhibPWD <- "../ANTIFUNGAL/inhibitory_metadata_MANUAL.txt"
# Minimum reads in otu table for OTU to be retained
minOTUTab <- 100
# Minimum reads in sample for OTU to NOT be turned into zero
minOTUSample <- 5
###### Load data files #######
dm <- read.delim(dmPWD, header=1, row.names=1, as.is=TRUE)
mf <- read.delim(mfPWD, header=1, as.is=TRUE)
otu <- read.delim(otuPWD, skip=1, header=1, as.is=TRUE)
inhib <- read.delim(inhibPWD, header=FALSE, as.is=TRUE)
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
               ,"CONT_NOT_DOSED_BD_UNKNOWNIVE"
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
              , "secondexp_con")
mf.tb <- mf %>%
    as_tibble() %>% # make into tibble
    dplyr::select(one_of(toKeepMF)) %>% # filter to only relevant variables
    rename_at(vars(toKeepMF), ~ newNames) %>%#rename variable names
    filter(prepost == "Pre" | prepost == "Pos") %>% # get rid of things in "tre", which is a different experiment
    separate(timepoint, into=c("exposure","time"), sep="_", convert=TRUE) %>% # create a time variable
    mutate( time = ifelse(exposure == "Pre", time, time + 5)) %>% # time variable "restarts" at BD exposure point, but I want it to be continuous
    filter(species != "None")
mf.tre <- mf %>%
    as_tibble() %>% # make into tibble
    dplyr::select(one_of(toKeepMF)) %>% # filter to only relevant variables
    rename_at(vars(toKeepMF), ~ newNames) %>%#rename variable names
    filter(prepost == "Tre" | prepost == "Tre") %>% # get rid of things in "tre", which is a different experiment
    separate(timepoint, into=c("exposure","time"), sep="_", convert=TRUE) %>% # create a time variable
    mutate( time = ifelse(exposure == "Pre", time, time + 5)) %>% # time variable "restarts" at BD exposure point, but I want it to be continuous
    filter(species != "None")
#+ include=FALSE
mf.tb$shannon <- as.numeric(mf.tb$shannon)
mf.tb$observed_otus <- as.numeric(mf.tb$observed_otus)
mf.tre$shannon <- as.numeric(mf.tre$shannon)
mf.tre$observed_otus <- as.numeric(mf.tre$observed_otus)
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
#Samples to keep in otu table
OTU_names <- otu$X.OTU.ID
otu.tb <- otu %>%
    as_tibble() %>%
    dplyr::select(one_of(c(mf.tb$SampleID))) %>%
    replace(.<minOTUSample, 0) %>%
    mutate(rowsums=rowSums(.)) %>%
    mutate(OTUID = OTU_names) %>%
    filter(rowsums > minOTUTab) %>%
    dplyr::select(-rowsums)
#Samples to keep in otu table for "Tre" experiment
otu.tre <- otu %>%
    as_tibble() %>%
    dplyr::select(one_of(c(mf.tre$SampleID))) %>%
    replace(.<minOTUSample, 0) %>%
    mutate(rowsums=rowSums(.)) %>%
    mutate(OTUID = OTU_names) %>%
    filter(rowsums > minOTUTab) %>%
    dplyr::select(-rowsums)
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
mf.tb[,c("BD1","BD2","BD3")] <- allBd
# Now, make an "average" bd load
mf.tb <- mf.tb %>%
    rowwise() %>%
    mutate(aveBD = mean(c(BD1,BD2,BD3)), aveBD = max(c(BD1,BD2,BD3))) %>%
    # Finally, make BD presence absence
    mutate(PABD = ifelse(aveBD>0, 1, 0), logRich=log(observed_otus))

# Repeat for Tre
allBd <- mf.tre %>% 
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
mf.tre[,c("BD1","BD2","BD3")] <- allBd
# Now, make an "average" bd load
mf.tre <- mf.tre %>%
    rowwise() %>%
    mutate(aveBD = mean(c(BD1,BD2,BD3)), aveBD = max(c(BD1,BD2,BD3))) %>%
    # Finally, make BD presence absence
    mutate(PABD = ifelse(aveBD>0, 1, 0), logRich=log(observed_otus))


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

# Repeat for tre
mf.tre$distance <- NA
for ( i in 1:nrow(mf.tre) ) {
    currentSample <- mf.tre$SampleID[i]
    current <- mf.tre[i, c("toadID","time")]
    prevSample <- mf.tre$SampleID[mf.tre$toadID==as.character(current[1,1]) & mf.tre$time == as.numeric(current[1,2]-1)]
    if ( length(prevSample) > 0) {
        distTemp <- dm[currentSample,as.character(prevSample)]
        if ( length(distTemp) > 0) {
            mf.tre[i,"distance"] <- distTemp
            
        }
    }
    
    # print(paste0("done",i," out of ",nrow(mf.tre)))
}

### Get proportion and count of inhibitory otus from inhibitory otu metadata
inhib.tb <- inhib %>%
    as_tibble() %>%
    rename(Name=V1, inhib=V2, num=V3, seq=V4)

otu.tb.inhib <- otu.tb %>%
    mutate(inhib = (inhib.tb[match(otu.tb$OTUID, inhib.tb$seq),"inhib"])$inhib)

# for tre
otu.tb.inhib.tre <- otu.tre %>%
    mutate(inhib = (inhib.tb[match(otu.tre$OTUID, inhib.tb$seq),"inhib"])$inhib)

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

## now do for tre
# total counts
sampleCounts <- otu.tb.inhib.tre %>%
    dplyr::select(-c("inhib","OTUID")) %>%
    colSums()
# counts of total inhibitory sequences
inhibCounts <- otu.tb.inhib.tre %>%
    filter(inhib=="inhibitory") %>%
    dplyr::select(-c("inhib","OTUID")) %>%
    colSums()
# richness of inhibitory sequences
inhibRich <- otu.tb.inhib.tre %>%
    filter(inhib=="inhibitory") %>%
    dplyr::select(-c("inhib","OTUID")) %>%
    replace(.>0, 1) %>%
    colSums()
## Now, add this information into the mf
mf.tre.inhib <- mf.tre %>%
    mutate(n=sampleCounts[match(SampleID, names(sampleCounts))]
           ,inhibCounts =inhibCounts[match(SampleID, names(inhibCounts))]
           , inhibRich=inhibRich[match(SampleID, names(inhibRich))]) %>%
    mutate(percInhib = inhibCounts/n)

### Let's delete all intermediate files to save environmental space
# OTUt ables are an insane amount of space...
remove(otu)
remove(otu.tb.inhib)
remove(otu.tb)
remove(mf.tb)
remove(mf.tre)
remove(otu.tre)
remove(otu.tre.inhib)

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
grid.arrange(gg_logRich, gg_shannon, gg_distance, gg_inhibRich, gg_percInhib, nrow=5)

quartz()
grid.arrange(gg_logRich,gg_inhibRich, gg_percInhib, nrow=3)
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
mf.indiv5 <- mf.tb.inhib %>%
    dplyr::select(toadID, species, BD_infected, toadID, shannon, logRich, distance, inhibRich, percInhib, time) %>%
    arrange(toadID) %>%
    filter(time == 5) %>%
    full_join(mf.indiv)
p2 <- mf.indiv5 %>%
    filter(BD_infected=="y") %>%
    ggplot(aes(y=log(Pos_maxBD+0.1)))

#----
# Between species
quartz()
p2 + geom_point(aes(x=species, color=species, group=toadID), position=position_jitter(width=0.1, height=0), cex=5)

# So it looks like Bubo is particularily susceptible to BD infection once it takes hold...
#----

# as a factor of community diversity 
# p2 + geom_point(aes(x=Pre_div, color=species), cex=5) +
#     geom_smooth(aes(x=Pre_div), method="lm")
quartz()
p2 + geom_point(aes(x=Pre_div, color=species), cex=5) +
    labs(x="Shannon Diversity", y="Log max BD +0.1")
# as a factor of richness
# p2 + geom_point(aes(x=Pre_rich, color=species), cex=5) + 
#     geom_smooth(aes(x=Pre_rich), method="lm")
quartz()
p2 + geom_point(aes(x=Pre_rich, color=species), cex=5) + 
    labs(x="OTU richness", "Log max BD +0.1")
#----
# as a factor of inhibitory richness or percent of microbiome?
# p2 + geom_point(aes(x=Pre_inhibRich, color=species),position=position_jitter(width=0.1, height=0.1), cex=5) +
#     geom_smooth(aes(x=Pre_inhibRich), method="lm")
quartz()
p2 + geom_point(aes(x=Pre_inhibRich, color=species),position=position_jitter(width=0.1, height=0.1), cex=5) +
    labs(x="Richness of BD-inhibitory bacteria", y="Log max BD +0.1")
quartz()
p2 + geom_point(aes(x=Pre_percInhib, color=species),position=position_jitter(width=0.1, height=0.1), cex=5) +
    labs(x="Percent BD-inhibitory bacteria (of community)", y="Log max BD +0.1")
# p2 + geom_point(aes(x=Pre_percInhib, color=species),position=position_jitter(width=0.1, height=0.1), cex=5) +
#     geom_smooth(aes(x=Pre_percInhib), method="lm")

# There seems to be a correlation between the richness and proportion of microbiome that is inhibitory, but not actually microiome richness

# What about distance travelled in "multidimensional space" over time?
p2 + geom_point(aes(x=Pre_dist, color=species), cex=5) +
    geom_smooth(aes(x=Pre_dist), method="lm")

##### STEP TWO (B): Look at all individuals who were exposed, infected is categorical and not continuous #####

# jitter <- position_jitter(width=0, height=1)
p2b <- mf.indiv5 %>%
    filter(BD_infected == "y") %>%
    ggplot(aes(y=Pos_PABD))
mf.indiv5
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
#+ warnings=FALSE, results="hide", message=FALSE, cache=TRUE
glm_toad_inhibRich <- stan_glm(Pos_PABD ~ -1 + species + Pre_inhibRich, data=mf.exposed.t5, family=binomial(link="logit"))
#+ warnings=FALSE, results="hide", message=FALSE, cache=TRUE
glm_toad_percInhib <- stan_glm(Pos_PABD ~ -1 + species + Pre_percInhib, data=mf.exposed.t5, family=binomial(link="logit"))
#+ warnings=FALSE, results="hide", message=FALSE, cache=TRUE
glm_toad_inhib <- stan_glm(Pos_PABD ~ -1 + species + Pre_inhibRich:Pre_percInhib, data=mf.exposed.t5, family=binomial(link="logit"))

samps1 <- extract(glm_toad_inhibRich$stanfit)
samps2 <- extract(glm_toad_percInhib$stanfit)
samps3 <- extract(glm_toad_inhib$stanfit)
df1 <- data.frame(samp = as.vector(samps1$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("XInhib",4000))
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
rbind(df1,df2,df3) %>%
    filter(parameter=="XInhib") %>%
    ggplot() +
    geom_histogram(aes(x=samp), bins=20) +
    geom_vline(mapping=aes(xintercept=0), col="red")+
    geom_text(data=overZero, aes(x=x,y=y,label=label), col="blue") +
    facet_wrap(~model, nrow=1, scales="free")

#' It looks like the multiplicative model is actually the best model? The samples suggest that there are the fewest
#' parameter estimates that are over 0 than the other two models.
#' Let's look at the model fit for parameter number three:

# Let's look at the model fit
samps_jointinhibmodel <- samps3$beta %>%
    as_tibble() %>%
    rename("Bubo"=V1,"Buma"=V2,"Osse"=V3,"Raca"=V4,"Rapi"=V5,"Inhib"=V6)

# get predictor x vector
inhib_jointeffect <- mf.exposed.t5 %>%
    transmute(joint_eff = Pre_inhibRich*Pre_percInhib) %>%
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
    geom_point(mapping=aes(x=Pre_inhibRich*Pre_percInhib, y=Pos_PABD), position=position_jitter(width=0, height=0.05)) +
    geom_line(data=summary_all, aes(x=x, y=mean, color=species)) +
    geom_ribbon(data=summary_all, aes(x=x, ymin=lwr, ymax=upr, bg=species), alpha=0.2) +
    facet_wrap(~species, nrow=1) +
    labs(x="Interaction of richness and percent composition of inhibitory bacteria", y="Presence/Absence BD")

#' It looks like there are certain species who are more likely to get infected than others, and there is a lower chance of infection
#' if your microbiome is RICHER (but not necessarily more proportionally inhibitory bacteria).
#' The effect of a slightly richer microbiome is very subtle, but it is fairly precise\


#' ##### Does richness or proportion of BD-inhibitory bacteria affect INTENSITY of infection? #####

mf.infected.t5 <- mf.indiv5 %>%
    filter(BD_infected == "y", Pos_PABD == 1) %>%
    mutate(logmaxBD = log(Pos_maxBD +0.1))

#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glm_toad_inhibRich_intensity <- stan_glm(logmaxBD ~ -1 + species + Pre_inhibRich, data=mf.infected.t5)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glm_toad_percInhib_intensity <- stan_glm(logmaxBD ~ -1 + species + Pre_percInhib, data=mf.infected.t5)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glm_toad_inhib_intensity <- stan_glm(logmaxBD ~ -1 + species + Pre_percInhib:Pos_inhibRich, data=mf.infected.t5)

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
rbind(df1,df2,df3) %>%
    filter(parameter=="XInhib") %>%
    as.data.frame() %>%
    ggplot(aes(x=samp)) +
    geom_histogram(aes(y=..count..),bins=50) +
    facet_grid(~model, scale="free")+
    geom_vline(aes(xintercept=0), col="red") +
    geom_text(data=overZero, aes(x=x,y=y,label=label), col="blue")


# Let's look at the model fit
samps_inhibRichmodel <- samps1$beta %>%
    as_tibble() %>%
    rename("Bubo"=V1,"Buma"=V2,"Osse"=V3,"Raca"=V4,"Inhib"=V5)

# get predictor x vector
jointeffect_x <- seq(max(c(0,min(mf.infected.t5$Pre_inhibRich))), max(mf.infected.t5$Pre_inhibRich), length.out = 100)

# run for all species
#Bubo
pred_Bubo <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_inhibRichmodel$Bubo + jointeffect_x[i]*samps_inhibRichmodel$Inhib
    pred_Bubo[,i] <- temp
}
summary_Bubo <- data.frame(mean=colMeans(pred_Bubo), upr=apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), lwr=-apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), x=jointeffect_x, species=rep("Bubo",100))
#Buma
pred_Buma <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_inhibRichmodel$Buma + jointeffect_x[i]*samps_inhibRichmodel$Inhib
    pred_Buma[,i] <- temp
}
summary_Buma <- data.frame(mean=colMeans(pred_Buma), upr=apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), lwr=-apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), x=jointeffect_x, species=rep("Buma",100))
#Osse
pred_Osse <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_inhibRichmodel$Osse + jointeffect_x[i]*samps_inhibRichmodel$Inhib
    pred_Osse[,i] <- temp
}
summary_Osse <- data.frame(mean=colMeans(pred_Osse), upr=apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), lwr=-apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), x=jointeffect_x, species=rep("Osse",100))

#Raca
pred_Raca <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_inhibRichmodel$Raca + jointeffect_x[i]*samps_inhibRichmodel$Inhib
    pred_Raca[,i] <- temp
}
summary_Raca <- data.frame(mean=colMeans(pred_Raca), upr=apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), lwr=-apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), x=jointeffect_x, species=rep("Raca",100))

# make into one
summary_all <- rbind(summary_Bubo, summary_Buma, summary_Osse, summary_Raca)
#+ fig.width=12, fig.height=5
quartz()
ggplot(data=mf.infected.t5) +
    geom_point(mapping=aes(x=Pre_inhibRich, y=logmaxBD), position=position_jitter(width=0, height=0.05)) +
    geom_line(data=summary_all, aes(x=x, y=mean, color=species)) +
    geom_ribbon(data=summary_all, aes(x=x, ymin=lwr, ymax=upr, bg=species), alpha=0.2) +
    facet_wrap(~species, nrow=1) +
    labs(x="Richness of inhibitory bacteria", y="Log max BD")


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
lmer_controls_inhibRich_logRich <- stan_lmer(inhibRich ~ -1 + species + logRich + time + (1|toadID), data=mf.notinfect)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
lmer_controls_inhibRich_div <- stan_lmer(inhibRich ~ -1 + species + shannon + time + (1|toadID), data=mf.notinfect)


samps1 <- extract(lmer_controls_inhibRich_logRich$stanfit)
samps2 <- extract(lmer_controls_inhibRich_div$stanfit)
df1 <- data.frame(samp = as.vector(samps1$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("Xdiv",4000), rep("xTime",4000))
                  , model=rep("logRich",4000*7)
) 
byToad1 <- ranef(lmer_controls_inhibRich_logRich)$toadID

df2 <- data.frame(samp = as.vector(samps2$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("Xdiv",4000), rep("xTime",4000))
                  , model=rep("div",4000*7)
) 
byToad2 <- ranef(lmer_controls_inhibRich_div)$toadID

#+ fig.height=8, fig.width=10
rbind(df1,df2)%>%
    filter(parameter=="Xdiv") %>%
    as.data.frame() %>%
    ggplot(aes(x=samp)) +
    geom_histogram(aes(y=..count..),bins=50) +
    geom_vline(aes(xintercept=0), col="red")+
    facet_grid(~model, scale="free")

# Let's look at the model fit of the best model
samps_shannonmodel <- samps2$beta %>%
    as_tibble() %>%
    rename("Bubo"=V1,"Buma"=V2,"Osse"=V3,"Raca"=V4,"Rapi"=V5,"shannon"=V6, "Time"=V7)

# get predictor x vector
inhib_effect <- mf.notinfect %>%
    dplyr::select(shannon) %>%
    pull()
jointeffect_x <- seq(max(0,min(inhib_effect,na.rm=TRUE)), max(inhib_effect, na.rm=TRUE), length.out = 100)

# run for all species
#Bubo
pred_Bubo <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_shannonmodel$Bubo + jointeffect_x[i]*samps_shannonmodel$shannon
    pred_Bubo[,i] <- (temp)
}
summary_Bubo <- data.frame(mean=colMeans(pred_Bubo), upr=apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), lwr=-apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), x=jointeffect_x, species=rep("Bubo",100))
#Buma
pred_Buma <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_shannonmodel$Buma + jointeffect_x[i]*samps_shannonmodel$shannon
    pred_Buma[,i] <-(temp)
}
summary_Buma <- data.frame(mean=colMeans(pred_Buma), upr=apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), lwr=-apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), x=jointeffect_x, species=rep("Buma",100))
#Osse
pred_Osse <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_shannonmodel$Osse + jointeffect_x[i]*samps_shannonmodel$shannon
    pred_Osse[,i] <- (temp)
}
summary_Osse <- data.frame(mean=colMeans(pred_Osse), upr=apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), lwr=-apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), x=jointeffect_x, species=rep("Osse",100))
#Raca
pred_Raca <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_shannonmodel$Raca + jointeffect_x[i]*samps_shannonmodel$shannon
    pred_Raca[,i] <- (temp)
}
summary_Raca <- data.frame(mean=colMeans(pred_Raca), upr=apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), lwr=-apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), x=jointeffect_x, species=rep("Raca",100))
#Rapi
pred_Rapi <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x)) {
    temp <- samps_shannonmodel$Rapi + jointeffect_x[i]*samps_shannonmodel$shannon
    pred_Rapi[,i] <- (temp)
}
summary_Rapi <- data.frame(mean=colMeans(pred_Rapi), upr=apply(pred_Rapi,MARGIN=2,sd) + colMeans(pred_Rapi), lwr=-apply(pred_Rapi,MARGIN=2,sd) + colMeans(pred_Rapi), x=jointeffect_x, species=rep("Rapi",100))

# make into one
summary_all <- rbind(summary_Bubo, summary_Buma, summary_Osse, summary_Raca, summary_Rapi)
# Adjust the shannon values to account for time and individual
mean_Timeeff <- as.vector(colMeans(samps_shannonmodel[,"Time"]))
mf.notinfect <- mf.notinfect %>%
    mutate(adjinhibRich = inhibRich-time*mean_Timeeff-byToad1[toadID,1])
#+ fig.width=12, fig.height=5
# filter out everything except first time point
ggplot(data=mf.notinfect) +
    geom_point(mapping=aes(x=shannon, y=adjinhibRich)) +
    geom_line(data=summary_all, aes(x=x, y=(mean), color=species)) +
    geom_ribbon(data=summary_all, aes(x=x, ymin=(lwr), ymax=(upr), bg=species), alpha=0.2) +
    facet_wrap(~species, nrow=1)



#' Overall richness doesn't influence richness of inibitory microbes as much as "diversity"-- 
#' so balanced microbiota that are more even tend to have more inhibitory bacteria. (In addition, time seems to affect diversity slightly as well). In general, however,
#' there is actually not that strong of an effect of either shannon or OTU richness on inhibitory bacteria richness... or at least, it doesn't appear that way.
#' 
#' 
#' \
#' \
#' Let's check to see if log percent inhibitory is normal-looking
ggplot(mf.notinfect, aes(x=log(percInhib))) +
    geom_histogram(bins=20) +
    facet_grid(~species)
# Yup; looks normal now.

#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
lmer_controls_percInhib_logRich <- stan_glmer(percInhib ~ -1 +species + logRich + time + (1|toadID), data=mf.notinfect, family=gaussian(link="log"))
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
lmer_controls_percInhib_div <- stan_glmer(percInhib ~ -1 +species + shannon + time + (1|toadID), data=mf.notinfect,, family=gaussian(link="log"))

samps1 <- extract(lmer_controls_percInhib_logRich$stanfit)
samps2 <- extract(lmer_controls_percInhib_div$stanfit)
df1 <- data.frame(samp = as.vector(samps1$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("Xdiv",4000), rep("xTime",4000))
                  , model=rep("logRich",4000*7)
) 
bytoad1 <- ranef(lmer_controls_percInhib_logRich)$toadID
df2 <- data.frame(samp = as.vector(samps2$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("Xdiv",4000), rep("xTime",4000))
                  , model=rep("shannon",4000*7)
) 
bytoad2 <- ranef(lmer_controls_percInhib_div)$toadID

#+ fig.height=8, fig.width=10
rbind(df1,df2)%>%
    filter(parameter=="Xdiv") %>%
    as.data.frame() %>%
    ggplot(aes(x=samp)) +
    geom_histogram(aes(y=..count..),bins=50) +
    geom_vline(aes(xintercept=0), col="red")+
    facet_grid(~model, scale="free")

# Let's look at the model fit of the best model (they're both good, so let's look at both)
samps_shannonmodel <- samps2$beta %>%
    as_tibble() %>%
    rename("Bubo"=V1,"Buma"=V2,"Osse"=V3,"Raca"=V4,"Rapi"=V5,"shannon"=V6, "Time"=V7)
samps_obsotunmodel <- samps1$beta %>%
    as_tibble() %>%
    rename("Bubo"=V1,"Buma"=V2,"Osse"=V3,"Raca"=V4,"Rapi"=V5,"logRich"=V6, "Time"=V7)

# get predictor x vector
inhib_effect_shannon <- mf.notinfect %>%
    dplyr::select(shannon) %>%
    pull()
jointeffect_x_shannon <- seq(max(c(0,min(inhib_effect_shannon, na.rm=TRUE))), max(inhib_effect_shannon, na.rm=TRUE), length.out = 100)
inhib_effect_obsotus <- mf.notinfect %>%
    dplyr::select(logRich) %>%
    pull()
jointeffect_x_obsotus <- seq(max(c(0,min(inhib_effect_obsotus, na.rm=TRUE))), max(inhib_effect_obsotus, na.rm=TRUE), length.out = 100)

# run for all species (shannon)
#Bubo
pred_Bubo <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x_shannon)) {
    temp <- samps_shannonmodel$Bubo + jointeffect_x_shannon[i]*samps_shannonmodel$shannon
    pred_Bubo[,i] <- exp(temp)
}
summary_Bubo <- data.frame(mean=colMeans(pred_Bubo), upr=apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), lwr=-apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), x=jointeffect_x_shannon, species=rep("Bubo",100))
#Buma
pred_Buma <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x_shannon)) {
    temp <- samps_shannonmodel$Buma + jointeffect_x_shannon[i]*samps_shannonmodel$shannon
    pred_Buma[,i] <-exp(temp)
}
summary_Buma <- data.frame(mean=colMeans(pred_Buma), upr=apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), lwr=-apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), x=jointeffect_x_shannon, species=rep("Buma",100))
#Osse
pred_Osse <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x_shannon)) {
    temp <- samps_shannonmodel$Osse + jointeffect_x_shannon[i]*samps_shannonmodel$shannon
    pred_Osse[,i] <- exp(temp)
}
summary_Osse <- data.frame(mean=colMeans(pred_Osse), upr=apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), lwr=-apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), x=jointeffect_x_shannon, species=rep("Osse",100))
#Raca
pred_Raca <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x_shannon)) {
    temp <- samps_shannonmodel$Raca + jointeffect_x_shannon[i]*samps_shannonmodel$shannon
    pred_Raca[,i] <- exp(temp)
}
summary_Raca <- data.frame(mean=colMeans(pred_Raca), upr=apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), lwr=-apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), x=jointeffect_x_shannon, species=rep("Raca",100))
#Rapi
pred_Rapi <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x_shannon)) {
    temp <- samps_shannonmodel$Rapi + jointeffect_x_shannon[i]*samps_shannonmodel$shannon
    pred_Rapi[,i] <- exp(temp)
}
summary_Rapi <- data.frame(mean=colMeans(pred_Rapi), upr=apply(pred_Rapi,MARGIN=2,sd) + colMeans(pred_Rapi), lwr=-apply(pred_Rapi,MARGIN=2,sd) + colMeans(pred_Rapi), x=jointeffect_x_shannon, species=rep("Rapi",100))
# make into one
summary_all_shannon <- rbind(summary_Bubo, summary_Buma, summary_Osse, summary_Raca, summary_Rapi)

# run for all species (obs otus)
#Bubo
pred_Bubo <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x_obsotus)) {
    temp <- samps_obsotunmodel$Bubo + jointeffect_x_obsotus[i]*samps_obsotunmodel$logRich
    pred_Bubo[,i] <- exp(temp)
}
summary_Bubo <- data.frame(mean=colMeans(pred_Bubo), upr=apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), lwr=-apply(pred_Bubo,MARGIN=2,sd) + colMeans(pred_Bubo), x=jointeffect_x_obsotus, species=rep("Bubo",100))
#Buma
pred_Buma <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x_obsotus)) {
    temp <- samps_obsotunmodel$Buma + jointeffect_x_obsotus[i]*samps_obsotunmodel$logRich
    pred_Buma[,i] <-exp(temp)
}
summary_Buma <- data.frame(mean=colMeans(pred_Buma), upr=apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), lwr=-apply(pred_Buma,MARGIN=2,sd) + colMeans(pred_Buma), x=jointeffect_x_obsotus, species=rep("Buma",100))
#Osse
pred_Osse <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x_obsotus)) {
    temp <- samps_obsotunmodel$Osse + jointeffect_x_obsotus[i]*samps_obsotunmodel$logRich
    pred_Osse[,i] <- exp(temp)
}
summary_Osse <- data.frame(mean=colMeans(pred_Osse), upr=apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), lwr=-apply(pred_Osse,MARGIN=2,sd) + colMeans(pred_Osse), x=jointeffect_x_obsotus, species=rep("Osse",100))
#Raca
pred_Raca <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x_obsotus)) {
    temp <- samps_obsotunmodel$Raca + jointeffect_x_obsotus[i]*samps_obsotunmodel$logRich
    pred_Raca[,i] <- exp(temp)
}
summary_Raca <- data.frame(mean=colMeans(pred_Raca), upr=apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), lwr=-apply(pred_Raca,MARGIN=2,sd) + colMeans(pred_Raca), x=jointeffect_x_obsotus, species=rep("Raca",100))
#Rapi
pred_Rapi <- matrix(ncol=100, nrow=4000)
for ( i in 1:length(jointeffect_x_obsotus)) {
    temp <- samps_obsotunmodel$Rapi + jointeffect_x_obsotus[i]*samps_obsotunmodel$logRich
    pred_Rapi[,i] <- exp(temp)
}
summary_Rapi <- data.frame(mean=colMeans(pred_Rapi), upr=apply(pred_Rapi,MARGIN=2,sd) + colMeans(pred_Rapi), lwr=-apply(pred_Rapi,MARGIN=2,sd) + colMeans(pred_Rapi), x=jointeffect_x_obsotus, species=rep("Rapi",100))
# make into one
summary_all_obsotus <- rbind(summary_Bubo, summary_Buma, summary_Osse, summary_Raca, summary_Rapi)


# Adjust the shannon values to account for time
mean_Timeeff <- as.vector(colMeans(samps_shannonmodel[,"Time"]))
mf.notinfect <- mf.notinfect %>%
    mutate(adjpercInhib2 = log(percInhib)-time*mean_Timeeff-bytoad2[toadID,])
mean_Timeeff <- as.vector(colMeans(samps_obsotunmodel[,"Time"]))
mf.notinfect <- mf.notinfect %>%
    mutate(adjpercInhib1 = log(percInhib)-time*mean_Timeeff-bytoad1[toadID,])
#+ fig.width=12, fig.height=5, warning=FALSE, message=FALSE
# filter out everything except first time point
ggplot(data=mf.notinfect) +
    geom_point(mapping=aes(x=shannon, y=exp(adjpercInhib2))) +
    geom_line(data=summary_all_shannon, aes(x=x, y=(mean), color=species)) +
    geom_ribbon(data=summary_all_shannon, aes(x=x, ymin=(lwr), ymax=(upr), bg=species), alpha=0.2) +
    facet_wrap(~species, nrow=1)
ggplot(data=mf.notinfect) +
    geom_point(mapping=aes(x=logRich, y=exp(adjpercInhib1))) +
    geom_line(data=summary_all_obsotus, aes(x=x, y=(mean), color=species)) +
    geom_ribbon(data=summary_all_obsotus, aes(x=x, ymin=(lwr), ymax=(upr), bg=species), alpha=0.2) +
    facet_wrap(~species, nrow=1)

#' Huh. There is actually a negative relationship between overall richness and diversity and proportion of the microbiome that is inhibitory
#' That means there IS a decoupling of microbiome richness vs richness of inhibitory bacteria
#' Which could explain why there is so much variation in studies of microbiome effects on BD inhibition?
#' Or perhaps what is actually means is that when microbiomes are diverse, we should have proportionally characterized less of them
#' It might SEEM like there are less inhibitory microbes, proportionally, but maybe the  uncharacterized microbes in high diversity samples
#' are driving this pattern...

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

#filter to just post-exposure
mf.pos <- mf.exposed %>%
    filter(prepost == "Pos") 
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glmer_PABD_affects_betadiv <- stan_lmer(distance ~ -1 + species + PABD + (1|toadID), data=mf.exposed)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glmer_PABD_affects_shannon <- stan_lmer(shannon ~ -1 + species + PABD + (1|toadID), data=mf.exposed)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glmer_PABD_affects_rich <- stan_lmer(logRich ~ -1 + species + PABD + (1|toadID), data=mf.exposed)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glmer_PABD_affects_inhibRich <- stan_lmer(inhibRich ~ -1 + species + PABD + (1|toadID), data=mf.exposed)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glmer_PABD_affects_percInhib <- stan_glmer(percInhib ~ -1 + species + PABD + (1|toadID), data=mf.exposed, family=gaussian(link="log"))

samps1 <- extract(glmer_PABD_affects_betadiv$stanfit)
samps2 <- extract(glmer_PABD_affects_shannon$stanfit)
samps3 <- extract(glmer_PABD_affects_rich$stanfit)
samps4 <- extract(glmer_PABD_affects_inhibRich$stanfit)
samps5 <- extract(glmer_PABD_affects_percInhib$stanfit)

df1 <- data.frame(samp = as.vector(samps1$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("XPABD",4000))
                  , model=rep("betadiv",4000*6)
) 
df2 <- data.frame(samp = as.vector(samps2$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("XPABD",4000))
                  , model=rep("shannon",4000*6)
) 
df3 <- data.frame(samp = as.vector(samps3$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("XPABD",4000))
                  , model=rep("rich",4000*6)
)
df4 <- data.frame(samp = as.vector(samps4$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("XPABD",4000))
                  , model=rep("inhibRich",4000*6)
)
df5 <- data.frame(samp = as.vector(samps5$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("Rapi",4000),rep("XPABD",4000))
                  , model=rep("percInhib",4000*6)
)
df1.filt <- df1 %>%
    filter(parameter == "XPABD")
df2.filt <- df2 %>%
    filter(parameter == "XPABD")
df3.filt <- df3 %>%
    filter(parameter == "XPABD")
df4.filt <- df4 %>%
    filter(parameter == "XPABD")
df5.filt <- df5 %>%
    filter(parameter == "XPABD")
#+ fig.height=10, fig.width=8
as.data.frame(rbind(df1.filt,df2.filt,df3.filt,df4.filt,df5.filt)) %>%
    ggplot(aes(x=samp)) +
    geom_histogram(aes(y=..density..),bins=100) +
    facet_grid(model ~ parameter, scales="free_y") +
    geom_vline(aes(xintercept=0), col="red")


#' Presence and absence of BD does NOT affect percent of inhibitory bacteria within community. 
#' BD has a highly variable effect on inhibitory bacterial richness-- and it doesn't seem to be "significantly" affecting it one way or another
#' BD seems to be associated with a slight increase in OTU richness (perhaps not significant) but a slight decrease in diversity. This probably means the community is becoming less even, but there are more OTUs.
#' Finally, BD seems to actually make distance travelled smaller-- the community is more stable when it is infected with BD than when not
#' But I wonder if this is because:
#' (a) things that are sick all look the same
#' (b) there are particular taxa that thrive in the presence of BD and these make it look all similar


# Now let's see with BD intensity.
mf.pos.infect <- mf.exposed %>%
    filter(prepost == "Pos", PABD > 0) 
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glmer_aveBD_affects_betadiv <- stan_lmer(distance ~ -1 + species + log(aveBD) + (1|toadID), data=mf.pos.infect)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glmer_aveBD_affects_shannon <- stan_lmer(shannon ~ -1 + species + log(aveBD) + (1|toadID), data=mf.pos.infect)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glmer_aveBD_affects_rich <- stan_lmer(logRich ~ -1 + species + log(aveBD) + (1|toadID), data=mf.pos.infect)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glmer_aveBD_affects_inhibRich <- stan_lmer(inhibRich ~ -1 + species + log(aveBD) + (1|toadID), data=mf.pos.infect)
#+ message=FALSE, results="hide", warnings=FALSE, cache=TRUE
glmer_aveBD_affects_percInhib <- stan_glmer(percInhib ~ -1 + species + log(aveBD) + (1|toadID), data=mf.pos.infect, family=gaussian(link="log"))

samps1 <- extract(glmer_aveBD_affects_betadiv$stanfit)
samps2 <- extract(glmer_aveBD_affects_shannon$stanfit)
samps3 <- extract(glmer_aveBD_affects_rich$stanfit)
samps4 <- extract(glmer_aveBD_affects_inhibRich$stanfit)
samps5 <- extract(glmer_aveBD_affects_percInhib$stanfit)

df1 <- data.frame(samp = as.vector(samps1$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("XaveBD",4000))
                  , model=rep("betadiv",4000*4)
) 
df2 <- data.frame(samp = as.vector(samps2$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("XaveBD",4000))
                  , model=rep("shannon",4000*5)
) 
df3 <- data.frame(samp = as.vector(samps3$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("XaveBD",4000))
                  , model=rep("rich",4000*5)
)
df4 <- data.frame(samp = as.vector(samps4$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("XaveBD",4000))
                  , model=rep("inhibRich",4000*5)
)
df5 <- data.frame(samp = as.vector(samps5$beta)
                  , parameter=c(rep("Bubo",4000),rep("Buma",4000),rep("Osse",4000),rep("Raca",4000),rep("XaveBD",4000))
                  , model=rep("percInhib",4000*5)
)
df1.filt <- df1 %>%
    filter(parameter == "XaveBD")
df2.filt <- df2 %>%
    filter(parameter == "XaveBD")
df3.filt <- df3 %>%
    filter(parameter == "XaveBD")
df4.filt <- df4 %>%
    filter(parameter == "XaveBD")
df5.filt <- df5 %>%
    filter(parameter == "XaveBD")
#+ fig.height=10, fig.width=8
as.data.frame(rbind(df1.filt,df2.filt,df3.filt,df4.filt,df5.filt)) %>%
    ggplot(aes(x=samp)) +
    geom_histogram(aes(y=..density..),bins=100) +
    facet_grid(model ~ parameter, scales="free_y") +
    geom_vline(aes(xintercept=0), col="red")

#' Whoa! The infection intensity of BD is associated with a decrease in percent inhibitory bacteria within the microbiome! Unsure if this is cause or effect...
#' \
#' Now let's see if the same trend holds for second round of treatment
#' 

View(mf.tb.inhib)

compare_richness_2 <- mf.tre.inhib %>%
    mutate(experiment = "SecondTreatment") %>%
    dplyr::select(toadID, percInhib, inhibRich, observed_otus, shannon, distance, prepost, BD_infected, aveBD, PABD, experiment, secondexp_con)
compare_richness_1 <- mf.tb.inhib %>%
    filter(species == "Bubo") %>%
    mutate(experiment = "FirstTreatment") %>%
    dplyr::select(toadID, percInhib, inhibRich, observed_otus, shannon, distance, prepost, BD_infected, aveBD, PABD, experiment,secondexp_con)
bubo_all <- rbind(compare_richness_1, compare_richness_2)

bubo_all %>%
ggplot() +
    geom_violin(aes(x=experiment, y=percInhib, col=toadID)
                # , position=position_jitter(width=0.1, height=0)
                ) +
    facet_wrap(~toadID, ncol=2)

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
quartz()
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
quartz()
mf.tb.inhib %>%
    mutate(prepost_full = ifelse(prepost=="Pos","Post-exposure","Pre-exposure"), BD_inf_full = ifelse(BD_infected == "y","Exposed","Not Exposed")) %>%
    mutate(prepost_full = factor(prepost_full, levels=c("Pre-exposure","Post-exposure"))) %>%
    ggplot() +
    geom_point(aes(x=inhibRich, y=percInhib, col=species)) +
    facet_grid(BD_inf_full ~ prepost_full) +
    labs(x="Richness of BD-inhibitory bacteria", y="Percent of OTUs that are inhibitory toward BD") + 
    geom_smooth(aes(x=inhibRich, y=percInhib, col=species), method = "lm")


