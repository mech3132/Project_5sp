#' ---
#' title: Exploratory Data Analysis for 5Sp Dataset
#' author: Melissa Chen
#' output: github_document
#' 
#' ---
#' 

#' ## First, load required packages\
#### Load packages ####
library(tidyverse) # for data manipulation
library(MASS) # For fitting distributions and NMDS
library(rstanarm) # For calculated expected alpha and beta diversity
library(gridExtra) # For arranging ggplots
library(mgcv) # For beta distribution (beta diversity)
library(vegan) # for permanova
library(car) # for type III Anova
library(RColorBrewer) # colors for barplots
library(grid) # for text grobs in gridExtra
#' Define pathways for input files

#### Pathways ####
# OTU table in text format; rarefied to 'rared' OTUs
oturarePWD="../MF_and_OTU_edited/otu_table_r5000.txt"
# rarefaction depth
rared <- 5000
# OTU table in text format; not rarefied
otuPWD = "../MF_and_OTU_edited/otu_table_text.txt"
# FULL Mapping file; with alpha diversity added-- must have all columns described in aMetrics
mfPWD = "../MF_and_OTU_edited/MF_withalpha.txt"
# aMetrics
aMetrics <- c("observed_otus","shannon")
# Get header names for all alpha metrics
aHeader <- paste0(aMetrics,"_even_",rared,"_alpha")

# Inhibitory metadata; columns must be "host","inhib status","isolate number","OTU ID", "fulldescrp","OTU taxonomy"
inhibPWD <- "../ANTIFUNGAL/inhibitory_metadata_MANUAL.txt"

# folder with distance matrices
dmPWD <- "../beta_div/"
# distance metrics
dmMetrics <- c("bray_curtis")

minOTUSample = 5 # min number of reads per sample for the OTU to be non-zero
minOTUTab = 100 # min number of OTUs found in table for OTU to be kept in table
otuCutoff = 5000 # cutoff to discard sample; 



#### FILTERING BD NOTES:
# First, checks all numbers to see if less than indivTHRESH. If less than indivTHRESH, it is changed to '0'
# Then, it see is if at least 2 are NOT zero and the third is more than 50. 
# If the third is less than 50 AND the other two measurements are zero, they are all changed to zeros./


#### LOAD DATA ####

# Load mappingfile 
mf <- read.delim(paste0(mfPWD), header=TRUE, as.is=TRUE)
# Load full OTU table
otu <- read.delim(paste0(otuPWD), header=TRUE, as.is=TRUE, skip=1)
# Load rarefied OTU table
otu_rare <- read.delim(paste0(oturarePWD), header=TRUE, as.is=TRUE, skip=1)
# Loop through DMs
dm_all <- list()
for ( d in dmMetrics) {
    dm_all[[d]] <-  read.delim(paste0(dmPWD,dmMetrics,"_dm.txt"), header=TRUE)
}
# inhib
inhib <- read.delim(paste0(inhibPWD), header=FALSE, as.is=TRUE)

#' ## Determining BD infection loads ## \
#' \
#' One of the problems with the BD qPCR results is that we get very irregular results. Thus, each individual measurement
#' is unreliable. Here, I use a parameterized model to predict the "true" Bd load given the measurements taken.
#'  I would expect BD load to be modelled by an approximately poisson process; here, we check if this is true.

# Filter out all information except BD-positive scenarios
BD <- mf %>%
    as_tibble() %>%
    dplyr::select(X.SampleID, Bd_Run_1,Bd_Run_2, Bd_Average_Run_3) %>%
    filter(!is.na(Bd_Average_Run_3), !is.na(Bd_Run_1), !is.na(Bd_Run_2))
samples_bd <- BD$X.SampleID
BD <- BD %>%
    dplyr::select(-X.SampleID)
BD12 <- c(BD$Bd_Run_1, BD$Bd_Run_2, BD$Bd_Average_Run_3)
# Get rid of zeros due to overinflation in poisson model
BD12 <- BD12[BD12!=0]
# Fit a poisson model to log BD
poisfit <- MASS::fitdistr(round(log(BD12)), densfun = "Poisson")
# Set new range of X's to test
xfit <- seq(0,(max(log(BD12))))
# Predict y
pred.y <- data.frame(y.pred=dpois(x=xfit, poisfit$estimate), xfit=xfit)
# Plot histogram with poisson distribution fit
BD12 %>%
    as_tibble() %>%
    rename(BDload=value)%>%
    ggplot(aes(x=log(BDload))) +
    geom_histogram(aes(y=..density..), bins=8) +
    geom_line(data=pred.y,aes(x=xfit, y=y.pred), col="red")

#' What we see above is that BD load (when not zero) is, in fact, modelled well by a poisson process. This means we could fit a poisson
#' model to the aPCR results to estimate the "true" infection load (lambda) for each toad.

# Model true infection load for each toad
BD_lambda_est <- t(apply(round(log(BD+1)), MARGIN=1, FUN=function(x) {
    temp <- fitdistr(x, densfun="Poisson")
    return(c(temp$estimate,temp$sd))
})) %>%
    as_tibble() %>%
    rename(sd=lambda1) %>% 
    mutate(pval=dnorm(0,mean=lambda,sd=sd)) %>% # To estimate if lambda is significantly different from zero, we see if pval for parameter estimate different from zero
    mutate(sig=pval<=0.10)
final_BD <- cbind(BD,BD_lambda_est)

# To compare, let's see what a different filtering method yields:
# If 2 BD samples are 0 and the third is less than 50, then set to zero.
# Also, if anything is less than 5, make it zero anyway.

# BD filtering notes and thresholds:
indivTHRESH = 5 # BD individual threshold
thirdTHRESH = 50 # BD 3rd sample threshold
BD_alt <- BD
for ( r in 1:nrow(BD_alt) ) {
    for ( c in 1:ncol(BD_alt) ) {
        if ( (BD_alt[r,c] < indivTHRESH) | is.na(BD_alt[r,c]) ) {
            BD_alt[r,c] <- 0
        }
    }
    if ( (sum(BD_alt[r,] == 0) == 2) & (max(BD_alt[r,],na.rm=TRUE) < thirdTHRESH) ) {
        BD_alt[r,] <- c(0,0,0)
    }
    
}
BD_alt$infected <- rowSums(BD_alt) >0

cbind(final_BD[,c("Bd_Run_1","Bd_Run_2","Bd_Average_Run_3","pval","sig")], alt_infect =BD_alt$infected)

#' What we see is that by using a significant threshold of p=0.1 (which is fairly relaxed), the poisson model method is less stringent
#' than the straight threshold method. I believe the poisson model method is likely more reliable since it is able to detect cases where infeciton
#' is truely low, but consistent. For example, there are some cases where all 3/3 samples were BD-positive but at low abundances; however, the
#' straight threshold model does not recognize it as a 'positive' since the abundances are below the individual threshold. I believe the poisson model
#' method uses more of the information in the qPCR methods than the straight threshold method.\
#' From here, we will add in the "expected" BD load (lambda) into the mapping file for each toad.

# Insert expected BD loads for each sample
final_BD <- cbind(samples_bd, final_BD)
mf$eBD <- 0
mf[match(final_BD$samples_bd, mf$X.SampleID),"eBD"] <- final_BD$lambda

##### ADJUSTING DATA: Make same order, make 'time' variable ######

#### Edit metadata/mapping file #####
toKeepMF <-  c("X.SampleID"
               ,"eBD" # infection levels
               ,"SPEC" # Species code
               ,"TREATMENT_GROUP" # whether it was pre or post; should rename
               ,"BD100KZSP_3122011" # whether or not individual was exposed to BD: will rename
               ,"ANONYMIZED_NAME" # individual frog ID
               , "PRE_POST_BD_NUM" # timepoint
               , aHeader
               
)
newNames <- c("SampleID"
              , "eBD_log"
              , "species"
              , "prepost"
              , "BD_infected"
              , "toadID"
              , "timepoint"
              , aMetrics
)

mf <- mf %>% # mapping file with controls still in it
    as_tibble() %>% # make into tibble
    dplyr::select(one_of(toKeepMF)) %>% # filter to only relevant variables
    rename_at(vars(toKeepMF), ~ newNames) %>% #rename variable names
    filter(prepost == "Pre" | prepost == "Pos") %>% # get rid of things in "tre", which is a different experiment
    separate(timepoint, into=c("exposure","time"), sep="_", convert=TRUE) %>% # create a time variable
    mutate( time = ifelse(exposure == "Pre", time, time + 5)) %>% # time variable "restarts" at BD exposure point, but I want it to be continuous
    mutate(species=ifelse(species=="Bubo","Anbo",ifelse(species=="Buma","Anma",ifelse(species=="Raca","Lica",ifelse(species=="Rapi","Lipi",species))))) %>%
    separate(toadID, into=c("todelete","ID"), sep="_",remove=TRUE) %>%
    unite(toadID, species,ID, sep = "_", remove = FALSE) %>%
    dplyr::select(-todelete, -ID)

#### Adjusting values in mf for Osse and Anbo ####

#' OSSE QUIRK\
#' \
#' So it turns out that Osse has this weird quirk where they weren't sampled in the first timepoint. Then,
#' in the mapping file they were given time points 1-4 instead of 2-5, which screws up how the sampling
#' lines up. So, I'm going to 
#' change the "pre" numbeers so that they line up nicely with the rest of the samples.

# rows to change
addTime1 <- which((mf$species=="Osse") & (mf$time %in% c(1,2,3,4,5)))
mf[addTime1, "time"] <- mf[addTime1, "time"] +1 # add 1 to all pre time points so that time point 5 lines up with all other species

#' BUBO QUIRK\
#' \
#' Additionally, there is this strange quirk with the Anbo dataset (or rather, strange quirk with every OTHER dataset) where nothing but
#' Anbo was sampled at time point 10. Thus, let's add +1 to all time point 10's for all samples except Anbo

addTime1_notanbo <- which((mf$species!="Anbo") & (mf$time == 15))
mf[addTime1_notanbo,"time"] <- mf[addTime1_notanbo,"time"] + 1

#### filtering mapping files
mf.wcon <- mf 
mf.tb <- mf.wcon %>%
    filter(species != "None") # get rid of controls

# Number of controls lost
c(nrow(mf.wcon)-nrow(mf.tb))

### Check to make sure there are no duplicates
#+ message=FALSE, echo=FALSE, warnings=FALSE, include=FALSE
any(table(dplyr::select(mf.tb, toadID, time)) > 1)
# Oh no! There is a duplicate! Check which one it is.
#+ message=FALSE, echo=FALSE, warnings=FALSE, include=FALSE
table(dplyr::select(mf.tb, toadID, time))>1
# Lica10 has 2 timepoint 3's
# Anma11 has 2 time point 3's
# I think this is just a mistake in the data entry. I've gone through the rest of the metadata and determined that 
#   it is probably mis-labeled so:
# I'm going to manually change:
# I think mck115Pre.3Buma.11.452405 Anma11 timepoint 3 is actually Anma9, time 3
# I think mck105Pre.3Raca.10.451190 Lica10 timepoint 3 is actually Lica9, time 3
mf.tb[mf.tb$SampleID == "mck115Pre.3Buma.11.452405",c("toadID")] <- "Anma_9"
mf.tb[mf.tb$SampleID == "mck105Pre.3Raca.10.451190", c("toadID")] <- "Lica_9"

#### Adjusting for individuals who were contaminated at the beginning

#' Finally, there were certain individuals who actually tested positive for BD when they arrived.
#' These were unknown until later because the PCR process takes a while, so they were included in the experiment.
#' However, let us identify these individuals and remove them. The list below was manually curated from the spreadsheet
#' provided by Val.
#' 

BD_contam_upon_arrival <- c("Anma_4"
                            , "Anma_6"
                            , "Anma_7"
                            , "Anma_10"
                            , "Anma_11"
                            , "Lipi_1"
                            , "Lipi_2"
                            , "Lipi_5"
                            , "Lipi_7"
                            , "Lipi_8"
                            , "Lipi_9"
                            , "Lipi_10"
                            , "Lipi_11"
                            , "Lipi_12")

# Add new column for individuals that were originally contamined
mf.tb$orig_contam <- 0
mf.tb[which(mf.tb$toadID %in% BD_contam_upon_arrival),"orig_contam"] <- 1

#### Other adjustments

# Make a PABD column as well as a "raw" BD counts column
mf.tb <- mf.tb %>%
    mutate(PABD=eBD_log>0, eBD_raw=exp(eBD_log)-1)

# # make diversity metrics numeric
mf.tb$shannon <- as.numeric(mf.tb$shannon)
mf.tb$observed_otus <- as.numeric(mf.tb$observed_otus)
mf.tb$logRich <- log(mf.tb$observed_otus)
# View(mf.tb %>%
         # dplyr::select(SampleID, shannon))

# filter to ONLY include those things in otu table (Raw)
mf.raw <- mf.tb %>%
    filter(SampleID %in% colnames(otu))
# How many things did we lose? Hopefully none, so that means the OTU table is complete.
c(nrow(mf.raw), nrow(mf.tb))
# Note: This ALSO means that controls were filtered out before making the OTU table.


# Samples to keep in otu table (raw)
keepSamples <- colnames(otu_rare)[-match(c("X.OTU.ID","taxonomy"),colnames(otu_rare))]
### ** Note, this assumes that there are no samples in OTU table that are NOT in mapping file!
# If you get random warnings later on, you should come back here.
OTU_names <- otu$X.OTU.ID
otu.filt <- otu %>%
    as_tibble() %>%
    dplyr::select(one_of(keepSamples)) %>% # Use rare because it naturally cuts off for <5000
    replace(.<minOTUSample, 0) %>%
    mutate(rowsums=rowSums(.)) %>%
    mutate(OTUID = OTU_names) %>%
    filter(rowsums > minOTUTab) %>%
    dplyr::select(-rowsums) 
# Samples to keep in otu table (rarefied)
OTU_names_rare <- otu_rare$X.OTU.ID
otu.filt_rare <- otu %>%
    as_tibble() %>%
    dplyr::select(one_of(keepSamples)) %>% # Use rare
    replace(.<minOTUSample, 0) %>%
    mutate(rowsums=rowSums(.)) %>%
    mutate(OTUID = OTU_names) %>%
    filter(rowsums > minOTUTab) %>%
    dplyr::select(-rowsums) 
# View(otu.filt_rare)
# Check that each have the same number of samples-- should be the same.
ncol(otu.filt_rare)
ncol(otu.filt)


#### The mapping file already has diversity in it, so I need to add beta diversity and inihibitory bacterial diversity

##### Beta diversity turnover #####
# Get distance to sample directly before for every sample
# AKA: How similar was the current sample to the sample JUST before that time point for that individual toad?
# Note; tried to do this with dplyr but it is VERY slow. Just used base R subsetting instead

for ( name in names(dm_all) ) {
    mf.raw[,paste0("distance_",name)] <- NA
    dm <- dm_all[[name]]
    dm <- data.frame(dm, row.names = 1)
    for ( i in 1:nrow(mf.raw) ) {
        currentSample <- mf.raw$SampleID[i]
        current <- mf.raw[i, c("toadID","time")]
        prevSample <- mf.raw$SampleID[mf.raw$toadID==as.character(current[1,1]) & mf.raw$time == as.numeric(current[1,2]-1)]
        if ( length(prevSample) > 0 & (currentSample %in% rownames(dm))) {
            distTemp <- dm[currentSample,as.character(prevSample)]
            if ( length(distTemp) > 0) {
                mf.raw[i,paste0("distance_",name)] <- distTemp
                
            }
        }
        
        # print(paste0("done",i," out of ",nrow(mf.tb)))
    }
}

## 


### Get proportion and count of inhibitory otus from inhibitory otu metadata
inhib.tb <- inhib %>%
    as_tibble() %>%
    rename(Name=V1, inhib=V2, num=V3, seq=V4,sampleID=V5, taxa=V6)
otu.tb.inhib <- otu.filt %>%
    mutate(inhib = (inhib.tb[match(otu.filt$OTUID, inhib.tb$seq),"inhib"])$inhib)
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
mf.raw <- mf.raw %>%
    mutate(n=sampleCounts[match(SampleID, names(sampleCounts))]
           ,inhibCounts =inhibCounts[match(SampleID, names(inhibCounts))]
           , inhibRich=inhibRich[match(SampleID, names(inhibRich))]) %>%
    mutate(percInhib = inhibCounts/n)
# Now, filter mf to rarefied OTU table
mf.rare <- mf.raw %>%
    filter(SampleID %in% colnames(otu_rare))


#### Adding NMDS ####
# MAKE NMDS
dm.filt <- dm[mf.rare$SampleID,mf.rare$SampleID]
# Make NMDS of distance matrix
set.seed(1017937)
nmds.all <- isoMDS(dist(dm.filt), k = 2)
nmds <- nmds.all$points
colnames(nmds) <- c("NMDS1","NMDS2")
mf.rare <- cbind(mf.rare, nmds)


#### Created filtered mapping files with all extra information
# Mapping file with only controls (training dataset)
# Save Rdata
save(mf.raw, file="mf.raw.RData")
save(mf.rare, file="mf.rare.RData")
# How many samples did we lose due to rarefaction?
c(nrow(mf.rare), nrow(mf.tb))
# Which samples did we lose due to rarefaction?
mf.tb %>%
    filter(!(SampleID %in% colnames(otu_rare))) %>%
    dplyr::select(SampleID) %>%
    pull()
# View(mf_con)
mf_con <- mf.rare %>%
    filter(BD_infected=="n")
# were any controls contaminated at any point?
con_contam <- mf_con %>%
    filter(PABD==TRUE) %>%
    dplyr::select(toadID, time) %>%
    group_by(toadID) %>%
    summarise(firstinfect=min(time))
# get all time points after that
toDel <- c()
for ( n in 1:length(con_contam$toadID)) {
    toad <- con_contam$toadID[n]
    minTime <- con_contam$firstinfect[n]
    temp <- mf_con %>%
        filter(toadID==toad, time>=minTime) %>%
        dplyr::select(SampleID) %>%
        pull()
    toDel <- c(toDel, temp)
}

# Get rid of post-infected ones that weren't supposed to be infected
mf_con_without_init_infect <- mf_con %>%
    filter(orig_contam == 0, !(SampleID %in% toDel))

# Mapping file with only treatment (testing dataset)
mf_treat <- mf.rare %>%
    filter(BD_infected=="y")
mf_treat_without_init_infect <- mf_treat %>%
    filter(orig_contam == 0)
save(mf_con_without_init_infect, file="mf_con_without_init_infect.RData")
save(mf_treat_without_init_infect, file="mf_treat_without_init_infect.RData")

#### Make OTU table of JUST inhibitory bacteria ####
# Should double check this at some point
temp_char_otu <- otu.tb.inhib %>%
    dplyr::select(c(OTUID, inhib))
# Make con and treat with just inhibitory with taxa names
otu.temp <- t(otu.tb.inhib %>%
                  dplyr::select(-c(OTUID, inhib)))
otu.temp <- cbind(t(otu.temp/sampleCounts), temp_char_otu)

otu.inhibOnly.con <- otu.temp %>%
    filter(inhib=="inhibitory") %>%
    dplyr::select(mf_con_without_init_infect$SampleID, "OTUID")
otu.inhibOnly.treat <- otu.temp %>%
    filter(inhib=="inhibitory") %>%
    dplyr::select(mf_treat_without_init_infect$SampleID, "OTUID") 

# save otu tables
save(otu.inhibOnly.con, file="otu.inhibOnly.con.RData")
save(otu.inhibOnly.treat, file="otu.inhibOnly.treat.RData")


#### PLOTTING EXP DESIGN ####
#+ fig.height=12, fig.width=5
mf.rare %>%
    separate(toadID, into=c("sp2", "indiv"), remove = FALSE) %>%
    mutate(indiv = factor(indiv, levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))) %>%
    mutate(Treatment=ifelse(BD_infected=="y","Bd-exposed","Control"), "LnBd_load" = eBD_log) %>%
    ggplot(aes(x=time, y=indiv)) +
    geom_line(aes(group=toadID, col=Treatment)) +
    geom_point(aes(group=toadID,bg=LnBd_load), cex=4, pch=21)+
    scale_color_manual(values=c("blue","orange")) +
    scale_fill_gradient(low = "white", high = "red") +
    geom_vline(aes(xintercept=5.5), col="orange")+
    facet_wrap(~species, nrow=5) +
    xlab("Time") +
    ylab("Individual Toad")

# Here, we show only infected individuals and their BD load; the individuals who were already infected are removed.
mf_treat_without_init_infect %>%
    separate(toadID, into=c("sp2", "indiv"), remove = FALSE) %>%
    mutate(indiv = factor(indiv, levels=c("1","2","3","4","5","6","7","8","9","10","11","12"))) %>%
    ggplot(aes(x=time, y=indiv)) +
    theme_bw() +
    geom_tile(aes(group=toadID,fill=eBD_log), lwd=0.5)+
    scale_color_manual(values=c("grey","black")) +
    geom_vline(aes(xintercept=5.5), col="orange")+
    facet_wrap(~species, nrow=5)

# What we learn is that no control individuals were infected at any point, after removing "pre-infected" individuals"

#### PLOTTING BETA PLOTS ####

#+
# What do different species look like? (Control only)
mf_con_without_init_infect %>%
    ggplot(aes(x=NMDS1, y=NMDS2)) +
    geom_point(aes(col=species), cex=3)
# What do different species look like? (ALL data)
mf.rare %>%
    ggplot(aes(x=NMDS1,y=NMDS2)) +
    geom_point(aes(col=species), cex=2)
# Color dots by time
mf_con_without_init_infect %>%
    ggplot(aes(x=NMDS1, y=NMDS2)) +
    geom_point(aes(col=time), cex=3)
# Filtering dm to include only controls
dm.filt.con <- dm.filt[mf_con_without_init_infect$SampleID,mf_con_without_init_infect$SampleID ]
save(dm.filt.con, file="dm.filt.con.RData")

#' There is a significant effect of species AND time AND interaction on COMPOSITION

rbind(mf_con_without_init_infect, mf_treat_without_init_infect) %>%
    mutate(BD_infected=ifelse(BD_infected=="y","Treatment","Control")
           , exposure = factor(exposure, levels=c("Pre","Post"))) %>%
    ggplot(aes(x=NMDS1,y=NMDS2)) +
    geom_point(aes(bg=exposure, col=PABD), cex=2, alpha=0.8, pch=21) +
    facet_grid(BD_infected ~ species) +
    scale_color_manual(values=c("white","black")) +
    scale_fill_manual(values=c("blue","red"))

mf_treat_without_init_infect %>%
    filter(exposure=="Post") %>%
    ggplot(aes(x=NMDS1, y=NMDS2)) +
    geom_point(aes(bg=time, col=PABD), cex=3, pch=21) +
    scale_color_manual(values=c("white","red"))+
    facet_wrap(~species, nrow=1)
mf_treat_without_init_infect_post <- mf_treat_without_init_infect %>%
    filter(exposure=="Post")
dm.filt.treat <- dm.filt[mf_treat_without_init_infect_post$SampleID, mf_treat_without_init_infect_post$SampleID]
save(dm.filt.treat, file="dm.filt.treat.RData")

#### ALPHA DIVERISTY ####

#' ### ALPHA DIVERSITY PLOTTING
#' \
#' First thing to do is calculated expected alpha diversity for each individual, and then
#' calculate how for that particular sample was from the "expected" diversity. Before
#' we do that though, let's plot alpha diversity to see how it looks.

# Fit normal and lognormal distributions to each of these
shannon_normal_param <- fitdistr(mf_con_without_init_infect$shannon, densfun="Normal")
x.fit.shannon <- seq(min(mf_con_without_init_infect$shannon)-sd(mf_con_without_init_infect$shannon)
                     , max(mf_con_without_init_infect$shannon)+sd(mf_con_without_init_infect$shannon)
                     , length.out = 100)
y.pred.shannon <- dnorm(x.fit.shannon, mean = shannon_normal_param$estimate[1], sd = shannon_normal_param$estimate[2])

obsotu_lognormal_param <- fitdistr(mf_con_without_init_infect$logRich, densfun="Normal")
x.fit.obsotu <- seq(min(mf_con_without_init_infect$logRich)-sd(mf_con_without_init_infect$logRich)
                     , max(mf_con_without_init_infect$logRich)+sd(mf_con_without_init_infect$logRich)
                     , length.out = 100)
y.pred.obsotu <- dnorm(x.fit.obsotu, mean = (obsotu_lognormal_param$estimate[1]), sd = (obsotu_lognormal_param$estimate[2]))

gg_shannon_all <- ggplot(data=mf_con_without_init_infect, aes(x=shannon)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.shannon, y=y.pred.shannon), aes(x=x, y=y), col="red")
gg_obsotu_all <- ggplot(data=mf_con_without_init_infect, aes(x=logRich)) +
    geom_histogram(aes(y=..density..), bins=20) +
    geom_line(data=data.frame(x=x.fit.obsotu, y=y.pred.obsotu), aes(x=x, y=y), col="red")
grid.arrange(gg_shannon_all, gg_obsotu_all, nrow=1)

# Check to see if turnover is changing with time significantly
gg_divtime_con <- mf_con_without_init_infect %>%
    filter(!is.na(shannon)) %>%
    ggplot(aes(x=time, y=shannon)) + 
    geom_line(aes(group=toadID)) +
    geom_point(aes(col=PABD)) +
    scale_color_manual(values="blue")+
    facet_grid(~species)
gg_divtime_treat <- mf_treat_without_init_infect %>%
    filter(!is.na(shannon)) %>%
    ggplot(aes(x=time, y=shannon)) + 
    geom_line(aes(group=toadID)) +
    geom_point(aes(group=toadID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5)) +
    facet_grid(~species)
grid.arrange(gg_divtime_con, gg_divtime_treat, nrow=2)

#' Does diversity change over time in control individuals?
# Type I ANOVA to test for interaction-- (AB | A, B)
anova(lm(shannon ~ species*time, data=mf_con_without_init_infect))
# Use Type II ANOVA (no interaction present)
Anova(lm(shannon ~ species + time, data=mf_con_without_init_infect), type = 2)
# There is a significant effect of species but not time or interaction

#' Does richness change over time in treatment individuals?
# Type I ANOVA to test for interaction (AB | A,B)
anova(lm(shannon ~ species*time, data=mf_treat_without_init_infect))
# Type III ANOVA (valid in presence of interaction)
Anova(lm(shannon ~ species * time, data=mf_treat_without_init_infect, contrasts=list(species=contr.sum)), type=3)


gg_richtime_con <- mf_con_without_init_infect %>%
    filter(!is.na(logRich)) %>%
    ggplot(aes(x=time, y=logRich)) + 
    geom_line(aes(group=toadID)) +
    geom_point(aes(group=toadID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_richtime_treat <- mf_treat_without_init_infect %>%
    filter(!is.na(logRich)) %>%
    ggplot(aes(x=time, y=logRich)) + 
    geom_line(aes(group=toadID)) +
    geom_point(aes(group=toadID, col=PABD)) +
    geom_vline(aes(xintercept=5.5)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
grid.arrange(gg_richtime_con, gg_richtime_treat, nrow=2)

#' Does richness change over time in control individuals?
# Type I ANOVA to test for interaction-- (AB | A, B)
anova(lm(logRich ~ species*time, data=mf_con_without_init_infect))
# Use Type II ANOVA (no interaction present)
Anova(lm(logRich ~ species + time, data=mf_con_without_init_infect), type = 2)
# There is a significant effect of species but not time or interaction

#' Does richness change over time in treatment individuals?
# Type I ANOVA to test for interaction (AB | A,B)
anova(lm(logRich ~ species*time, data=mf_treat_without_init_infect))
# Type III ANOVA (valid in presence of interaction)
Anova(lm(logRich ~ species * time, data=mf_treat_without_init_infect, contrasts=list(species=contr.sum)), type=3)


#### SHANNON ####
#' We see that shannon and log of observed otus look approximately normal. Great!
#' Now, let's fit some models to this data. For shannon and logRich, we should use a LMM
#' to find out the average richness for each species and the average richness for each individual through time.
#' \
#' u ~ N(u_i, sigma_i)\
#' u_i = a_j\
#' a_j ~ N(u_sp, sigma_sp)\
#' where i = sample, j = individual, sp = species
#' Below, we use the dataset with JUST the controls.

if ( FALSE ) {
    lmer_shannon <- stan_lmer(shannon ~ -1 + species + (1|toadID), data=mf_con_without_init_infect
                              , prior = normal(0, 10, autoscale = TRUE)
                              , seed = 98374)
    save(lmer_shannon, file="lmer_shannon.RData")
} else {
    load("lmer_shannon.RData")
}
prior_summary(lmer_shannon)
# Look at distributions according to models
samps_lmer_shannon <- rstan::extract(lmer_shannon$stanfit)
pre_test_set <- mf_treat_without_init_infect %>%
    filter(time<6) 
## PLOT OF EXPECTED RCIHNESS FOR EACH SPECIES
samps_lmer_shannon$beta %>%
    as.data.frame() %>%
    rename(Anbo=V1, Anma=V2, Lica=V3, Lipi=V4, Osse=V5) %>%
    dplyr::select(Anbo,Anma,Lica,Lipi,Osse) %>%
    gather(key=species, value=shannon) %>%
    ggplot(mapping=aes(x=species, y=shannon))+
    geom_violin() +
    geom_point(data=mf_con_without_init_infect, aes(y=shannon, x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
    geom_point(data=pre_test_set, aes(y=shannon, x=species), position=position_jitter(width = 0.1, height=0), col="red")

# Get standard deviation between toad individuals and samples
samp_sigma <- samps_lmer_shannon$aux
toadID_sigma <- sd(samps_lmer_shannon$b[,ncol(samps_lmer_shannon$b)])

# Now, we can calculate the probability that the "test" dataset values come from this distribution
# List of individuals
treat_indiv <- unique(mf_treat_without_init_infect$toadID)
# List of each species
species_list <- levels(factor(mf_con_without_init_infect$species))

exp_distr <- as.data.frame(matrix(ncol=length(species_list), nrow=4000, dimnames = list(1:4000, species_list)))
for ( num_sp in 1:length(species_list)) {
    exp_distr[,num_sp] <- rnorm(length(samps_lmer_shannon$beta[,num_sp]), mean=samps_lmer_shannon$beta[,num_sp], sd=toadID_sigma)
}

# Loop through and calculate probability of having diversity at that level
pre_exp_indiv <- data.frame(toadID=treat_indiv, exp_shan=rep(NA, length(treat_indiv)), p_shan=rep(NA, length(treat_indiv)), infect=rep(NA, length(treat_indiv)))
for ( i in treat_indiv ) {
    n_row <- match(i, treat_indiv)
    sp <- unlist(strsplit(i,"_"))
    num_sp <- match(sp[1], levels(factor(mf_con_without_init_infect$species)))
    temp_shan <- mf_treat_without_init_infect %>%
        filter(toadID==i, time <=5 ) %>%
        dplyr::select(shannon) %>%
        pull()
    if ( length(temp_shan)>1 ) {
        exp_shan <-  fitdistr(temp_shan, "normal")$estimate[1]
        
    } else {
        exp_shan <- temp_shan
    }
    
    # pred_distr <- rnorm(length(samps_lmer_shannon$beta[,num_sp]), mean=rnorm(length(samps_lmer_shannon$beta[,num_sp]), mean=samps_lmer_shannon$beta[,num_sp], sd=toadID_sigma), sd=samp_sigma)
    p_shan <- sum(exp_distr[,sp[1]]<exp_shan)/length(exp_distr[,sp[1]])
    
    ### Did they get infected?
    infect <- max(mf_treat_without_init_infect %>%
        filter(toadID==i) %>%
        dplyr::select(eBD_raw) %>%
        pull()
        )
    
    pre_exp_indiv[n_row,c("exp_shan","p_shan","infect")] <- c(exp_shan, p_shan, infect)
    
}
# create species column
pre_exp_indiv <- pre_exp_indiv %>%
    separate(toadID, into=c("species","indiv"), remove = FALSE)

# Plot results 
gg_shan_p <- ggplot(pre_exp_indiv, aes(x=p_shan, y=log(infect+1))) +
    geom_smooth(aes(color=species),method=lm, se = FALSE) +
    geom_point(aes(color=species), cex=4) +
    geom_smooth(method=lm, se=FALSE, col="black")
# if we'd JUST plotted raw values
gg_shan_raw <- ggplot(pre_exp_indiv, aes(x=exp_shan, y=log(infect+1)))+
    geom_point(aes(color=species), cex=4) +
    geom_smooth(method=lm, se=FALSE) +
    geom_smooth(aes(color=species),method=lm, se = FALSE) 
grid.arrange(gg_shan_p, gg_shan_raw, nrow=1)
exp_distr %>%
    gather(key=species, value=shannon) %>%
    ggplot(aes(x=species, y=shannon)) +
    geom_violin() +
    geom_point(data=pre_exp_indiv, aes(x=species, y=exp_shan, col=log(infect+1)), cex=4, position=position_jitter(height=0, width=0.1))

all_p <- pre_exp_indiv %>%
    dplyr::select(toadID, exp_shan, p_shan,infect)

#### RICHNESS (observed OTUs) ####

if (FALSE) {
    lmer_rich <- stan_lmer(logRich ~ -1 + species + (1|toadID), data=mf_con_without_init_infect
                           , prior = normal(0, 10, autoscale = TRUE)
                           , seed = 98374)
    save(lmer_rich, file="lmer_rich.RData")
} else {
    load(file="lmer_rich.RData")
}
prior_summary(lmer_rich)

# Look at distributions according to models
samps_lmer_rich <- rstan::extract(lmer_rich$stanfit)
pre_test_set <- mf_treat_without_init_infect %>%
    filter(time<6) 
# Plot of expected by species, and real samples
samps_lmer_rich$beta %>%
    as.data.frame() %>%
    rename(Anbo=V1, Anma=V2, Lica=V3, Lipi=V4, Osse=V5) %>%
    dplyr::select(Anbo,Anma,Lica,Lipi,Osse) %>%
    gather(key=species, value=logRich) %>%
    ggplot(mapping=aes(x=species, y=logRich))+
    geom_violin() +
    geom_point(data=mf_con_without_init_infect, aes(y=logRich, x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
    geom_point(data=pre_test_set, aes(y=logRich, x=species), position=position_jitter(width = 0.1, height=0), col="red")

# Get standard deviation between toad individuals and samples
samp_sigma <- samps_lmer_rich$aux
toadID_sigma <- sd(samps_lmer_rich$b[,ncol(samps_lmer_rich$b)])

# Get estimates for control toads
rich_species_exp <- fixef(lmer_rich)  %>%
    as_tibble() %>%
    mutate(species=c("Anbo", "Anma","Lica","Lipi","Osse"))
con_toad_est_rich <- ranef(lmer_rich)$toadID %>%
    rename(sp_est="(Intercept)") %>%
    mutate(toadID=rownames(ranef(lmer_rich)$toadID)) %>%
    separate(toadID, into=c("species","num"), sep="_",remove=FALSE) %>%
    dplyr::select(-num) %>%
    left_join(rich_species_exp, by = "species") %>%
    mutate(est_rich=sp_est+value)

# Now, we can calculate the probability that the "test" dataset values come from this distribution
# List of individuals
treat_indiv <- unique(mf_treat_without_init_infect$toadID)
# List of each species
species_list <- levels(factor(mf_con_without_init_infect$species))

exp_distr <- as.data.frame(matrix(ncol=length(species_list), nrow=4000, dimnames = list(1:4000, species_list)))
for ( num_sp in 1:length(species_list)) {
    exp_distr[,num_sp] <- rnorm(length(samps_lmer_rich$beta[,num_sp]), mean=samps_lmer_rich$beta[,num_sp], sd=toadID_sigma)
    # exp_distr[,num_sp] <- rnorm(length(samps_lmer_rich$beta[,num_sp]), mean=rnorm(length(samps_lmer_rich$beta[,num_sp]), mean=samps_lmer_rich$beta[,num_sp], sd=toadID_sigma), sd=samp_sigma)
}

# Loop through and calculate probability of having diversity at that level
pre_exp_indiv <- data.frame(toadID=treat_indiv, exp_rich=rep(NA, length(treat_indiv)), p_rich=rep(NA, length(treat_indiv)), infect=rep(NA, length(treat_indiv)))
for ( i in treat_indiv ) {
    n_row <- match(i, treat_indiv)
    sp <- unlist(strsplit(i,"_"))
    num_sp <- match(sp[1], levels(factor(mf_con_without_init_infect$species)))
    temp_rich <- mf_treat_without_init_infect %>%
        filter(toadID==i, time <=5 ) %>%
        dplyr::select(logRich) %>%
        pull()
    if ( length(temp_rich)>1 ) {
        exp_rich <-  fitdistr(temp_rich, "normal")$estimate[1]
        
    } else {
        exp_rich <- temp_rich
    }
    
    # pred_distr <- rnorm(length(samps_lmer_rich$beta[,num_sp]), mean=rnorm(length(samps_lmer_rich$beta[,num_sp]), mean=samps_lmer_rich$beta[,num_sp], sd=toadID_sigma), sd=samp_sigma)
    p_rich <- sum(exp_distr[,sp[1]]<exp_rich)/length(exp_distr[,sp[1]])
    
    ### Did they get infected?
    infect <- max(mf_treat_without_init_infect %>%
                      filter(toadID==i) %>%
                      dplyr::select(eBD_raw) %>%
                      pull()
    )
    
    pre_exp_indiv[n_row,c("exp_rich","p_rich","infect")] <- c(exp_rich, p_rich, infect)
    
}

con_exp_indiv_rich <- data.frame(toadID=con_toad_est_rich$toadID, exp_rich=rep(NA, length(con_toad_est_rich$toadID)), p_rich=rep(NA, length(con_toad_est_rich$toadID)))
for ( i in 1:nrow(con_toad_est_rich) ) {
    s <- con_toad_est_rich[i,"species"]
    
    exp_rich <- con_toad_est_rich[i,"est_rich"]
    p_rich <- sum(exp_distr[,s]<exp_rich)/length(exp_distr[,s])
    
    con_exp_indiv_rich[i,c("exp_rich","p_rich")] <- c(exp_rich, p_rich)
}
# create species column
pre_exp_indiv <- pre_exp_indiv %>%
    separate(toadID, into=c("species","indiv"), remove = FALSE)

# Plot results 
gg_logRich_p <- ggplot(pre_exp_indiv, aes(x=p_rich, y=log(infect+1))) +
    geom_smooth(aes(color=species),method=lm, se = FALSE) +
    geom_point(aes(color=species), cex=4) +
    geom_smooth(method=lm, se=FALSE, col="black")
# if we'd JUST plotted raw values
gg_logRich_raw <- ggplot(pre_exp_indiv, aes(x=exp_rich, y=log(infect+1)))+
    geom_point(aes(color=species), cex=4) +
    geom_smooth(aes(color=species),method=lm, se = FALSE) +
    geom_smooth(method=lm, se=FALSE)
grid.arrange(gg_logRich_p, gg_logRich_raw, nrow=1)
exp_distr %>%
    gather(key=species, value=logRich) %>%
    ggplot(aes(x=species, y=logRich)) +
    geom_violin() +
    geom_point(data=pre_exp_indiv, aes(x=species, y=exp_rich, col=log(infect+1)), cex=4, position=position_jitter(height=0, width=0.1))

all_p <- pre_exp_indiv %>%
    dplyr::select(toadID, exp_rich, p_rich) %>%
    full_join(all_p, by="toadID")

#### BETA DIVERSTY ####

#### Beta plot
# First, fit a beta distribution
x.fit.beta <- seq(min(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)-sd(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)
                  , max(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)+sd(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)
                  , length.out = 100)
beta.fit <- fitdistr(mf_con_without_init_infect$distance_bray_curtis[!is.na(mf_con_without_init_infect$distance_bray_curtis)]
, densfun="beta", start=list(shape1=6,shape2=6))
y.pred.beta <- dbeta(x.fit.beta, shape1 = beta.fit$estimate[1], shape2 = beta.fit$estimate[2])
# Try a normal too?
beta.fit.norm <- fitdistr(mf_con_without_init_infect$distance_bray_curtis[!is.na(mf_con_without_init_infect$distance_bray_curtis)]
                     , densfun="Normal")
y.pred.beta.norm <- dnorm(x.fit.beta, mean = beta.fit.norm$estimate[1], sd = beta.fit.norm$estimate[2])

mf_con_without_init_infect %>%
    filter(!is.na(distance_bray_curtis)) %>%
    ggplot(aes(x=distance_bray_curtis)) + 
    geom_histogram(aes(y=..density..), bins=25) +
    geom_line(data=data.frame(x=x.fit.beta, y=y.pred.beta), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.fit.beta, y=y.pred.beta.norm), aes(x=x, y=y), col="blue") 

# Check to see if turnover is changing with time significantly
gg_betatime_con <- mf_con_without_init_infect %>%
    filter(!is.na(distance_bray_curtis)) %>%
    ggplot(aes(x=time, y=distance_bray_curtis)) + 
    geom_line(aes(group=toadID)) +
    geom_point(aes(group=toadID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_betatime_treat <- mf_treat_without_init_infect %>%
    filter(!is.na(distance_bray_curtis)) %>%
    ggplot(aes(x=time, y=distance_bray_curtis)) + 
    geom_line(aes(group=toadID)) +
    geom_point(aes(group=toadID, col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    geom_vline(aes(xintercept=5.5))+
    facet_grid(~species)
grid.arrange(gg_betatime_con, gg_betatime_treat, nrow=2)

# Is there an effect of species and time on controls?
# Type I ANOVA (to check for interaction) (AB | A,B)
anova(lm(distance_bray_curtis ~ species*time, data=mf_con_without_init_infect))
# Type II ANOVA with no interaction
Anova(lm(distance_bray_curtis ~ species + time, data=mf_con_without_init_infect), type = 2)

# Is there an effect of species and time on treatment??
# Type I ANOVA (to check for interaction) (AB | A,B)
anova(lm(distance_bray_curtis ~ species*time, data=mf_treat_without_init_infect))
# Type II ANOVA with no interaction
Anova(lm(distance_bray_curtis ~ species + time, data=mf_treat_without_init_infect), type = 2)


#' We see that beta diversity is fairly normal, but we probably want to use a beta distribution since it's continuous and bound between
#' 0 and 1.
#' Now, let's fit some models to this data. We should use a GLMM with beta distribution as the response variable
#' to find out the average beta diversity turnover for each species and for each individual through time.
#' \
#' u ~ beta(u_i, sigma_i)\
#' u_i = a_j\
#' a_j ~ N(u_sp, sigma_sp)\
#' where i = sample, j = individual, sp = species
#' Below, we use the dataset with JUST the controls.

# inhibBin <- cbind(inhibCount=mf_con_without_init_infect$inhibCounts, noninhibCount=mf_con_without_init_infect$n)

    
if ( FALSE) {
    glmer_BC <- stan_glmer(distance_bray_curtis ~ -1 + species + (1|toadID)
                           , data=mf_con_without_init_infect
                           , family =mgcv::betar
                           , prior_intercept = normal(location = 0.5,scale = 0.1, autoscale = TRUE)
                           , prior = normal(location=0.5, scale=0.1, autoscale=TRUE)
                           , seed= 623445
    )
    save(glmer_BC, file="glmer_BC.RData")
    } else {
        load("glmer_BC.RData")
    }
prior_summary(glmer_BC)

# rbeta has a strange parameterization using a nd b so need to convert mu and phi to this.
a <- function(mu,phi){
    mu*phi
}
b <- function(mu,phi) {
    phi-mu*phi
}
mu <- function(a,phi) {
    a/phi
}
inv_logit <- function(x) {
    exp(x)/(exp(x)+1)
}
logit <- function(p) {
    log(p/(1-p))
}


# Look at distributions according to models
samps_glmer_BC<- rstan::extract(glmer_BC$stanfit)
pre_test_set <- mf_treat_without_init_infect %>%
    filter(time<=5) 
samps_glmer_BC$beta %>%
    as.data.frame() %>%
    rename(Anbo=V1, Anma=V2, Lica=V3, Lipi=V4, Osse=V5) %>%
    mutate(Anbo=inv_logit(mu(Anbo, samps_glmer_BC$aux))
           ,Anma=inv_logit(mu(Anma, samps_glmer_BC$aux))
           ,Osse=inv_logit(mu(Osse, samps_glmer_BC$aux))
           ,Lica=inv_logit(mu(Lica, samps_glmer_BC$aux))
           ,Lipi=inv_logit(mu(Lipi, samps_glmer_BC$aux))) %>%
    dplyr::select(Anbo,Anma,Osse,Lica,Lipi) %>%
    gather(key=species, value=distance_bray_curtis) %>%
    ggplot(mapping=aes(x=species, y=distance_bray_curtis))+
    geom_violin() +
    geom_point(data=mf_con_without_init_infect, aes(y=distance_bray_curtis, x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
    geom_point(data=pre_test_set, aes(y=distance_bray_curtis, x=species), position=position_jitter(width = 0.1, height=0), col="red")

# Get standard deviation between toad individuals and samples
# samp_sigma <- sigma(glmer_BC)
toadID_sigma <- sd(samps_glmer_BC$b[,ncol(samps_glmer_BC$b)])
phi <- samps_glmer_BC$aux

# Now, we can calculate the probability that the "test" dataset values come from this distribution
# List of individuals
treat_indiv <- unique(mf_treat_without_init_infect$toadID)
# List of each species
species_list <- levels(factor(mf_con_without_init_infect$species))

exp_distr <- as.data.frame(matrix(ncol=length(species_list), nrow=4000, dimnames = list(1:4000, species_list)))
for ( num_sp in 1:length(species_list)) {
    mu <- inv_logit(rnorm(4000, mean=samps_glmer_BC$beta[,num_sp], sd=toadID_sigma))
    # exp_distr[,nump_sp] <- mu
    exp_distr[,num_sp] <- rbeta(length(samps_glmer_BC$beta[,num_sp])
                                ,shape1=a(mu,samps_glmer_BC$aux)
                                ,shape2=b(mu,samps_glmer_BC$aux))
    
    # exp_distr[,num_sp] <- rnorm(length(samps_glmer_BC$beta[,num_sp]), mean=rnorm(length(samps_glmer_BC$beta[,num_sp]), mean=samps_glmer_BC$beta[,num_sp], sd=toadID_sigma), sd=samp_sigma)
}


# Now, we can calculate the probability that the "test" dataset values come from this distribution
# List of individuals
treat_indiv <- unique(mf_treat_without_init_infect$toadID)
# Loop through and calculate probability of having diversity at that level
pre_exp_indiv <- data.frame(toadID=treat_indiv, exp_mu=rep(NA, length(treat_indiv)), p_mu=rep(NA, length(treat_indiv)), infect=rep(NA, length(treat_indiv)))
for ( i in treat_indiv ) {
    n_row <- match(i, treat_indiv)
    sp <- unlist(strsplit(i,"_"))
    num_sp <- match(sp[1], levels(factor(mf_con_without_init_infect$species)))
    temp_bc <- mf_treat_without_init_infect %>%
        filter(toadID==i, time <=5 ) %>%
        filter(!is.na(distance_bray_curtis))%>%
        filter(!is.na(n))%>%
        dplyr::select(distance_bray_curtis) %>%
        pull()
    
    if ( length(temp_bc) > 1) {
        exp_mu <-  fitdistr(temp_bc, "normal")$estimate[1]
        # exp_mu <-  fitdistr(temp_bc, "beta", start=list(shape1=0.1, shape2=0.1))$estimate[1]
    } else if ( length(temp_bc) == 1) {
        exp_mu <- temp_bc
    } else {
        exp_mu <- NA
    }
    
    # dm is distance matrix; larger exp_mu means more dissimilar. We want to know if MORE dissimilar == MORE infection
    p_mu <- sum(exp_distr[,sp[1]]<exp_mu, na.rm=TRUE)/length(exp_distr[,sp[1]])

    ### Did they get infected?
    infect <- max(mf_treat_without_init_infect %>%
                      filter(toadID==i) %>%
                      dplyr::select(eBD_raw) %>%
                      pull()
    )

    pre_exp_indiv[n_row,c("exp_mu","p_mu","infect")] <- c(exp_mu, p_mu, infect)
    
}
# create species column
pre_exp_indiv <- pre_exp_indiv %>%
    separate(toadID, into=c("species","indiv"), remove = FALSE)

# Plot results 
gg_beta_p <- pre_exp_indiv %>%
    filter(!is.na(exp_mu)) %>%
    ggplot(aes(x=p_mu, y=log(infect+1))) +
    geom_smooth(method=lm, se=FALSE) +
    geom_smooth(aes(col=species), method=lm, se=FALSE)+
    geom_point(aes(color=species), cex=4) 
# Raw numbers
gg_beta_raw <- pre_exp_indiv %>%
    filter(!is.na(exp_mu)) %>%
    ggplot(aes(x=exp_mu, y=log(infect+1))) +
    geom_smooth(method=lm, se=FALSE) +
    geom_smooth(aes(col=species), method=lm, se=FALSE)+
    geom_point(aes(color=species), cex=4)
grid.arrange(gg_beta_p, gg_beta_raw, nrow=1)
all_p <- pre_exp_indiv %>%
    dplyr::select(toadID, exp_mu, p_mu) %>%
    full_join(all_p, by="toadID")

#### PERCENT INHIBITORY ####
x.fit.percInhib <- round(seq(max(min(mf_con_without_init_infect$percInhib, na.rm=TRUE)-sd(mf_con_without_init_infect$percInhib, na.rm=TRUE),0)
                       , max(mf_con_without_init_infect$percInhib, na.rm=TRUE)+sd(mf_con_without_init_infect$percInhib, na.rm=TRUE), length.out = 100), digits = 2)
percInhib.fit <- fitdistr(mf_con_without_init_infect$percInhib[!is.na(mf_con_without_init_infect$percInhib)], densfun = "beta", start=list(shape1=1, shape2=6))
# percInhib.fit <- fitdistr(mf_con_without_init_infect$percInhib[!is.na(mf_con_without_init_infect$percInhib)], densfun = "gamma", start=list(shape=1, rate=8))
y.pred.percInhib <- dbeta(x.fit.percInhib, shape1=percInhib.fit$estimate[1], shape2=percInhib.fit$estimate[2])
# y.pred.percInhib <- dgamma(x.fit.percInhib, shape=percInhib.fit$estimate[1], rate=percInhib.fit$estimate[2])

mf_con_without_init_infect %>%
    filter(!is.na(percInhib)) %>%
    ggplot(aes(y=..density.., x=percInhib)) + 
    geom_histogram(bins=20) +
    geom_line(data=data.frame(x=x.fit.percInhib, y=y.pred.percInhib), aes(x=x,y=y), col="red")
mf_con_without_init_infect %>%
    filter(!is.na(percInhib)) %>%
    ggplot(aes(x=percInhib)) + 
    geom_histogram(aes(color=species), bins=20,show.legend = FALSE) +
    facet_grid(~species)
# Check to see if turnover is changing with time significantly
gg_perctime_con <- mf_con_without_init_infect %>%
    filter(!is.na(percInhib)) %>%
    ggplot(aes(x=time, y=percInhib)) + 
    geom_line(aes(group=toadID)) +
    geom_point(aes(col=PABD)) +
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_perctime_treat <- mf_treat_without_init_infect %>%
    filter(!is.na(percInhib)) %>%
    ggplot(aes(x=time, y=percInhib)) + 
    geom_line(aes(group=toadID)) +
    geom_point(aes(col=PABD)) +
    geom_vline(aes(xintercept=5.5))+
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
grid.arrange(gg_perctime_con, gg_perctime_treat, nrow=2)

# Does percent inhibitory change with species or time?
# Type I ANOVA (to test for interaction) in control group?
anova(glm(percInhib ~ species*time, family = binomial(), data=mf_con_without_init_infect, weights=mf_con_without_init_infect$n), test = "Chisq")
# Type III ANOVA (to test for main effects, given interaction) in control group?
Anova(glm(percInhib ~ species*time, family = binomial(), data=mf_con_without_init_infect, weights=mf_con_without_init_infect$n, contrasts=list(species=contr.sum)), type=3)

# Does percent inhibitory change with species or time in treatment group?
# Type I ANOVA (to test for interaction) in control group?
anova(glm(percInhib ~ species*time, family = binomial(), data=mf_treat_without_init_infect, weights=mf_treat_without_init_infect$n), test = "Chisq")
# Type III ANOVA (to test for main effects, given interaction) in control group?
Anova(glm(percInhib ~ species*time, family = binomial(), data=mf_treat_without_init_infect, weights=mf_treat_without_init_infect$n, contrasts = list(species=contr.sum)), type=3)

#' We see that beta diversity is fairly normal, but we probably want to use binomial since it's a proportion.
#' Now, let's fit some models to this data. We should use a GLMM with binomial as the response variable
#' to find out the average beta diversity turnover for each species and for each individual through time.
#' \
#' u ~ BIN(u_i, sigma_i)\
#' u_i = a_j\
#' a_j ~ N(u_sp, sigma_sp)\
#' where i = sample, j = individual, sp = species
#' Below, we use the dataset with JUST the controls.

inhibBin <- cbind(inhibCount=mf_con_without_init_infect$inhibCounts, noninhibCount=mf_con_without_init_infect$n)

if ( FALSE) {
    glmer_percInhib <- stan_glmer(inhibBin ~ -1 + species + (1|toadID)
                           , data=mf_con_without_init_infect
                           , family = 'binomial'
                           , seed= 9837423)
    save(glmer_percInhib, file="glmer_percInhib.RData")
} else {
    load("glmer_percInhib.RData")
}
prior_summary(glmer_percInhib)

# Look at distributions according to models
samps_glmer_percInhib<- rstan::extract(glmer_percInhib$stanfit)
pre_test_set <- mf_treat_without_init_infect %>%
    filter(time<=5) 
samps_glmer_percInhib$beta %>%
    as.data.frame() %>%
    rename(Anbo=V1, Anma=V2, Lica=V3, Lipi=V4, Osse=V5) %>%
    mutate(Anbo=inv_logit(Anbo)
           ,Anma=inv_logit(Anma)
           ,Osse=inv_logit(Osse)
           ,Lica=inv_logit(Lica)
           ,Lipi=inv_logit(Lipi)) %>%
    dplyr::select(Anbo,Anma,Lica,Lipi, Osse) %>%
    gather(key=species, value=percInhib) %>%
    ggplot(mapping=aes(x=species, y=percInhib))+
    geom_violin() +
    geom_point(data=mf_con_without_init_infect, aes(y=percInhib, x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
    geom_point(data=pre_test_set, aes(y=percInhib, x=species), position=position_jitter(width = 0.1, height=0), col="red")
# legend("topright", legend=c("Control (all)","Treatment (Pre)"), pch=21, col=c("blue","red"))



# indiv_mu <- ranef(glmer_percInhib)$toadID
# sp_mu <- fixef(glmer_BC)
toadID_sigma <- sd(samps_glmer_percInhib$b[,ncol(samps_glmer_percInhib$b)])

# instead of predicted distr, we get predicted mu 
# List of each species
species_list <- levels(factor(mf_con_without_init_infect$species))
# Get predicted distribuion
mu_exp_distr <- as.data.frame(matrix(ncol=length(species_list), nrow=4000, dimnames = list(1:4000, species_list)))
for ( num_sp in 1:length(species_list)) {
    # mu_exp_distr[,num_sp] <- rbinom(4000, mean=samps_glmer_percInhib$beta[,num_sp], sd=toadID_sigma)
    mu_exp_distr[,num_sp] <- rnorm(4000, mean=samps_glmer_percInhib$beta[,num_sp], sd=toadID_sigma)
}


# Loop through and calculate probability of having diversity at that level
pre_exp_indiv <- data.frame(toadID=treat_indiv, exp_pinhib=rep(NA, length(treat_indiv)), p_pinhib=rep(NA, length(treat_indiv)), infect=rep(NA, length(treat_indiv)))
for ( i in treat_indiv ) {
    n_row <- match(i, treat_indiv)
    sp <- unlist(strsplit(i,"_"))
    num_sp <- match(sp[1], levels(factor(mf_con_without_init_infect$species)))
    temp_p <- mf_treat_without_init_infect %>%
        filter(toadID==i, time <=5 ) %>%
        filter(!is.na(percInhib))%>%
        filter(!is.na(n))%>%
        dplyr::select(percInhib) %>%
        pull()
    temp_size <- mf_treat_without_init_infect %>%
        filter(toadID==i, time <=5 ) %>%
        filter(!is.na(percInhib))%>%
        filter(!is.na(n))%>%
        dplyr::select(n) %>%
        pull()
    
        exp_pinhib <- glm(cbind(temp_p*temp_size, temp_size-temp_p*temp_size) ~ 1, family="binomial")$coefficients
        p_pinhib <- sum(mu_exp_distr[,sp[1]]<exp_pinhib)/length(mu_exp_distr[,sp[1]])
        
        ### Did they get infected?
        infect <- max(mf_treat_without_init_infect %>%
                          filter(toadID==i) %>%
                          dplyr::select(eBD_raw) %>%
                          pull()
        )
        
        pre_exp_indiv[n_row,c("exp_pinhib","p_pinhib","infect")] <- c(exp_pinhib, p_pinhib, infect)

    
}
# create species column
pre_exp_indiv <- pre_exp_indiv %>%
    separate(toadID, into=c("species","indiv"), remove = FALSE)

# Plot results 
gg_percInhib_p <- ggplot(pre_exp_indiv, aes(x=p_pinhib, y=log(infect+1))) +
    geom_smooth(aes(col=species),method=lm, se=FALSE) +
    geom_point(aes(color=species), cex=4)+
    geom_smooth(method=lm, se=FALSE, col="black")
gg_percInhib_raw <- ggplot(pre_exp_indiv, aes(x=inv_logit(exp_pinhib), y=log(infect+1)))+
    geom_smooth(aes(col=species),method=lm, se=FALSE) +
    geom_point(aes(color=species), cex=4)+
    geom_smooth(method=lm, se=FALSE, col="black")
grid.arrange(gg_percInhib_p, gg_percInhib_raw, nrow=1)

mu_exp_distr %>%
    gather(key=species, value=percInhib) %>%
    ggplot(aes(x=species, y=percInhib)) +
    geom_violin() +
    geom_point(data=pre_exp_indiv, aes(x=species, y=exp_pinhib, col=log(infect+1)), cex=4, position=position_jitter(height=0, width=0.1))

all_p <- pre_exp_indiv %>%
    dplyr::select(toadID, exp_pinhib, p_pinhib) %>%
    full_join(all_p, by="toadID")

#### INHIB RICHNESS ####

### Fit a poisson and log poisson to see fit
x.fit.inhibRich <- round(seq(0,max((mf_con_without_init_infect$inhibRich), na.rm=TRUE)+sd((mf_con_without_init_infect$inhibRich), na.rm=TRUE), length.out = 100))

inhibRich.fit <- fitdistr((mf_con_without_init_infect$inhibRich), densfun = "Poisson")
y.pred.inhibRich <- dpois(x.fit.inhibRich, lambda = (inhibRich.fit$estimate[1]))

ggplot(data=mf_con_without_init_infect, aes(x=(inhibRich))) +
    geom_histogram(aes(y=..density..), bins=20, show.legend = FALSE) +
    geom_line(data=data.frame(x=(x.fit.inhibRich), y=y.pred.inhibRich), aes(x=x,y=y), col="red")


# Check to see if turnover is changing with time significantly
gg_inhibtime_con <- mf_con_without_init_infect %>%
    filter(!is.na(inhibRich)) %>%
    ggplot(aes(x=time, y=inhibRich)) + 
    geom_line(aes(group=toadID)) +
    geom_point(aes(col=PABD))+
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
gg_inhibtime_treat <- mf_treat_without_init_infect %>%
    filter(!is.na(inhibRich)) %>%
    ggplot(aes(x=time, y=inhibRich)) + 
    geom_line(aes(group=toadID)) +
    geom_point(aes(col=PABD))+
    geom_vline(aes(xintercept=5.5))+
    scale_color_manual(values=c("blue","red"))+
    facet_grid(~species)
grid.arrange(gg_inhibtime_con, gg_inhibtime_treat, nrow=2)

# Does proportion of inhibitory bacteria differ betwen species and time points?
# Type I ANOVA to test for interactions in control
anova(glm(inhibRich ~ species*time, data=mf_con_without_init_infect, family=poisson()), test="Chisq")
# TYpe III ANOVA to test for main effects with interactions in control
Anova(glm(inhibRich ~ species*time, data=mf_con_without_init_infect, family=poisson(), contrasts=list(species=contr.sum)),type=3)

# Type I ANOVA to test for interactions
anova(glm(inhibRich ~ species*time, data=mf_treat_without_init_infect, family=poisson()), test="Chisq")
# TYpe III ANOVA to test for main effects with interactions
Anova(glm(inhibRich ~ species*time, data=mf_treat_without_init_infect, family=poisson(), contrasts = list(species=contr.sum)),type=3)

if (FALSE) {
    glmer_inhibRich <- stan_glmer(inhibRich ~ species + (1|toadID), data=mf_con_without_init_infect
                                  , prior = normal(0, 10, autoscale = TRUE)
                                  , family= poisson(link="identity")
                                  , seed = 5423409)
    save(glmer_inhibRich, file="glmer_inhibRich.RData")
} else {
    load(file="glmer_inhibRich.RData")

}
prior_summary(glmer_inhibRich)

# Look at distributions according to models
samps_glmer_inhibRich <- rstan::extract(glmer_inhibRich$stanfit)
pre_test_set <- mf_treat_without_init_infect %>%
    filter(time<6) 

new_samps_beta <- samps_glmer_inhibRich$beta %>%
    as.data.frame() %>%
    mutate(V0=samps_glmer_inhibRich$alpha[,1]) %>%
    transmute(Anbo=V0, Anma=V0+V1, Lica=V0+V2, Lipi=V0+V3, Osse=V0+V4) %>%
    dplyr::select(Anbo,Anma,Lica,Lipi,Osse)
new_samps_beta %>% 
    gather(key=species, value=inhibRich) %>%
    ggplot(mapping=aes(x=species, y=inhibRich))+
    geom_violin() +
    geom_point(data=mf_con_without_init_infect, aes(y=(inhibRich), x=species), position = position_jitter(width = 0.1, height=0.05), col="blue") +
    geom_point(data=pre_test_set, aes(y=(inhibRich), x=species), position=position_jitter(width = 0.1, height=0.05), col="red")

# Get standard deviation between toad individuals and samples
toadID_sigma <- sd(samps_glmer_inhibRich$b[,ncol(samps_glmer_inhibRich$b)])

# Get estimates for control toads
inhib_species_exp <- fixef(glmer_inhibRich)  %>%
    as_tibble() %>%
    mutate(species=c("Anbo", "Anma","Lica","Lipi","Osse"))
inhib_species_exp$value[2:5] <- unlist(c(inhib_species_exp[1,1]+inhib_species_exp[2,1]
                                         , inhib_species_exp[1,1]+inhib_species_exp[3,1]
                                         , inhib_species_exp[1,1]+inhib_species_exp[4,1]
                                         , inhib_species_exp[1,1]+inhib_species_exp[5,1]))

con_toad_est_inhib <- ranef(glmer_inhibRich)$toadID %>%
    rename(sp_est="(Intercept)") %>%
    mutate(toadID=rownames(ranef(glmer_inhibRich)$toadID)) %>%
    separate(toadID, into=c("species","num"), sep="_",remove=FALSE) %>%
    dplyr::select(-num) %>%
    left_join(inhib_species_exp, by = "species") %>%
    mutate(est_inhib=sp_est+value)

# Now, we can calculate the probability that the "test" dataset values come from this distribution
# List of individuals
treat_indiv <- unique(mf_treat_without_init_infect$toadID)
# List of each species
species_list <- levels(factor(mf_con_without_init_infect$species))

exp_distr <- as.data.frame(matrix(ncol=length(species_list), nrow=4000, dimnames = list(1:4000, species_list)))
for ( num_sp in 1:length(species_list)) {
    exp_distr[,num_sp] <- rnorm(length(new_samps_beta[,num_sp]), mean=new_samps_beta[,num_sp], sd=toadID_sigma)
}

# Loop through and calculate probability of having diversity at that level
pre_exp_indiv <- data.frame(toadID=treat_indiv, exp_inhibRich=rep(NA, length(treat_indiv)), p_inhibRich=rep(NA, length(treat_indiv)), infect=rep(NA, length(treat_indiv)))
for ( i in treat_indiv ) {
    n_row <- match(i, treat_indiv)
    sp <- unlist(strsplit(i,"_"))
    num_sp <- match(sp[1], levels(factor(mf_con_without_init_infect$species)))
    temp_rich <- mf_treat_without_init_infect %>%
        filter(toadID==i, time <=5 ) %>%
        dplyr::select(inhibRich) %>%
        pull()
    if ( length(temp_rich)>1 ) {
        exp_inhibRich <-  fitdistr((temp_rich), "Poisson")$estimate[1]
        
    } else {
        exp_inhibRich <- log(temp_rich)
    }
    
    p_inhibRich <- sum(exp_distr[,sp[1]]<exp_inhibRich)/length(exp_distr[,sp[1]])
    
    ### Did they get infected? 
    infect <- max(mf_treat_without_init_infect %>%
                      filter(toadID==i) %>%
                      dplyr::select(eBD_raw) %>%
                      pull()
    )
    
    pre_exp_indiv[n_row,c("exp_inhibRich","p_inhibRich","infect")] <- c(exp_inhibRich, p_inhibRich, infect)
    
}

con_exp_indiv_inhib <- data.frame(toadID=con_toad_est_inhib$toadID, exp_inhibRich=rep(NA, length(con_toad_est_inhib$toadID)), p_inhibRich=rep(NA, length(con_toad_est_inhib$toadID)))
for ( i in 1:nrow(con_toad_est_inhib) ) {
    s <- con_toad_est_inhib[i,"species"]
    
    exp_rich <- con_toad_est_inhib[i,"est_inhib"]
    p_rich <- sum(exp_distr[,s]<exp_rich)/length(exp_distr[,s])
    
    con_exp_indiv_inhib[i,c("exp_inhibRich","p_inhibRich")] <- c(exp_rich, p_rich)
}

# create species column
pre_exp_indiv <- pre_exp_indiv %>%
    separate(toadID, into=c("species","indiv"), remove = FALSE)

# Plot results 
gg_inhibRich_p <- ggplot(pre_exp_indiv, aes(x=p_inhibRich, y=log(infect+1))) +
    geom_smooth(aes(color=species),method=lm, se = FALSE) +
    geom_point(aes(color=species), cex=4) +
    geom_smooth(method=lm, se=FALSE, col="black")
gg_inhibRich_raw <- ggplot(pre_exp_indiv, aes(x=exp_inhibRich, y=log(infect+1))) +
    geom_smooth(aes(color=species),method=lm, se = FALSE) +
    geom_point(aes(color=species), cex=4) +
    geom_smooth(method=lm, se=FALSE, col="black")
grid.arrange(gg_inhibRich_p, gg_inhibRich_raw, nrow=1)

exp_distr %>%
    gather(key=species, value=loginhibRich) %>%
    ggplot(aes(x=species, y=loginhibRich)) +
    geom_violin() +
    geom_point(data=pre_exp_indiv, aes(x=species, y=exp_inhibRich, col=log(infect+1)), cex=4, position=position_jitter(height=0, width=0.1))

all_p <- pre_exp_indiv %>%
    dplyr::select(toadID, exp_inhibRich, p_inhibRich) %>%
    full_join(all_p, by="toadID")

##### Save all work #####
save(all_p, file="all_p.RData")

##### PART II: Now, how does infection itself affect the microbiome? #####
mf_all_noinfect <-rbind(mf_treat_without_init_infect, mf_con_without_init_infect) %>%
    filter(!(BD_infected=="y" & exposure=="Post")) 

##### SHANNON ####

if ( FALSE ) {
    lmer_shannon_all <- stan_lmer(shannon ~ -1 + species + (1|toadID), data=mf_all_noinfect
                              , prior = normal(0, 10, autoscale = TRUE)
                              , seed = 98374)
    save(lmer_shannon_all, file="lmer_shannon_all.RData")
} else {
    load("lmer_shannon_all.RData")
}
prior_summary(lmer_shannon_all)
# Look at distributions according to models
samps_lmer_shannon_all <- rstan::extract(lmer_shannon_all$stanfit)
post_test_set <- mf_treat_without_init_infect %>%
    filter(time>=6) 
## PLOT OF EXPECTED RCIHNESS FOR EACH SPECIES
samps_lmer_shannon_all$beta %>%
    as.data.frame() %>%
    rename(Anbo=V1, Anma=V2, Lica=V3, Lipi=V4, Osse=V5) %>%
    dplyr::select(Anbo,Anma,Lica,Lipi,Osse) %>%
    gather(key=species, value=shannon) %>%
    ggplot(mapping=aes(x=species, y=shannon))+
    geom_violin() +
    geom_point(data=mf_con_without_init_infect, aes(y=shannon, x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
    geom_point(data=post_test_set, aes(y=shannon, x=species), position=position_jitter(width = 0.1, height=0), col="red")

# Get standard deviation between toad individuals and samples
samp_sigma <- samps_lmer_shannon_all$aux
toad_intercept <- ranef(lmer_shannon_all)$toadID
# toadID_sigma <- sd(samps_lmer_shannon_all$b[,ncol(samps_lmer_shannon_all$b)])
samp_toad <- samps_lmer_shannon_all$b[,1:nrow(toad_intercept)]
colnames(samp_toad) <- rownames(toad_intercept)

# Now, we can calculate the probability that the "test" dataset values come from this distribution
# List of individuals
treat_indiv <- unique(mf_treat_without_init_infect$toadID)
# Make key of species for each individual?
species_key <- treat_indiv %>%
    as_tibble() %>%
    rename(toadID=value) %>%
    separate(toadID,into=c("species","indiv"), remove=FALSE)
species_order <- levels(as.factor(mf_all_noinfect$species))

exp_distr <- as.data.frame(matrix(ncol=length(treat_indiv), nrow=4000, dimnames = list(1:4000, treat_indiv)))
for ( num_indiv in 1:length(treat_indiv)) {
    indiv <- treat_indiv[num_indiv]
    sp <- pull(species_key[num_indiv,"species"])
    num_sp <- match(sp, species_order)
    exp_distr[,num_indiv] <- rnorm(4000, mean=rnorm(4000, mean=(samps_lmer_shannon_all$beta[,num_sp]+samp_toad[,indiv]), sd=samp_sigma))
}

pos_exp_indiv <- mf_treat_without_init_infect %>%
    filter(time>5) %>%
    dplyr::select(toadID, time, species, shannon, eBD_log) %>%
    mutate(p_shan=NA)
for ( r in 1:nrow(pos_exp_indiv)) {
    pos_exp_indiv[r,"p_shan"] <- sum(exp_distr[,pos_exp_indiv$toadID[r]]<pos_exp_indiv[r,"shannon"])/4000
}

gg_shan_pos_p <- pos_exp_indiv %>%
    ggplot(aes(x=p_shan, y=eBD_log)) +
    geom_point(aes(col=species), cex=3, position=position_jitter(height=0.15)) +
    geom_smooth(aes(col=species), method="lm",se=FALSE) +
    geom_smooth(method="lm", se=FALSE, col="black")
gg_shan_pos_raw <- pos_exp_indiv %>%
    ggplot(aes(x=shannon, y=eBD_log)) +
    geom_point(aes(col=species), cex=3, position=position_jitter(height=0.15)) +
    geom_smooth(aes(col=species), method="lm",se=FALSE) +
    geom_smooth(method="lm", se=FALSE, col="black")
grid.arrange(gg_shan_pos_p, gg_shan_pos_raw, nrow=1)

exp_distr %>%
    gather(key=toadID, value=shannon) %>%
    ggplot(aes(x=toadID, y=shannon)) +
    geom_violin() +
    geom_point(data=pos_exp_indiv, aes(x=toadID, y=shannon, col=eBD_log),cex=4, position=position_jitter(height=0, width=0.1))

all_p_infected <- pos_exp_indiv


#### RICHNESS (observed otus) #####

if (FALSE) {
    lmer_rich_all <- stan_lmer(logRich ~ -1 + species + (1|toadID), data=mf_all_noinfect
                           , prior = normal(0, 10, autoscale = TRUE)
                           , seed = 98374)
    save(lmer_rich_all, file="lmer_rich_all.RData")
} else {
    load(file="lmer_rich_all.RData")
}
prior_summary(lmer_rich_all)

# Look at distributions according to models
samps_lmer_rich_all <- rstan::extract(lmer_rich_all$stanfit)
pos_test_set <- mf_treat_without_init_infect %>%
    filter(time>=6) 
# Plot of expected by species, and real samples
samps_lmer_rich_all$beta %>%
    as.data.frame() %>%
    rename(Anbo=V1, Anma=V2, Lica=V3, Lipi=V4, Osse=V5) %>%
    dplyr::select(Anbo,Anma,Lica,Lipi,Osse) %>%
    gather(key=species, value=logRich) %>%
    ggplot(mapping=aes(x=species, y=logRich))+
    geom_violin() +
    geom_point(data=mf_all_noinfect, aes(y=logRich, x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
    geom_point(data=pos_test_set, aes(y=logRich, x=species), position=position_jitter(width = 0.1, height=0), col="red")

# Get standard deviation between toad individuals and samples
samp_sigma <- samps_lmer_rich_all$aux
toad_intercept <- ranef(lmer_rich_all)$toadID
samp_toad <- samps_lmer_rich_all$b[,1:nrow(toad_intercept)]
colnames(samp_toad) <- rownames(toad_intercept)

# Now, we can calculate the probability that the "test" dataset values come from this distribution
# List of individuals
treat_indiv <- unique(mf_treat_without_init_infect$toadID)
# Make key of species for each individual?
species_key <- treat_indiv %>%
    as_tibble() %>%
    rename(toadID=value) %>%
    separate(toadID,into=c("species","indiv"), remove=FALSE)
species_order <- levels(as.factor(mf_all_noinfect$species))

exp_distr <- as.data.frame(matrix(ncol=length(treat_indiv), nrow=4000, dimnames = list(1:4000, treat_indiv)))
for ( num_indiv in 1:length(treat_indiv)) {
    indiv <- treat_indiv[num_indiv]
    sp <- pull(species_key[num_indiv,"species"])
    num_sp <- match(sp, species_order)
    exp_distr[,num_indiv] <- rnorm(4000, mean=rnorm(4000, mean=(samps_lmer_rich_all$beta[,num_sp]+samp_toad[,indiv]), sd=samp_sigma))
}

pos_exp_indiv <- mf_treat_without_init_infect %>%
    filter(time>5) %>%
    dplyr::select(toadID, time, species, logRich, eBD_log) %>%
    mutate(p_rich=NA)
for ( r in 1:nrow(pos_exp_indiv)) {
    pos_exp_indiv[r,"p_rich"] <- sum(exp_distr[,pos_exp_indiv$toadID[r]]<pos_exp_indiv[r,"logRich"])/4000
}

gg_logRich_pos_p <- pos_exp_indiv %>%
    ggplot(aes(x=p_rich, y=eBD_log)) +
    geom_point(aes(col=species), cex=3, position=position_jitter(height=0.15))+
    geom_smooth(aes(col=species), method="lm",se=FALSE) +
    geom_smooth(method="lm", se=FALSE, col="black")
gg_logRich_pos_raw <- pos_exp_indiv %>%
    ggplot(aes(x=logRich, y=eBD_log)) +
    geom_point(aes(col=species), cex=3, position=position_jitter(height=0.15))+
    geom_smooth(aes(col=species), method="lm",se=FALSE) +
    geom_smooth(method="lm", se=FALSE, col="black")
grid.arrange(gg_logRich_pos_p, gg_logRich_pos_raw, nrow=1)

exp_distr %>%
    gather(key=toadID, value=logRich) %>%
    ggplot(aes(x=toadID, y=logRich)) +
    geom_violin() +
    geom_point(data=pos_exp_indiv, aes(x=toadID, y=logRich, col=eBD_log),cex=4, position=position_jitter(height=0, width=0.1))

all_p_infected <- pos_exp_indiv %>%
    dplyr::select(toadID, time, logRich, p_rich) %>%
    full_join(all_p_infected, by=c("toadID","time"))


##### BETA DIVERSITY ####
if ( FALSE) {
    glmer_BC_all <- stan_glmer(distance_bray_curtis ~ -1 + species + (1|toadID)
                           , data=mf_all_noinfect
                           , family =mgcv::betar
                           , prior_intercept = normal(location = 0.5,scale = 0.1, autoscale = TRUE)
                           , prior = normal(location=0.5, scale=0.1, autoscale=TRUE)
                           , seed= 623445
    )
    save(glmer_BC_all, file="glmer_BC_all.RData")
} else {
    load("glmer_BC_all.RData")
}
prior_summary(glmer_BC_all)

# rbeta has a strange parameterization using a nd b so need to convert mu and phi to this.
a <- function(mu,phi){
    mu*phi
}
b <- function(mu,phi) {
    phi-mu*phi
}
mu <- function(a,phi) {
    a/phi
}
inv_logit <- function(x) {
    exp(x)/(exp(x)+1)
}
logit <- function(p) {
    log(p/(1-p))
}

# Look at distributions according to models
samps_glmer_BC_all<- rstan::extract(glmer_BC_all$stanfit)
pos_test_set <- mf_treat_without_init_infect %>%
    filter(time>5) 
samps_glmer_BC_all$beta %>%
    as.data.frame() %>%
    rename(Anbo=V1, Anma=V2, Lica=V3, Lipi=V4, Osse=V5) %>%
    mutate(Anbo=inv_logit(mu(Anbo, samps_glmer_BC_all$aux))
           ,Anma=inv_logit(mu(Anma, samps_glmer_BC_all$aux))
           ,Lica=inv_logit(mu(Lica, samps_glmer_BC_all$aux))
           ,Lipi=inv_logit(mu(Lipi, samps_glmer_BC_all$aux))
           ,Osse=inv_logit(mu(Osse, samps_glmer_BC_all$aux))) %>%
    dplyr::select(Anbo,Anma,Osse,Lica,Lipi) %>%
    gather(key=species, value=distance_bray_curtis) %>%
    ggplot(mapping=aes(x=species, y=distance_bray_curtis))+
    geom_violin() +
    geom_point(data=mf_all_noinfect, aes(y=distance_bray_curtis, x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
    geom_point(data=pos_test_set, aes(y=distance_bray_curtis, x=species), position=position_jitter(width = 0.1, height=0), col="red")

toad_intercept <- ranef(glmer_BC_all)$toadID
samp_toad <- samps_glmer_BC_all$b[,1:nrow(toad_intercept)]
colnames(samp_toad) <- rownames(toad_intercept)
# toadID_sigma <- sd(samps_glmer_BC_all$b[,ncol(samps_glmer_BC_all$b)])
phi <- samps_glmer_BC_all$aux

# Now, we can calculate the probability that the "test" dataset values come from this distribution
# List of individuals (different here bc some individuals don't have bray-curtis values)
treat_indiv <- colnames(samp_toad)

# Make key of species for each individual
species_key <- treat_indiv %>%
    as_tibble() %>%
    rename(toadID=value) %>%
    separate(toadID,into=c("species","indiv"), remove=FALSE)
species_order <- levels(as.factor(mf_all_noinfect$species))

exp_distr <- as.data.frame(matrix(ncol=length(treat_indiv), nrow=4000, dimnames = list(1:4000, treat_indiv)))
for ( num_indiv in 1:length(treat_indiv)) {
    indiv <- treat_indiv[num_indiv]
    sp <- pull(species_key[num_indiv,"species"])
    num_sp <- match(sp, species_order)
    
    mu <- inv_logit(samps_glmer_BC_all$beta[,num_sp]+samp_toad[,indiv])
    exp_distr[,num_indiv] <- rbeta(4000
                                ,shape1=a(mu,samps_glmer_BC_all$aux)
                                ,shape2=b(mu,samps_glmer_BC_all$aux))
}

pos_exp_indiv <- mf_treat_without_init_infect %>%
    filter(time>5, !is.na(distance_bray_curtis), toadID %in% colnames(samp_toad)) %>%
    dplyr::select(toadID, time, species, distance_bray_curtis, eBD_log) %>%
    mutate(p_BC=NA)
for ( r in 1:nrow(pos_exp_indiv)) {
    pos_exp_indiv[r,"p_BC"] <- sum(exp_distr[,pos_exp_indiv$toadID[r]]<pos_exp_indiv[r,"distance_bray_curtis"])/4000
}

gg_beta_pos_p <- pos_exp_indiv %>%
    ggplot(aes(x=p_BC, y=eBD_log)) +
    geom_point(aes(col=species), cex=3, position=position_jitter(height=0.15))+
    geom_smooth(aes(col=species), method="lm",se=FALSE) +
    geom_smooth(method="lm", se=FALSE, col="black")
gg_beta_pos_raw <- pos_exp_indiv %>%
    ggplot(aes(x=distance_bray_curtis, y=eBD_log)) +
    geom_point(aes(col=species), cex=3, position=position_jitter(height=0.15))+
    geom_smooth(aes(col=species), method="lm",se=FALSE) +
    geom_smooth(method="lm", se=FALSE, col="black")
grid.arrange(gg_beta_pos_p, gg_beta_pos_raw, nrow=1)

exp_distr %>%
    gather(key=toadID, value=distance_bray_curtis) %>%
    ggplot(aes(x=toadID, y=distance_bray_curtis)) +
    geom_violin() +
    geom_point(data=pos_exp_indiv, aes(x=toadID, y=distance_bray_curtis, col=eBD_log),cex=4, position=position_jitter(height=0, width=0.1))

all_p_infected <- pos_exp_indiv %>%
    dplyr::select(toadID, time, distance_bray_curtis, p_BC) %>%
    full_join(all_p_infected, by=c("toadID","time"))


##### PERCENT INHIBITORY #####

inhibBin <- cbind(inhibCount=mf_all_noinfect$inhibCounts, noninhibCount=mf_all_noinfect$n)

if ( FALSE) {
    glmer_percInhib_all <- stan_glmer(inhibBin ~ -1 + species + (1|toadID) +(1|SampleID)
                                  , data=mf_all_noinfect
                                  , family = 'binomial'
                                  , seed= 9837423)
    save(glmer_percInhib_all, file="glmer_percInhib_all.RData")
} else {
    load("glmer_percInhib_all.RData")
}
prior_summary(glmer_percInhib_all)

# Look at distributions according to models
samps_glmer_percInhib_all<- rstan::extract(glmer_percInhib_all$stanfit)
pos_test_set <- mf_treat_without_init_infect %>%
    filter(time>5) 
samps_glmer_percInhib_all$beta %>%
    as.data.frame() %>%
    rename(Anbo=V1, Anma=V2, Lica=V3, Lipi=V4, Osse=V5) %>%
    mutate(Anbo=inv_logit(Anbo)
           ,Anma=inv_logit(Anma)
           ,Lica=inv_logit(Lica)
           ,Lipi=inv_logit(Lipi)
           ,Osse=inv_logit(Osse)) %>%
    dplyr::select(Anbo,Anma,Osse,Lica,Lipi) %>%
    gather(key=species, value=percInhib) %>%
    ggplot(mapping=aes(x=species, y=percInhib))+
    geom_violin() +
    geom_point(data=mf_all_noinfect, aes(y=percInhib, x=species), position = position_jitter(width = 0.1, height=0), col="blue") +
    geom_point(data=pos_test_set, aes(y=percInhib, x=species), position=position_jitter(width = 0.1, height=0), col="red")


# indiv_mu <- ranef(glmer_percInhib_all)$toadID
# sp_mu <- fixef(glmer_BC)
toad_intercept <- ranef(glmer_percInhib_all)$toadID
sample_intercept <- ranef(glmer_percInhib_all)$SampleID

samp_toad <- samps_glmer_percInhib_all$b[,(nrow(sample_intercept)+2):(nrow(sample_intercept)+1+nrow(toad_intercept))]
colnames(samp_toad) <- rownames(toad_intercept)

sample_sigma <- sd(samps_glmer_percInhib_all$b[,nrow(sample_intercept)+1])

# Now, we can calculate the probability that the "test" dataset values come from this distribution
# List of individuals (different here bc some individuals don't have bray-curtis values)
treat_indiv <- colnames(samp_toad)

# Make key of species for each individual
species_key <- treat_indiv %>%
    as_tibble() %>%
    rename(toadID=value) %>%
    separate(toadID,into=c("species","indiv"), remove=FALSE)
species_order <- levels(as.factor(mf_all_noinfect$species))

mu_exp_distr <- as.data.frame(matrix(ncol=length(treat_indiv), nrow=4000, dimnames = list(1:4000, treat_indiv)))
for ( num_indiv in 1:length(treat_indiv)) {
    indiv <- treat_indiv[num_indiv]
    sp <- pull(species_key[num_indiv,"species"])
    num_sp <- match(sp, species_order)
    
    mu_exp_distr[,num_indiv] <- inv_logit(rnorm(4000, mean=samps_glmer_percInhib_all$beta[,num_sp] + samp_toad[,indiv], sd=sample_sigma))
    
}

pos_exp_indiv <- mf_treat_without_init_infect %>%
    filter(time>5) %>%
    dplyr::select(toadID, time, species, percInhib, inhibCounts,n, eBD_log) %>%
    mutate(p_percInhib=NA)
for ( r in 1:nrow(pos_exp_indiv)) {
    temp_samp <- rbinom(n=4000, size=pos_exp_indiv[r,"n"], prob = mu_exp_distr[,pos_exp_indiv$toadID[r]])
    
    pos_exp_indiv[r,"p_percInhib"] <- sum(temp_samp < pos_exp_indiv[r,"inhibCounts"])/4000
}



gg_percInhib_pos_p <- pos_exp_indiv %>%
    ggplot(aes(x=p_percInhib, y=eBD_log)) +
    geom_point(aes(col=species), cex=3, position=position_jitter(height=0.15))+
    geom_smooth(aes(col=species), method="lm",se=FALSE) +
    geom_smooth(method="lm", se=FALSE, col="black")
gg_percInhib_pos_raw <- pos_exp_indiv %>%
    ggplot(aes(x=percInhib, y=eBD_log)) +
    geom_point(aes(col=species), cex=3, position=position_jitter(height=0.15))+
    geom_smooth(aes(col=species), method="lm",se=FALSE) +
    geom_smooth(method="lm", se=FALSE, col="black")
grid.arrange(gg_percInhib_pos_p, gg_percInhib_pos_raw, nrow=1)

mu_exp_distr %>%
    gather(key=toadID, value=percInhib) %>%
    ggplot(aes(x=toadID, y=percInhib)) +
    geom_violin() +
    geom_point(data=pos_exp_indiv, aes(x=toadID, y=percInhib, col=eBD_log),cex=4, position=position_jitter(height=0, width=0.1))

all_p_infected <- pos_exp_indiv %>%
    dplyr::select(toadID, time, percInhib, p_percInhib) %>%
    full_join(all_p_infected, by=c("toadID","time"))



######## INHIB RICHNESS ############


if (FALSE) {
    glmer_inhibRich_all <- stan_glmer(inhibRich ~ -1 + species + (1|toadID) + (1|SampleID), data=mf_all_noinfect
                                  , prior = normal(0, 10, autoscale = TRUE)
                                  , family= poisson(link="log")
                                  , seed = 5423409)
    save(glmer_inhibRich_all, file="glmer_inhibRich_all.RData")
} else {
    load(file="glmer_inhibRich_all.RData")
}
prior_summary(glmer_inhibRich_all)

# Look at distributions according to models
samps_glmer_inhibRich_all <- rstan::extract(glmer_inhibRich_all$stanfit)
pos_test_set <- mf_treat_without_init_infect %>%
    filter(time>=6) 
samps_glmer_inhibRich_all$beta %>%
    as.data.frame() %>%
    rename(Anbo=V1, Anma=V2, Lica=V3, Lipi=V4, Osse=V5) %>%
    dplyr::select(Anbo,Anma,Lica,Lipi,Osse) %>%
    gather(key=species, value=inhibRich) %>%
    ggplot(mapping=aes(x=species, y=inhibRich))+
    geom_violin() +
    geom_point(data=mf_all_noinfect, aes(y=log(inhibRich), x=species), position = position_jitter(width = 0.1, height=0.05), col="blue") +
    geom_point(data=pos_test_set, aes(y=log(inhibRich), x=species), position=position_jitter(width = 0.1, height=0.05), col="red")

# Get standard deviation between toad individuals and samples
samp_intercept <- ranef(glmer_inhibRich_all)$SampleID
toad_intercept <- ranef(glmer_inhibRich_all)$toadID

sample_sigma <- sd(samps_glmer_inhibRich_all$b[,nrow(samp_intercept)+1])
toadID_sigma <- sd(samps_glmer_inhibRich_all$b[,ncol(samps_glmer_inhibRich_all$b)])

samp_toad <- samps_glmer_inhibRich_all$b[,(nrow(samp_intercept)+2):(nrow(samp_intercept)+1+nrow(toad_intercept))]
colnames(samp_toad) <- rownames(toad_intercept)

# Now, we can calculate the probability that the "test" dataset values come from this distribution
# List of individuals (different here bc some individuals don't have bray-curtis values)
treat_indiv <- colnames(samp_toad)

# Make key of species for each individual
species_key <- treat_indiv %>%
    as_tibble() %>%
    rename(toadID=value) %>%
    separate(toadID,into=c("species","indiv"), remove=FALSE)
species_order <- levels(as.factor(mf_all_noinfect$species))

exp_distr <- as.data.frame(matrix(ncol=length(treat_indiv), nrow=4000, dimnames = list(1:4000, treat_indiv)))
for ( num_indiv in 1:length(treat_indiv)) {
    indiv <- treat_indiv[num_indiv]
    sp <- pull(species_key[num_indiv,"species"])
    num_sp <- match(sp, species_order)
    
    exp_distr[,num_indiv] <- exp(rnorm(4000, mean=samps_glmer_inhibRich_all$beta[,num_sp]+samp_toad[,indiv], sd=sample_sigma))
    
}

pos_exp_indiv <- mf_treat_without_init_infect %>%
    filter(time>5) %>%
    dplyr::select(toadID, time, species, inhibRich, eBD_log) %>%
    mutate(p_inhibRich=NA)
for ( r in 1:nrow(pos_exp_indiv)) {
    temp_samp <- rpois(n=4000, lambda = exp_distr[,pos_exp_indiv[r,"toadID"]])
    
    pos_exp_indiv[r,"p_inhibRich"] <- sum(temp_samp < pos_exp_indiv[r,"inhibRich"])/4000
}



gg_inhibRich_pos_p <- pos_exp_indiv %>%
    ggplot(aes(x=p_inhibRich, y=eBD_log)) +
    geom_point(aes(col=species), cex=3, position=position_jitter(height=0.15))+
    geom_smooth(aes(col=species), method="lm",se=FALSE) +
    geom_smooth(method="lm", se=FALSE, col="black")
gg_inhibRich_pos_raw <- pos_exp_indiv %>%
    ggplot(aes(x=inhibRich, y=eBD_log)) +
    geom_point(aes(col=species), cex=3, position=position_jitter(height=0.15))+
    geom_smooth(aes(col=species), method="lm",se=FALSE) +
    geom_smooth(method="lm", se=FALSE, col="black")
grid.arrange(gg_inhibRich_pos_p, gg_inhibRich_pos_raw, nrow=1)

exp_distr %>%
    gather(key=toadID, value=inhibRich) %>%
    ggplot(aes(x=toadID, y=inhibRich)) +
    geom_violin() +
    geom_point(data=pos_exp_indiv, aes(x=toadID, y=inhibRich, col=eBD_log),cex=4, position=position_jitter(height=0, width=0.1))


all_p_infected <- pos_exp_indiv %>%
    dplyr::select(toadID, time, inhibRich, p_inhibRich) %>%
    full_join(all_p_infected, by=c("toadID","time"))

#### save work #####
save(all_p_infected, file="all_p_infected.RData")


#### Are inhibRich and overall rich correlated? ####

# Extract estimates for each control toad and see if inhibitory richness and overall richness is correlated
con_exp_indiv <- con_exp_indiv_rich %>%
    full_join(con_exp_indiv_inhib, by="toadID")
save(con_exp_indiv, file="con_exp_indiv.RData")

all_p %>%
    dplyr::select(toadID, exp_rich, p_rich, exp_inhibRich, p_inhibRich) %>%
    # full_join(con_exp_indiv, by = "toadID") %>%
    rbind(con_exp_indiv) %>%
    dplyr::select(p_inhibRich, p_rich, toadID) %>%
    separate(toadID, into=c("species","n"), remove=FALSE) %>%
    ggplot(aes(x=p_rich, y=p_inhibRich)) +
    geom_point(aes(col=species), cex=3)

# anova(lm(p_inhib ~ p_rich, data=con_exp_indiv))

#### Inhibitory Bacteria ####

#Exact same?
otu.inhibOnly.con$OTUID == otu.inhibOnly.treat$OTUID
# Good

# Look at "must abundant" otus
combined_otu_inhibOnly <- cbind(otu.inhibOnly.con, otu.inhibOnly.treat)[,-(ncol(otu.inhibOnly.con)+ncol(otu.inhibOnly.treat))]
OTUs_temp <- combined_otu_inhibOnly %>%
    dplyr::select(OTUID) %>%
    mutate(taxa = (inhib.tb[match(OTUID, inhib.tb$seq),"taxa"])$taxa ) 
    
OTUs_temp[order(combined_otu_inhibOnly %>%
                    dplyr::select(-OTUID) %>%
                    rowSums(.), decreasing = TRUE),]

# chosen colors
brewer.pal.info
set.seed(8984)
set_col <- unique(c(brewer.pal(12,"Set3"), brewer.pal(8,"Spectral"), brewer.pal(8,"Dark2"), brewer.pal(8,"Accent")))
set_col <- set_col[c(1,35,34,3,33,4,5,30,7,29,8,28,9,27,26,10,11,25,12,24,13,23,14,22,15,21,16,20,17,19,18)]

#### CONTROL ####
# Going to leave zeros in for now, bc only 3-5 in each control and con
temp <- otu.inhibOnly.con %>%
    dplyr::select(-c(OTUID)) %>%
    unlist()
# THE ABOVE LISTS BY STACKING COL ON TOP OF EACH OTHER (not by stacking rows beside each other)

otu_long_con <- cbind(rep(colnames(otu.inhibOnly.con)[-ncol(otu.inhibOnly.con)], each=nrow(otu.inhibOnly.con))
                        ,temp
                        ,rep(otu.inhibOnly.con$OTUID,times=ncol(otu.inhibOnly.con)-1)) %>%
    as_tibble() %>%
    rename(SampleID=V1, reads=temp, OTUID=V2) %>%
    mutate(reads=as.numeric(reads))
mf_con_with_inhibOTUs <- mf_con_without_init_infect %>%
    dplyr::select(SampleID, species, time, toadID, eBD_log, PABD, prepost, BD_infected) %>%
    left_join(otu_long_con) %>%
    mutate(taxa = (inhib.tb[match(OTUID, inhib.tb$seq),"taxa"])$taxa ) %>%
    separate(taxa, into = c("K","P","C","O","F","G"), sep = ";", remove = FALSE, fill="left")

mf_con_with_inhibOTUs %>%
    group_by(toadID, time, species, G) %>%
    summarise(AverageProportion = mean(reads)) %>%
    ggplot(aes(x=time, y=AverageProportion)) +
    geom_bar(aes(fill=G), stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.key.size =  unit(0.2, "cm") ) +
    facet_wrap(~species, nrow=2) +
    scale_fill_manual(values=set_col)


#### TREATMENT ####
# Going to leave zeros in for now, bc only 3-5 in each control and treatment
temp <- otu.inhibOnly.treat %>%
    dplyr::select(-OTUID) %>%
    unlist()
# THE ABOVE LISTS BY STACKING COL ON TOP OF EACH OTHER (not by stacking rows beside each other)

otu_long_treat <- cbind(rep(colnames(otu.inhibOnly.treat)[-ncol(otu.inhibOnly.treat)], each=nrow(otu.inhibOnly.treat))
                        ,temp
                        ,rep(otu.inhibOnly.treat$OTUID,times=ncol(otu.inhibOnly.treat)-1)) %>%
    as_tibble() %>%
    rename(SampleID=V1, reads=temp, OTUID=V2) %>%
    mutate(reads=as.numeric(reads))
mf_treat_with_inhibOTUs <- mf_treat_without_init_infect %>%
    dplyr::select(SampleID, species, time, toadID, eBD_log, PABD, prepost, BD_infected) %>%
    left_join(otu_long_treat)%>%
    mutate(taxa = (inhib.tb[match(OTUID, inhib.tb$seq),"taxa"])$taxa ) %>%
    separate(taxa, into = c("K","P","C","O","F","G"), sep = ";", remove = FALSE, fill="left")

mf_treat_with_inhibOTUs %>%
    group_by(toadID, time, species, G) %>%
    summarise(AverageProportion = mean(reads)) %>%
    ggplot(aes(x=time, y=AverageProportion)) +
    geom_bar(aes(fill=G), stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.key.size =  unit(0.2, "cm") ) +
    facet_wrap(~species, nrow=2) +
    scale_fill_manual(values=set_col) +
    geom_vline(aes(xintercept=5.5), col="grey", lty=2)

g_legend <- function(a.gplot){ # from stack overflow
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    return(legend)} 

### Goal: plot abundance of inhibitory bacteria over time for each individual of each species
anbo_con <- mf_con_with_inhibOTUs %>%
    filter(species=="Anbo") %>%
    group_by(species, toadID, G, time) %>%
    summarise(reads=mean(reads)) %>%
    ggplot(aes(x=time, y=reads)) +
    geom_line(aes(col=G), stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~toadID, nrow=1) +
    scale_fill_manual(values=set_col) + 
    scale_color_manual(values=set_col) + 
    ylab("Anbo") +
    xlab("")
anma_con <- mf_con_with_inhibOTUs %>%
    filter(species=="Anma") %>%
    group_by(species, toadID, G, time) %>%
    summarise(reads=mean(reads)) %>%
    ggplot(aes(x=time, y=reads)) +
    geom_line(aes(col=G), stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~toadID, nrow=1) +
    scale_fill_manual(values=set_col) + 
    scale_color_manual(values=set_col) + 
    ylab("Anma")+
    xlab("")
lica_con <- mf_con_with_inhibOTUs %>%
    filter(species=="Lica") %>%
    group_by(species, toadID, G, time) %>%
    summarise(reads=mean(reads)) %>%
    ggplot(aes(x=time, y=reads)) +
    geom_line(aes(col=G), stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~toadID, nrow=1) +
    scale_fill_manual(values=set_col) + 
    scale_color_manual(values=set_col) + 
    ylab("Lica")+
    xlab("")
lipi_con <- mf_con_with_inhibOTUs %>%
    filter(species=="Lipi") %>%
    group_by(species, toadID, G, time) %>%
    summarise(reads=mean(reads)) %>%
    ggplot(aes(x=time, y=reads)) +
    geom_line(aes(col=G), stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~toadID, nrow=1) +
    scale_fill_manual(values=set_col) + 
    scale_color_manual(values=set_col) + 
    ylab("Lipi")+
    xlab("")
osse_con <- mf_con_with_inhibOTUs %>%
    filter(species=="Osse") %>%
    group_by(species, toadID, G, time) %>%
    summarise(reads=mean(reads)) %>%
    ggplot(aes(x=time, y=reads)) +
    geom_line(aes(col=G), stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~toadID, nrow=1) +
    scale_fill_manual(values=set_col) + 
    scale_color_manual(values=set_col) + 
    ylab("Osse")+
    xlab("")
lay_con <- rbind(c(1,1,1,1)
                 ,c(2,2,2,6)
                 ,c(3,3,3,3)
                 , c(4,4,7,7)
                 , c(5,5,5,5))
grid.arrange(anbo_con, anma_con, lica_con, lipi_con, osse_con
             , layout_matrix=lay_con
             , left = textGrob("Average proportion of reads", rot = 90, vjust = 1)
             , bottom = textGrob("Time", vjust=1))


# Now for treatment 
### Goal: plot abundance of inhibitory bacteria over time for each individual of each species
anbo_treat <- mf_treat_with_inhibOTUs %>%
    filter(species=="Anbo") %>%
    group_by(species, toadID, G, time) %>%
    summarise(reads=mean(reads)) %>%
    ggplot(aes(x=time, y=reads)) +
    geom_line(aes(col=G), stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~toadID, nrow=1) +
    scale_fill_manual(values=set_col) + 
    scale_color_manual(values=set_col) + 
    ylab("Anbo") +
    xlab("")
anma_treat <- mf_treat_with_inhibOTUs %>%
    filter(species=="Anma") %>%
    group_by(species, toadID, G, time) %>%
    summarise(reads=mean(reads)) %>%
    ggplot(aes(x=time, y=reads)) +
    geom_line(aes(col=G), stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~toadID, nrow=1) +
    scale_fill_manual(values=set_col) + 
    scale_color_manual(values=set_col) + 
    ylab("Anma") +
    xlab("")
lica_treat <- mf_treat_with_inhibOTUs %>%
    filter(species=="Lica") %>%
    group_by(species, toadID, G, time) %>%
    summarise(reads=mean(reads)) %>%
    ggplot(aes(x=time, y=reads)) +
    geom_line(aes(col=G), stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~toadID, nrow=1) +
    scale_fill_manual(values=set_col) + 
    scale_color_manual(values=set_col) + 
    ylab("Lica") +
    xlab("")
lipi_treat <- mf_treat_with_inhibOTUs %>%
    filter(species=="Lipi") %>%
    group_by(species, toadID, G, time) %>%
    summarise(reads=mean(reads)) %>%
    ggplot(aes(x=time, y=reads)) +
    geom_line(aes(col=G), stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~toadID, nrow=1) +
    scale_fill_manual(values=set_col) + 
    scale_color_manual(values=set_col) + 
    ylab("Lipi") +
    xlab("")
osse_treat <- mf_treat_with_inhibOTUs %>%
    filter(species=="Osse") %>%
    group_by(species, toadID, G, time) %>%
    summarise(reads=mean(reads)) %>%
    ggplot(aes(x=time, y=reads)) +
    geom_line(aes(col=G), stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~toadID, nrow=1) +
    scale_fill_manual(values=set_col) + 
    scale_color_manual(values=set_col) + 
    ylab("Osse") +
    xlab("")
lay_treat <- rbind(c(1,1,1,1,1,1)
                 ,c(2,2,2,6,6,6)
                 ,c(3,3,3,3,3,3)
                 , c(4,7,7,7,7,7)
                 , c(5,5,5,5,5,5))
grid.arrange(anbo_treat, anma_treat, lica_treat, lipi_treat, osse_treat
             , layout_matrix=lay_treat
             , left = textGrob("Average proportion of reads", rot = 90, vjust = 1)
             , bottom = textGrob("Time", vjust=1))


### Summraizing before and after exposure ####

# CONTROL 
mf_con_statdiff <- mf_con_with_inhibOTUs %>%
    group_by(toadID, species, prepost, G) %>%
    summarize(meanReads=mean(reads), meanBd = mean(eBD_log)) %>%
    dplyr::select(-meanBd) %>%
    spread(key=prepost, value=meanReads) %>%
    mutate(fc=log((Pos-Pre)/Pre +1), significant = NA, p=NA) %>%
    mutate(temp= ifelse(is.finite(fc), fc, ifelse(is.infinite(fc), 10, 0)))

#+ warning=FALSE, message=FALSE
allInhibTaxa <- unique(mf_con_with_inhibOTUs$G)
allIndiv <- unique(mf_con_with_inhibOTUs$toadID)
for ( inhib in allInhibTaxa) {
    # inhib <- allInhibTaxa[1]
    # inhib <- "Arthrobacter"
    for ( indiv in allIndiv) {
        if (exists("stat_temp")) {
            remove(stat_temp)
        }
        # indiv <- allIndiv[1]
        # indiv <- "Anbo_2"
        mf_temp <- mf_con_with_inhibOTUs %>%
            filter(G==inhib, toadID==indiv) %>%
            group_by(toadID, G, prepost, time) %>%
            summarise(reads=mean(reads)) %>%
            spread(key=prepost, value=reads)
        # stat_temp <- anova(lm(reads ~ prepost, data=mf_temp))
        #+ warning=FALSE, message=FALSE
        stat_temp <- wilcox.test(mf_temp$Pre, mf_temp$Pos)
        mf_con_statdiff[which(mf_con_statdiff$toadID==indiv & mf_con_statdiff$G == inhib),"p"] <- (stat_temp$p.value)
        mf_con_statdiff[which(mf_con_statdiff$toadID==indiv & mf_con_statdiff$G == inhib),"significant"] <- (stat_temp$p.value<0.05)
        
    }
}

mf_con_statdiff %>%
    filter(!is.na(significant)) %>%
    ggplot(aes(x=G)) +
    geom_point(aes(y=temp, col=species, pch=significant), cex=2, position=position_jitter(width=0, height=0.25)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_shape_manual(values=c(21,19)) +
    geom_hline(aes(yintercept=0))


# TREAT

mf_treat_statdiff <- mf_treat_with_inhibOTUs %>%
    group_by(toadID, species, prepost, G) %>%
    summarize(meanReads=mean(reads), meanBd = mean(eBD_log)) %>%
    dplyr::select(-meanBd) %>%
    spread(key=prepost, value=meanReads) %>%
    # mutate(diff=Pos-Pre, significant = NA, p=NA) %>%
    mutate(fc=log((Pos-Pre)/Pre +1), significant = NA, p=NA) %>%
    mutate(temp= ifelse(is.finite(fc), fc, ifelse(is.infinite(fc), 10, 0)))

#+ warning=FALSE, message=FALSE
allInhibTaxa <- unique(mf_treat_with_inhibOTUs$G)
allIndiv <- unique(mf_treat_with_inhibOTUs$toadID)
for ( inhib in allInhibTaxa) {
    # inhib <- allInhibTaxa[1]
    # inhib <- "Arthrobacter"
    for ( indiv in allIndiv) {
        if (exists("stat_temp")) {
            remove(stat_temp)
        }
        # indiv <- allIndiv[1]
        # indiv <- "Anbo_2"
        mf_temp <- mf_treat_with_inhibOTUs %>%
            filter(G==inhib, toadID==indiv) %>%
            group_by(toadID, G, prepost, time) %>%
            summarise(reads=mean(reads)) %>%
            spread(key=prepost, value=reads)
        # stat_temp <- anova(lm(reads ~ prepost, data=mf_temp))
        #+ warning=FALSE, message=FALSE
        stat_temp <- wilcox.test(mf_temp$Pre, mf_temp$Pos)
        mf_treat_statdiff[which(mf_treat_statdiff$toadID==indiv & mf_treat_statdiff$G == inhib),"p"] <- (stat_temp$p.value)
        mf_treat_statdiff[which(mf_treat_statdiff$toadID==indiv & mf_treat_statdiff$G == inhib),"significant"] <- (stat_temp$p.value<0.05)
        
    }
}

mf_treat_statdiff %>%
    filter(!is.na(significant)) %>%
    ggplot(aes(x=G)) +
    geom_point(aes(y=temp, col=species, pch=significant), cex=2, position=position_jitter(width=0, height=0.25)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_shape_manual(values=c(21,19)) +
    geom_hline(aes(yintercept=0))


### TEsting presence absence BD

# TREAT
mf_treat_statdiff <- mf_treat_with_inhibOTUs %>%
    filter(prepost=="Pre" | (prepost == "Pos" & PABD == TRUE)) %>%
    group_by(toadID, species, prepost, G) %>%
    summarize(meanReads=mean(reads), meanBd = mean(eBD_log)) %>%
    dplyr::select(-meanBd) %>%
    spread(key=prepost, value=meanReads) %>%
    # mutate(diff=Pos-Pre, significant = NA, p=NA) %>%
    mutate(fc=log((Pos-Pre)/Pre +1), significant = NA, p=NA) %>%
    mutate(temp= ifelse(is.finite(fc), fc, ifelse(is.infinite(fc), 10, 0)))

#+ warning=FALSE, message=FALSE
allInhibTaxa <- unique(mf_treat_with_inhibOTUs$G)
allIndiv <- unique(mf_treat_with_inhibOTUs$toadID)
for ( inhib in allInhibTaxa) {
    # inhib <- allInhibTaxa[1]
    # inhib <- "Flavobacterium"
    for ( indiv in allIndiv) {
        if (exists("stat_temp")) {
            remove(stat_temp)
        }
        # indiv <- allIndiv[1]
        # indiv <- "Lipi_3"
        mf_temp <- mf_treat_with_inhibOTUs %>%
            filter(prepost=="Pre" | (prepost == "Pos" & PABD == TRUE)) %>%
            filter(G==inhib, toadID==indiv) %>%
            group_by(toadID, G, prepost, time) %>%
            summarise(reads=mean(reads)) %>%
            spread(key=prepost, value=reads)
        # stat_temp <- anova(lm(reads ~ prepost, data=mf_temp))
        if ( any(!is.na(mf_temp$Pre)) & any(!is.na(mf_temp$Pos)) ) {
            #+ warning=FALSE, message=FALSE
            stat_temp <- wilcox.test(mf_temp$Pre, mf_temp$Pos)
            mf_treat_statdiff[which(mf_treat_statdiff$toadID==indiv & mf_treat_statdiff$G == inhib),"p"] <- (stat_temp$p.value)
            mf_treat_statdiff[which(mf_treat_statdiff$toadID==indiv & mf_treat_statdiff$G == inhib),"significant"] <- (stat_temp$p.value<0.05)
        } else {
            mf_treat_statdiff[which(mf_treat_statdiff$toadID==indiv & mf_treat_statdiff$G == inhib),"p"] <- NA
            mf_treat_statdiff[which(mf_treat_statdiff$toadID==indiv & mf_treat_statdiff$G == inhib),"significant"] <- NA        
            }
        
        
    }
}

mf_treat_statdiff %>%
    # filter(!is.na(significant)) %>%
    ggplot(aes(x=G)) +
    geom_point(aes(y=temp, col=species, pch=significant), cex=2, position=position_jitter(width=0, height=0.25)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_shape_manual(values=c(21,19)) +
    geom_hline(aes(yintercept=0))



# mf_treat_with_inhibOTUs %>%
#     ggplot(aes(x=time, y=reads)) +
#     geom_bar(aes(fill=G), stat="identity") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.key.size =  unit(0.2, "cm") ) +
#     facet_wrap(~species, nrow=2) +
#     scale_fill_manual(values=set_col) +
#     geom_vline(aes(xintercept=5.5), col="grey", lty=2)

# mf_con_with_inhibOTUs %>%
#     group_by(G, toadID, time, species) %>%
#     summarize(aveReads = mean(reads)) %>%
#     filter(species=="Anbo") %>%
#     ggplot(aes(x=time, y=aveReads)) +
#     geom_bar(aes(fill=G), stat = "identity") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     facet_wrap(~toadID, nrow=10) +
#     scale_fill_manual(values=set_col)




# mf_con_with_inhibOTUs %>%
#     group_by(G, toadID, time, species) %>%
#     summarize(aveReads = mean(reads)) %>%
#     filter(species=="Anma") %>%
#     # unite(Sample, c(toadID,time), remove=FALSE) %>%
#     ggplot(aes(x=time, y=aveReads)) +
#     geom_bar(aes(fill=G), stat = "identity") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     facet_wrap(~toadID, nrow=10) +
#     scale_fill_manual(values=set_col)
# 
# mf_con_with_inhibOTUs %>%
#     group_by(G, toadID, time, species) %>%
#     summarize(aveReads = mean(reads)) %>%
#     filter(species=="Osse") %>%
#     ggplot(aes(x=time, y=aveReads)) +
#     geom_bar(aes(fill=G), stat = "identity") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     facet_wrap(~toadID, nrow=10) +
#     scale_fill_manual(values=set_col)
# 
# mf_con_with_inhibOTUs %>%
#     group_by(G, toadID, time, species) %>%
#     summarize(aveReads = mean(reads)) %>%
#     filter(species=="Lica") %>%
#     ggplot(aes(x=time, y=aveReads)) +
#     geom_bar(aes(fill=G), stat = "identity") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     facet_wrap(~toadID, nrow=10) +
#     scale_fill_manual(values=set_col)
# 
# mf_con_with_inhibOTUs %>%
#     group_by(G, toadID, time, species) %>%
#     summarize(aveReads = mean(reads)) %>%
#     filter(species=="Lipi") %>%
#     ggplot(aes(x=time, y=aveReads)) +
#     geom_bar(aes(fill=G), stat = "identity") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     facet_wrap(~toadID, nrow=10)+
#     scale_fill_manual(values=set_col) 




# 
# mf_treat_with_inhibOTUs %>%
#     group_by(G, toadID, time, species) %>%
#     summarize(aveReads = mean(reads)) %>%
#     filter(species=="Anbo") %>%
#     ggplot(aes(x=time, y=aveReads)) +
#     geom_bar(aes(fill=G), stat = "identity") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     facet_wrap(~toadID, nrow=10) +
#     scale_fill_manual(values=set_col)
# 
# mf_treat_with_inhibOTUs %>%
#     group_by(G, toadID, time, species) %>%
#     summarize(aveReads = mean(reads)) %>%
#     filter(species=="Anma") %>%
#     # unite(Sample, c(toadID,time), remove=FALSE) %>%
#     ggplot(aes(x=time, y=aveReads)) +
#     geom_bar(aes(fill=G), stat = "identity") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     facet_wrap(~toadID, nrow=10)+
#     scale_fill_manual(values=set_col)
# 
# mf_treat_with_inhibOTUs %>%
#     group_by(G, toadID, time, species) %>%
#     summarize(aveReads = mean(reads)) %>%
#     filter(species=="Osse") %>%
#     ggplot(aes(x=time, y=aveReads)) +
#     geom_bar(aes(fill=G), stat = "identity") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     facet_wrap(~toadID, nrow=10)+
#     scale_fill_manual(values=set_col)
# 
# mf_treat_with_inhibOTUs %>%
#     group_by(G, toadID, time, species) %>%
#     summarize(aveReads = mean(reads)) %>%
#     filter(species=="Lica") %>%
#     ggplot(aes(x=time, y=aveReads)) +
#     geom_bar(aes(fill=G), stat = "identity") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     facet_wrap(~toadID, nrow=10)+
#     scale_fill_manual(values=set_col)
# 
# mf_treat_with_inhibOTUs %>%
#     group_by(G, toadID, time, species) %>%
#     summarize(aveReads = mean(reads)) %>%
#     filter(species=="Lipi") %>%
#     ggplot(aes(x=time, y=aveReads)) +
#     geom_bar(aes(fill=G), stat = "identity") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     facet_wrap(~toadID, nrow=10)+
#     scale_fill_manual(values=set_col)


