#' ---
#' title: Statistical Analysis for 5Sp Dataset
#' author: Melissa Chen
#' output: github_document
#' 
#' ---
#' 

# Load packages
library(tidyverse)
library(rstanarm)
library(car) #Anova
library(vegan) # for permanova and bray curtis
library(gridExtra)
library(betareg) # for beta distr
library(lmtest) # for beta analysis
library(MASS) # for isMDS
# Mapping files
load("mf_con_without_init_infect.RData")
load("mf_treat_without_init_infect.RData")
load("mf.rare.RData")
# OTU table of inhibitory bacteria
load("otu.inhibOnly.treat.RData")
load("otu.inhibOnly.con.RData")
# Distance matrices
load("dm.filt.con.RData")
load("dm.filt.treat.RData")

# Previous analyses summaries
load("all_p.RData")
load("all_p_infected.RData")
load("all_p_withcon.RData")

# Inhibitory OTUs
load("mf_con_with_inhibOTUs.RData")
load("mf_treat_with_inhibOTUs.RData")

# OTU table for inhibitory OTUs
load("otu.inhibOnly.treat.RData")
load("otu.inhibOnly.con.RData")

# add a species column and PABD column
all_p <- all_p %>%
    mutate(PABD=ifelse(infect>0,1,0), infect = log(infect+1)) %>%
    rename(eBD_log=infect) %>%
    separate(toadID, into=c("species","indiv"), remove=FALSE) %>%
    mutate(p_dist =ifelse(is.na(exp_dist),NA,p_dist))
all_p_infected <- all_p_infected %>%
    mutate(PABD=ifelse(eBD_log>0,1,0)) %>%
    separate(toadID, into=c("species","indiv"), remove=FALSE)

### TESTING: Remove 0 values to test intensity
all_p <- all_p %>%
    mutate(eBD_log_infected = ifelse(eBD_log>0,eBD_log,NA))
all_p_infected <- all_p_infected %>%
    mutate(eBD_log_infected = ifelse(eBD_log>0,eBD_log,NA))

#### CURSORY GLANCE AT DATA ####
gg_NMDS <- mf_con_without_init_infect %>%
    ggplot(aes(x=NMDS1, y=NMDS2)) +
    geom_point(aes(col=species), cex=3, show.legend = FALSE)
gg_infect <- mf_treat_without_init_infect  %>%
    ggplot(aes(x=species, y=eBD_log)) +
    geom_point(aes(col=species), cex=3, position = position_jitter(width=0.1, height=0.05), show.legend = FALSE)

temp1 <- mf_con_without_init_infect %>%
    dplyr::select(species, logRich) %>%
    mutate(metric="log_OTU_Richness") %>%
    rename(value=logRich)
temp2 <- mf_con_without_init_infect %>%
    dplyr::select(species, inhibRich) %>%
    mutate(metric="Inhibitory_OTU_Richness")%>%
    rename(value=inhibRich)
temp3 <- mf_con_without_init_infect %>%
    dplyr::select(species, percInhib) %>%
    mutate(metric="Percent_Inhibitory")%>%
    rename(value=percInhib)
temp4 <- mf_con_without_init_infect %>%
    dplyr::select(species, disper_bray_curtis) %>%
    mutate(metric="Dispersion_from_centroid")%>%
    rename(value=disper_bray_curtis)
temp5 <- mf_con_without_init_infect %>%
    dplyr::select(species, distance_bray_curtis) %>%
    mutate(metric="Distance_from_previous_timepoint")%>%
    rename(value=distance_bray_curtis)


gg_all <- rbind(temp1,temp2,temp3,temp4, temp5) %>%
    rename(Species=species) %>%
    mutate(Metric = gsub("_"," ",metric, fixed=TRUE)) %>%
    mutate(Metric = factor(Metric, levels=c("log OTU Richness","Dispersion from centroid", "Distance from previous timepoint","Inhibitory OTU Richness","Percent Inhibitory"))) %>%
    ggplot(aes(x=Species, y=value)) +
    geom_boxplot() +
    geom_point(aes(col=Species), position = position_jitter(width=0.1, height=0), alpha=1/3)+
    facet_grid(Metric~., scales = "free", switch="y") +
    ylab("")+
    xlab("Species") 
lay <- rbind(c(1,2),
             c(3,2))

#+ fig.height=12, fig.width=10
grid.arrange(gg_NMDS, gg_all, gg_infect, layout_matrix = lay)

#### Stats ####
does_comp_differ_btwn_sp_and_across_time_con <- adonis2(dist(dm.filt.con) ~ species*time, data=mf_con_without_init_infect)
does_comp_differ_btwn_sp_and_across_time_con
beta_con_main_p <- does_comp_differ_btwn_sp_and_across_time_con$`Pr(>F)`[1:2]
beta_con_main_df1 <- does_comp_differ_btwn_sp_and_across_time_con$Df[1:2]
beta_con_main_df2 <- does_comp_differ_btwn_sp_and_across_time_con$Df[4]

beta_con_interaction_p <- does_comp_differ_btwn_sp_and_across_time_con$`Pr(>F)`[3]
beta_con_interaction_df1 <- does_comp_differ_btwn_sp_and_across_time_con$Df[3]
beta_con_interaction_df2 <- does_comp_differ_btwn_sp_and_across_time_con$Df[4]

beta_con_main_f <- does_comp_differ_btwn_sp_and_across_time_con$`F`[1:2]
beta_con_interaction_f <- does_comp_differ_btwn_sp_and_across_time_con$`F`[3]


mf_treat_without_init_infect_post <- mf_treat_without_init_infect %>%
    filter(prepost == "Pos")
does_comp_differ_btwn_sp_and_across_time_and_infect_treat <- adonis2(dist(dm.filt.treat) ~ species*time*PABD, data=mf_treat_without_init_infect_post)
does_comp_differ_btwn_sp_and_across_time_and_infect_treat
beta_treat_main_p <- does_comp_differ_btwn_sp_and_across_time_and_infect_treat$`Pr(>F)`[1:2]
beta_treat_main_df1 <- does_comp_differ_btwn_sp_and_across_time_and_infect_treat$Df[1:2]
beta_treat_main_df2 <- does_comp_differ_btwn_sp_and_across_time_and_infect_treat$Df[8]

beta_treat_interaction_p <- does_comp_differ_btwn_sp_and_across_time_and_infect_treat$`Pr(>F)`[4]
beta_treat_interaction_df1 <- does_comp_differ_btwn_sp_and_across_time_and_infect_treat$Df[4]
beta_treat_interaction_df2 <- does_comp_differ_btwn_sp_and_across_time_and_infect_treat$Df[8]

beta_treat_main_f <- does_comp_differ_btwn_sp_and_across_time_and_infect_treat$`F`[1:2]
beta_treat_interaction_f <- does_comp_differ_btwn_sp_and_across_time_and_infect_treat$`F`[4]

beta_con_time_eff <- ""
beta_treat_time_eff <- ""


#' There is a significant effect of species, time, and PABD; all of these things also significantly interact
#' EXCEPT species and PABD and all 3 together, which suggests species microbiomes change in the "same way" when infected


#### Preliminary stats on broad patterns ####

### RICHNESS AND TIME ###

#' Does richness change over time in control individuals?
# Type I ANOVA to test for interaction-- (AB | A, B)
rich_con_interaction_lm <- lm(logRich ~ species*time, data=mf_con_without_init_infect)
rich_con_interaction <- anova(rich_con_interaction_lm)
rich_con_interaction
# Use Type II ANOVA (no interaction present)
rich_con_main_lm <- lm(logRich ~ species + time, data=mf_con_without_init_infect)
rich_con_main <- Anova(rich_con_main_lm, type = 2)
rich_con_main
# There is a significant effect of species but not time or interaction

#' Does richness change over time in treatment individuals?
# Type I ANOVA to test for interaction (AB | A,B)
rich_treat_interaction_lm <- lm(logRich ~ species*time, data=mf_treat_without_init_infect)
rich_treat_interaction <- anova(rich_treat_interaction_lm)
rich_treat_interaction
# Type III ANOVA (valid in presence of interaction)
rich_treat_main_lm <- lm(logRich ~ species * time, data=mf_treat_without_init_infect, contrasts=list(species=contr.sum))
rich_treat_main <- Anova(rich_treat_main_lm, type=3)
rich_treat_main

### DISTANCE TO CENTROID AND TIME ####
#' Is there an effect of species and time on controls?
# Type I ANOVA (to check for interaction) (AB | A,B)
centroid_con_interaction_lm <- lm(log(disper_bray_curtis) ~ species*time, data=mf_con_without_init_infect)
centroid_con_interaction <- anova(centroid_con_interaction_lm)
centroid_con_interaction
# Type II ANOVA with no interaction
centroid_con_main_lm <- lm(log(disper_bray_curtis) ~ species + time, data=mf_con_without_init_infect)
centroid_con_main <- Anova(centroid_con_main_lm, type = 2)
centroid_con_main

#' Is there an effect of species and time on treatment??
# Type I ANOVA (to check for interaction) (AB | A,B)
centroid_treat_interaction_lm <- lm(log(disper_bray_curtis) ~ species*time, data=mf_treat_without_init_infect)
centroid_treat_interaction <- anova(centroid_treat_interaction_lm)
centroid_treat_interaction
# Type II ANOVA with no interaction
centroid_treat_main_lm <- lm(log(disper_bray_curtis) ~ species + time, data=mf_treat_without_init_infect)
centroid_treat_main <- Anova(centroid_treat_main_lm, type = 2)
centroid_treat_main

### DISPERSAL AND TIME ###
#' Is there an effect of species and time on controls?
# Type I ANOVA (to check for interaction) (AB | A,B)
# disp_con_interaction_glm <- betareg(distance_bray_curtis ~ species*time, data=mf_con_without_init_infect)
# disp_con_interaction_lm <- lm(distance_bray_curtis ~ species*time, data=mf_con_without_init_infect)

disp_con_interaction_glmmain <- betareg(distance_bray_curtis ~ species + time, data=mf_con_without_init_infect)
disp_con_interaction_glminter <- betareg(distance_bray_curtis ~ species*time, data=mf_con_without_init_infect)
disp_con_interaction <- lrtest(disp_con_interaction_glmmain,disp_con_interaction_glminter) # plus species

# Type II ANOVA with no interaction

disp_con_main_glmtonly <- betareg(distance_bray_curtis ~ time, data=mf_con_without_init_infect)
disp_con_main_glmsponly <- betareg(distance_bray_curtis ~ species, data=mf_con_without_init_infect)

disp_con_main_sp <- lrtest(disp_con_main_glmtonly,disp_con_interaction_glmmain) # plus species
disp_con_main_time <- lrtest(disp_con_main_glmsponly,disp_con_interaction_glmmain) # plus time

#' Is there an effect of species and time on treatment??
# Type I ANOVA (to check for interaction) (AB | A,B)
# disp_treat_interaction_lm <- lm(distance_bray_curtis ~ species*time, data=mf_treat_without_init_infect)
# disp_treat_interaction <- anova(disp_treat_interaction_lm)
# disp_treat_interaction
disp_treat_interaction_glmmain <- betareg(distance_bray_curtis ~ species + time, data=mf_treat_without_init_infect)
disp_treat_interaction_glminter <- betareg(distance_bray_curtis ~ species*time, data=mf_treat_without_init_infect)
disp_treat_interaction <- lrtest(disp_treat_interaction_glmmain,disp_treat_interaction_glminter) # plus species
disp_treat_interaction
# Type II ANOVA with no interaction

disp_treat_main_glmtonly <- betareg(distance_bray_curtis ~ time, data=mf_treat_without_init_infect)
disp_treat_main_glmsponly <- betareg(distance_bray_curtis ~ species, data=mf_treat_without_init_infect)

disp_treat_main_sp <- lrtest(disp_treat_main_glmtonly,disp_treat_interaction_glmmain) # plus species
disp_treat_main_time <- lrtest(disp_treat_main_glmsponly,disp_treat_interaction_glmmain) # plus time
disp_treat_main_sp
disp_treat_main_time

### PERCENT INHIB ###
#' Does percent inhibitory change with species or time?
# # Type I ANOVA (to test for interaction) in control group?
# pinhib_con_interaction_glm <- glm(percInhib ~ species*time, family = binomial(), data=mf_con_without_init_infect, weights=mf_con_without_init_infect$n)
# pinhib_con_interaction <- anova(pinhib_con_interaction_glm, test = "Chisq")
# pinhib_con_interaction
# # Type III ANOVA (to test for main effects, given interaction) in control group?
# pinhib_con_main_glm <- glm(percInhib ~ species*time, family = binomial(), data=mf_con_without_init_infect, weights=mf_con_without_init_infect$n, contrasts=list(species=contr.sum))
# pinhib_con_main <- Anova(pinhib_con_main_glm, type=3)
# pinhib_con_main
# 
# # Does percent inhibitory change with species or time in treatment group?
# # Type I ANOVA (to test for interaction) in control group?
# pinhib_treat_interaction_glm <- glm(percInhib ~ species*time, family = binomial(), data=mf_treat_without_init_infect, weights=mf_treat_without_init_infect$n)
# pinhib_treat_interaction <- anova(pinhib_treat_interaction_glm, test = "Chisq")
# pinhib_treat_interaction
# # Type III ANOVA (to test for main effects, given interaction) in control group?
# pinhib_treat_main_glm <- glm(percInhib ~ species*time, family = binomial(), data=mf_treat_without_init_infect, weights=mf_treat_without_init_infect$n, contrasts = list(species=contr.sum))
# pinhib_treat_main <- Anova(pinhib_treat_main_glm, type=3)
# pinhib_treat_main

# There is ONE zero... so need to transform data to be between 0 and 1. Super annoying. 
y.transf.betareg <- function(y){
    n.obs <- sum(!is.na(y))
    (y * (n.obs - 1) + 0.5) / n.obs
}

# Type I ANOVA (to check for interaction) (AB | A,B)
pinhib_con_interaction_glmmain <- betareg(y.transf.betareg(percInhib) ~ species + time, data=mf_con_without_init_infect)
pinhib_con_interaction_glminter <- betareg(y.transf.betareg(percInhib) ~ species*time, data=mf_con_without_init_infect)
pinhib_con_interaction <- lrtest(pinhib_con_interaction_glmmain,pinhib_con_interaction_glminter) # plus species
pinhib_con_interaction

# Type III ANOVA to test for main effects

pinhib_con_main_glmtonly <- betareg(y.transf.betareg(percInhib) ~ time , data=mf_con_without_init_infect)
pinhib_con_main_glmsponly <- betareg(y.transf.betareg(percInhib) ~ species , data=mf_con_without_init_infect)

pinhib_con_main_sp <- lrtest(pinhib_con_main_glmtonly,pinhib_con_interaction_glmmain) # plus species
pinhib_con_main_time <- lrtest(pinhib_con_main_glmsponly,pinhib_con_interaction_glmmain) # plus time
pinhib_con_main_sp
pinhib_con_main_time
#' Is there an effect of species and time on treatment??
# Type I ANOVA (to check for interaction) (AB | A,B)
# pinhib_treat_interaction_lm <- lm(percInhib ~ species*time, data=mf_treat_without_init_infect)
# pinhib_treat_interaction <- anova(pinhib_treat_interaction_lm)
# pinhib_treat_interaction
#

pinhib_treat_interaction_glmmain <- betareg(y.transf.betareg(percInhib) ~ species + time, data=mf_treat_without_init_infect)
pinhib_treat_interaction_glminter <- betareg(y.transf.betareg(percInhib) ~ species*time, data=mf_treat_without_init_infect)
pinhib_treat_interaction <- lrtest(pinhib_treat_interaction_glmmain,pinhib_treat_interaction_glminter) # plus species
pinhib_treat_interaction
# Type II ANOVA with no interaction

pinhib_treat_main_glmtonly <- betareg(y.transf.betareg(percInhib) ~ time, data=mf_treat_without_init_infect)
pinhib_treat_main_glmsponly <- betareg(y.transf.betareg(percInhib) ~ species, data=mf_treat_without_init_infect)

pinhib_treat_main_sp <- lrtest(pinhib_treat_main_glmtonly,pinhib_treat_interaction_glmmain) # plus species
pinhib_treat_main_time <- lrtest(pinhib_treat_main_glmsponly,pinhib_treat_interaction_glmmain) # plus time
pinhib_treat_main_sp
pinhib_treat_main_time

### INHIB RICH ###
# Does richness of inhibitory bacteria differ betwen species and time points?
# Type I ANOVA to test for interactions in control
inhibRich_con_interaction_glm <- glm(inhibRich ~ species*time, data=mf_con_without_init_infect, family=poisson())
inhibRich_con_interaction <- anova(inhibRich_con_interaction_glm, test="Chisq")
inhibRich_con_interaction
# TYpe III ANOVA to test for main effects with interactions in control
inhibRich_con_main_glm <- glm(inhibRich ~ species*time, data=mf_con_without_init_infect, family=poisson(), contrasts=list(species=contr.sum))
inhibRich_con_main <- Anova(inhibRich_con_main_glm,type=3)
inhibRich_con_main

# Type I ANOVA to test for interactions
inhibRich_treat_interaction_glm <- glm(inhibRich ~ species*time, data=mf_treat_without_init_infect, family=poisson())
inhibRich_treat_interaction <- anova(inhibRich_treat_interaction_glm, test="Chisq")
inhibRich_treat_interaction
# TYpe III ANOVA to test for main effects with interactions
inhibRich_treat_main_glm <- glm(inhibRich ~ species*time, data=mf_treat_without_init_infect, family=poisson(), contrasts = list(species=contr.sum))
inhibRich_treat_main <- Anova(inhibRich_treat_main_glm,type=3)
inhibRich_treat_main

#### Summarize overall trends into table ####

# RICHNESS
rich_con_main_p <- rich_con_main$`Pr(>F)`[1:2]
rich_con_main_df1 <- rich_con_main$Df[1:2]
rich_con_main_df2 <- rich_con_main$Df[3]

rich_con_interaction_p <- rich_con_interaction$`Pr(>F)`[3]
rich_con_interaction_df1 <- rich_con_interaction$Df[3]
rich_con_interaction_df2 <- rich_con_interaction$Df[4]

rich_treat_main_p <- rich_treat_main$`Pr(>F)`[2:3]
rich_treat_main_df1 <- rich_treat_main$Df[2:3]
rich_treat_main_df2 <- rich_treat_main$Df[5]

rich_treat_interaction_p <- rich_treat_interaction$`Pr(>F)`[3]
rich_treat_interaction_df1 <- rich_treat_interaction$Df[3]
rich_treat_interaction_df2 <- rich_treat_interaction$Df[4]

rich_con_main_f <- rich_con_main$`F value`[1:2]
rich_con_interaction_f <- rich_con_interaction$`F value`[3]

rich_treat_main_f <- rich_treat_main$`F value`[2:3]
rich_treat_interaction_f <- rich_treat_interaction$`F value`[3]

rich_con_time_eff <- ifelse(rich_con_main_lm$coefficients["time"]>0,"(+)","(-)")
rich_treat_time_eff <- ifelse(rich_treat_main_lm$coefficients["time"]>0,"(+)","(-)")

# CENTROID
centroid_con_main_p <- centroid_con_main$`Pr(>F)`[1:2]
centroid_con_main_df1 <- centroid_con_main$Df[1:2]
centroid_con_main_df2 <- centroid_con_main$Df[3]

centroid_con_interaction_p <- centroid_con_interaction$`Pr(>F)`[3]
centroid_con_interaction_df1 <- centroid_con_interaction$Df[3]
centroid_con_interaction_df2 <- centroid_con_interaction$Df[4]

centroid_treat_main_p <- centroid_treat_main$`Pr(>F)`[1:2]
centroid_treat_main_df1 <- centroid_treat_main$Df[1:2]
centroid_treat_main_df2 <- centroid_treat_main$Df[3]

centroid_treat_interaction_p <- centroid_treat_interaction$`Pr(>F)`[3]
centroid_treat_interaction_df1 <- centroid_treat_interaction$Df[3]
centroid_treat_interaction_df2 <- centroid_treat_interaction$Df[4]

centroid_con_main_f <- centroid_con_main$`F value`[1:2]
centroid_con_interaction_f <- centroid_con_interaction$`F value`[3]

centroid_treat_main_f <- centroid_treat_main$`F value`[1:2]
centroid_treat_interaction_f <- centroid_treat_interaction$`F value`[3]

centroid_con_time_eff <- ifelse(centroid_con_main_lm$coefficients["time"]>0,"(+)","(-)")
centroid_treat_time_eff <- ifelse(centroid_treat_main_lm$coefficients["time"]>0,"(+)","(-)")


# DISPERSION
disp_con_main_p <- c(disp_con_main_sp$`Pr(>Chisq)`[2], disp_con_main_time$`Pr(>Chisq)`[2])
# disp_con_main_p <- disp_con_main$`Pr(>F)`[1:2]
disp_con_main_df1 <- c(disp_con_main_sp$Df[2], disp_con_main_time$Df[2])
disp_con_main_n <- sum(!is.na(mf_con_without_init_infect$distance_bray_curtis))-1-7

disp_con_interaction_p <- disp_con_interaction$`Pr(>Chisq)`[2]
disp_con_interaction_df1 <- disp_con_interaction$Df[2]
disp_con_interaction_n <- sum(!is.na(mf_con_without_init_infect$distance_bray_curtis))-1-11

disp_treat_main_p <- c(disp_treat_main_sp$`Pr(>Chisq)`[2], disp_treat_main_time$`Pr(>Chisq)`[2])
# disp_treat_main_p <- disp_treat_main$`Pr(>F)`[1:2]
disp_treat_main_df1 <- c(disp_treat_main_sp$Df[2], disp_treat_main_time$Df[2])
disp_treat_main_n <- sum(!is.na(mf_treat_without_init_infect$distance_bray_curtis))-1-7

disp_treat_interaction_p <- disp_treat_interaction$`Pr(>Chisq)`[2]
disp_treat_interaction_df1 <- disp_treat_interaction$Df[2]
disp_treat_interaction_n <- sum(!is.na(mf_treat_without_init_infect$distance_bray_curtis))-1-11

disp_con_main_f <- c(disp_con_main_sp$Chisq[2], disp_con_main_time$Chisq[2])
disp_con_interaction_f <- disp_con_interaction$Chisq[2]

disp_treat_main_f <- c(disp_treat_main_sp$Chisq[2], disp_treat_main_time$Chisq[2])
disp_treat_interaction_f <- disp_treat_interaction$Chisq[2]

# disp_con_time_eff <- ifelse(disp_con_main_lm$coefficients["time"]>0,"(+)","(-)")
disp_con_time_eff <- ifelse(disp_con_main_glmtonly$coefficients$mean[2]>0,"(+)","(-)")
# disp_treat_time_eff <- ifelse(disp_treat_main_lm$coefficients["time"]>0,"(+)","(-)")
disp_treat_time_eff <- ifelse(disp_treat_main_glmtonly$coefficients$mean[2]>0,"(+)","(-)")


# INHIB RICH
inhibRich_con_main_p <- inhibRich_con_main$`Pr(>Chisq)`[1:2]
inhibRich_con_main_df1 <- inhibRich_con_main$Df[1:2]

inhibRich_con_interaction_p <- inhibRich_con_interaction$`Pr(>Chi)`[4]
inhibRich_con_interaction_df1 <- inhibRich_con_interaction$Df[4]
inhibRich_con_interaction_n <- inhibRich_con_interaction$`Resid. Df`[1]

inhibRich_treat_main_p <- inhibRich_treat_main$`Pr(>Chisq)`[1:2]
inhibRich_treat_main_df1 <- inhibRich_treat_main$Df[1:2]

inhibRich_treat_interaction_p <- inhibRich_treat_interaction$`Pr(>Chi)`[4]
inhibRich_treat_interaction_df1 <- inhibRich_treat_interaction$Df[4]
inhibRich_treat_interaction_n <- inhibRich_treat_interaction$`Resid. Df`[1]

inhibRich_con_main_f <- inhibRich_con_main$`LR Chisq`[1:2]
inhibRich_con_interaction_f <- inhibRich_con_interaction$Deviance[4]

inhibRich_treat_main_f <- inhibRich_treat_main$`LR Chisq`[1:2]
inhibRich_treat_interaction_f <- inhibRich_treat_interaction$Deviance[4]

inhibRich_con_time_eff <- ifelse(inhibRich_con_main_glm$coefficients["time"]>0,"(+)","(-)")
inhibRich_treat_time_eff <- ifelse(inhibRich_treat_main_glm$coefficients["time"]>0,"(+)","(-)")


# PERCENT
pinhib_con_main_p <- c(pinhib_con_main_sp$`Pr(>Chisq)`[2], pinhib_con_main_time$`Pr(>Chisq)`[2])
# pinhib_con_main_p <- pinhib_con_main$`Pr(>F)`[1:2]
pinhib_con_main_df1 <- c(pinhib_con_main_sp$Df[2], pinhib_con_main_time$Df[2])
pinhib_con_main_n <- sum(!is.na(mf_con_without_init_infect$distance_bray_curtis))-1-7

pinhib_con_interaction_p <- pinhib_con_interaction$`Pr(>Chisq)`[2]
pinhib_con_interaction_df1 <- pinhib_con_interaction$Df[2]
pinhib_con_interaction_n <- sum(!is.na(mf_con_without_init_infect$distance_bray_curtis))-1-11

pinhib_treat_main_p <- c(pinhib_treat_main_sp$`Pr(>Chisq)`[2], pinhib_treat_main_time$`Pr(>Chisq)`[2])
# pinhib_treat_main_p <- pinhib_treat_main$`Pr(>F)`[1:2]
pinhib_treat_main_df1 <- c(pinhib_treat_main_sp$Df[2], pinhib_treat_main_time$Df[2])
pinhib_treat_main_n <- sum(!is.na(mf_treat_without_init_infect$distance_bray_curtis))-1-7

pinhib_treat_interaction_p <- pinhib_treat_interaction$`Pr(>Chisq)`[2]
pinhib_treat_interaction_df1 <- pinhib_treat_interaction$Df[2]
pinhib_treat_interaction_n <- sum(!is.na(mf_treat_without_init_infect$distance_bray_curtis))-1-11

pinhib_con_main_f <- c(pinhib_con_main_sp$Chisq[2], pinhib_con_main_time$Chisq[2])
pinhib_con_interaction_f <- pinhib_con_interaction$Chisq[2]

pinhib_treat_main_f <- c(pinhib_treat_main_sp$Chisq[2], pinhib_treat_main_time$Chisq[2])
pinhib_treat_interaction_f <- pinhib_treat_interaction$Chisq[2]

# pinhib_con_time_eff <- ifelse(pinhib_con_main_lm$coefficients["time"]>0,"(+)","(-)")
pinhib_con_time_eff <- ifelse(pinhib_con_main_glmtonly$coefficients$mean[2]>0,"(+)","(-)")
# pinhib_treat_time_eff <- ifelse(pinhib_treat_main_lm$coefficients["time"]>0,"(+)","(-)")
pinhib_treat_time_eff <- ifelse(pinhib_treat_main_glmtonly$coefficients$mean[2]>0,"(+)","(-)")

stat_results <- as.data.frame(matrix(ncol=5, nrow=12, dimnames = list(NULL,c("Microbiome metric","Control or Treatment","Main effect: species","Main effect: time", "Interaction: species x time"))), check.names=FALSE)
stat_results$`Microbiome metric` <- c("Beta Diversity"
                                      , "Beta Diversity"
                                      , "OTU Richness"
                                      , "OTU RIchness"
                                      , "Distance to centroid"
                                      , "Distance to centroid"
                                      , "Stability (BC distance)"
                                      , "Stability (BC distance)"
                                      , "Percent Inhibitory"
                                      , "Percent Inhibitory"
                                      ,"Inhibitory Richness"
                                      ,"Inhibitory Richness"
                                      )
stat_results$`Control or Treatment` <- rep(c("Control","Treatment"), 6)
current_row <- 1

for ( test in c("beta","rich","centroid","disp","pinhib","inhibRich") ) {
    for ( ct in c("con","treat")) {
        if (test %in% c("beta","rich","centroid")) {
            df1_sp <- get(paste(test, ct, "main_df1", sep="_"))[1] # species
            df2_sp <- get(paste(test, ct, "main_df2", sep="_")) # species
            
            df1_t <- get(paste(test, ct, "main_df1", sep="_"))[2] # time
            df2_t <- get(paste(test, ct, "main_df2", sep="_")) # time
            
            df1_inter <- get(paste(test, ct, "interaction_df1", sep="_"))
            df2_inter <- get(paste(test, ct, "interaction_df2", sep="_"))
            
            stat_main_sp<- paste0(", F(",df1_sp,",",df2_sp,")=")
            stat_main_t<- paste0(", F(",df1_t,",",df2_t,")=")
            
            stat_interaction <- paste0(", F(",df1_inter,",",df2_inter,")=")
            
        } else if ( test %in% c("disp","pinhib")) {
            df1_sp <- get(paste(test, ct, "main_df1", sep="_"))[1] # species
            n_sp <- get(paste(test, ct, "main_n", sep="_")) # species
            
            df1_t <- get(paste(test, ct, "main_df1", sep="_"))[2] # time
            n_t <- get(paste(test, ct, "main_n", sep="_")) # time
            
            df1_inter <- get(paste(test, ct, "interaction_df1", sep="_"))
            n_inter <- get(paste(test, ct, "interaction_n", sep="_"))
            
            stat_main_sp<- paste0(", Chisq(",df1_sp,",N=",df2_sp,")=")
            stat_main_t<- paste0(", Chisq(",df1_t,",N=",df2_t,")=")
            
            stat_interaction <- paste0(", Chisq(",df1_inter,",N=",df2_inter,")=")
            
        } else {
            
            df1_sp <- get(paste(test, ct, "main_df1", sep="_"))[1] # species

            df1_t <- get(paste(test, ct, "main_df1", sep="_"))[2] # time

            df1_inter <- get(paste(test, ct, "interaction_df1", sep="_"))
            n_inter <- get(paste(test, ct, "interaction_n", sep="_"))
            
            stat_main_sp <- paste0(", Chisq(",df1_sp,")=")
            stat_main_t <- paste0(", Chisq(",df1_t,")=")
            stat_interaction <- paste0(", LRChi(",df1_inter,",N=",n_inter,")=")
        }
        
            stat_results[current_row, 3:5] <- c(paste0("p=", signif(get(paste(test, ct, "main_p", sep="_"))[1],3), stat_main_sp, signif(get(paste(test, ct, "main_f", sep="_"))[1],2))
              , paste0("p=", signif(get(paste(test, ct, "main_p", sep="_"))[2],3), stat_main_t, signif(get(paste(test, ct, "main_f", sep="_"))[2],3),get(paste(test, ct, "time_eff", sep="_")) )
              , paste0("p=", signif(get(paste(test, ct, "interaction_p", sep="_")),3), stat_interaction, signif(get(paste(test, ct, "interaction_f", sep="_"))[1],2))
            )
        
        current_row <- current_row+1
    }
}


stat_results
write_csv(stat_results, path = "stats_table.csv")


#### PART I ####
#' # PART I 
#' ## Part I: Microbiome state and effect on infection risk and intensity\
#' \

# CREATE TABLE for statistics
all_stats <- matrix(ncol=7, nrow=6, dimnames = list(c("","ASV Richness","Inhibitory ASV Richness","Proportion Inhibitory","Dispersion","Instability")
                                       ,c("","","Effect on infection risk","","","Effect of infection","")))
all_stats[1,] <- c("","Infection presence/absence","Infection intensity (with zeros)","Infection intensity (no zeros)","Infection presence/absence","Infection intensity (with zeros)","Infection intensity (no zeros)")
# all_stats[,1] <- c("","ASV Richness","Inhibitory ASV Richness","Proportion Inhibitory","Dispersion","Instability")
#### PABD and OTU Richness ####

#' ### (1a) Does overall diversity of microbiome influence BD infection rate?\
#' The first thing we would like to know is whether microbiome richness of an individual influences 
#' its risk of becoming infected by BD. The most simple way to look at this would be
#' to plot OTU richness VS presence/absence of BD
#' Below, we fitted normal and lognormal distributions, respectively, to diversity (shannon) and
#' otu richness to individuals prior to BD infection. Now, we fit a binomial general
#' linearized model to see if there is a relationship between diversity and infection rate.
#' \

glm_PABD_prich <- glm(PABD ~ species*p_logRich, data=all_p, family=binomial(link="logit"))
PABD_rich_inter <- anova(glm_PABD_prich, test="Chisq") # test for interaction
PABD_rich_inter
PABD_rich_main <- Anova(glm_PABD_prich, type=2) # test for main effects
PABD_rich_main
all_p %>%
    ggplot(aes(x=p_logRich, y=PABD)) +
    geom_point(aes(col=species), cex=3)  

#' If anything, it looks like increased diversity and richness might increase infection risk


#'\

#### eBD and OTU Richness ####

#' (1b) Does overall diversity of microbiome influence BD infection intensity?\
#' \
#' The next thing we would like to know is if richness of the microbiome
#' influences infection intensity.
#' \

#' Now let's do richness-- without zeros
lm_eBD_prich <- lm(eBD_log_infected ~ species*p_logRich, data=all_p)
eBD_rich_inter1 <- anova(lm_eBD_prich)
eBD_rich_inter1
eBD_rich_main1 <- Anova(lm_eBD_prich, type=2)
eBD_rich_main1
all_p %>%
    ggplot(aes(x=p_logRich, y=eBD_log_infected)) +
    geom_point(aes(col=species), cex=3) 

#' Try a version where we keep zeros
lm_eBD_prich_nozeros <- lm(eBD_log ~ species*p_logRich, data=all_p)
eBD_rich_inter <- anova(lm_eBD_prich_nozeros)
eBD_rich_inter
eBD_rich_main <- Anova(lm_eBD_prich_nozeros, type=2)
eBD_rich_main
all_p %>%
    ggplot(aes(x=p_logRich, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 

all_stats["ASV Richness",1:4] <- c("ASV Richness"
                                   , paste0("p=",signif(PABD_rich_main$`Pr(>Chisq)`[2],3), " LR Chisq(",PABD_rich_main$Df[2],",N=",nrow(all_p),")=",signif(PABD_rich_main$`LR Chisq`[2],3))
                                   , paste0("p=",signif(eBD_rich_main$`Pr(>F)`[2],3), ", F(",eBD_rich_main$Df[2],",",eBD_rich_main$Df[4],")=",signif(eBD_rich_main$`F value`[2],3))
                                   , paste0("p=",signif(eBD_rich_main1$`Pr(>F)`[2],3), ", F(",eBD_rich_main1$Df[2],",",eBD_rich_main1$Df[4],")=",signif(eBD_rich_main1$`F value`[2],3))
)


#### PABD and Instability  ####

#' (2a) Does instability of microbiome influence BD infection rate?\
#' Here we look at average distance travelled (bray-curtis) between samples
#' prior to being infected. We see if it is correlated to infection risk.\

glm_PABD_pbc <- glm(PABD ~ species*p_dist, data=all_p, family=binomial)
PABD_dist_inter <- anova(glm_PABD_pbc, test="Chisq")
PABD_dist_inter
PABD_dist_main <- Anova(glm_PABD_pbc, type=2)
PABD_dist_main
all_p %>%
    ggplot(aes(x=p_dist, y=PABD)) +
    geom_point(aes(col=species), cex=3) 

#### eBD and Instability ####
 
#' \
#' (2b) Does instability of microbiome influence BD infection intensity?\

#' A version without zeros
lm_BD_pbc <- lm(eBD_log_infected ~ species*p_dist, data=all_p)
eBD_dist_inter1 <- anova(lm_BD_pbc)
eBD_dist_inter1
eBD_dist_main1 <- Anova(lm_BD_pbc, type=2)
eBD_dist_main1
all_p %>%
    ggplot(aes(x=p_dist, y=eBD_log_infected)) +
    geom_point(aes(col=species), cex=3) 

#' Try a version where we keep zeros
lm_BD_pbc_nozeros <- lm(eBD_log ~ species*p_dist, data=all_p)
eBD_dist_inter <- anova(lm_BD_pbc_nozeros)
eBD_dist_inter
eBD_dist_main <- Anova(lm_BD_pbc_nozeros, type=2)
eBD_dist_main
all_p %>%
    ggplot(aes(x=p_dist, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 




all_stats["Instability",1:4] <- c("Instability"
                                  , paste0("p=",signif(PABD_dist_main$`Pr(>Chisq)`[2],3), " LR Chisq(",PABD_dist_main$Df[2],",N=",sum(!is.na(all_p$p_dist)),")=",signif(PABD_dist_main$`LR Chisq`[2],3))
                                  , paste0("p=",signif(eBD_dist_main$`Pr(>F)`[2],3), ", F(",eBD_dist_main$Df[2],",",eBD_dist_main$Df[4],")=",signif(eBD_dist_main$`F value`[2],3))
                                  , paste0("p=",signif(eBD_dist_main1$`Pr(>F)`[2],3), ", F(",eBD_dist_main1$Df[2],",",eBD_dist_main1$Df[4],")=",signif(eBD_dist_main1$`F value`[2],3))
)


#### PABD and Dispersion ####

#' (2a) Does dispersion of microbiome influence BD infection rate?\
#' Here we look at average distance to centroid (bray-curtis) between samples
#' prior to being infected at same time point. We see if it is correlated to infection risk.

glm_PABD_pdist <- glm(PABD ~ species*p_disper, data=all_p, family=binomial)
PABD_disper_inter <- anova(glm_PABD_pdist, test="Chisq")
PABD_disper_inter
PABD_disper_main <- Anova(glm_PABD_pdist, type=2)
PABD_disper_main
all_p %>%
    ggplot(aes(x=p_disper, y=PABD)) +
    geom_point(aes(col=species), cex=3) 


#### eBD and Dispersion ####

#' 
#' \
#' (2b) Does dispersion of microbiome influence BD infection intensity?\
#' 
#' A version wtihout zeros
lm_BD_pdist <- lm(eBD_log_infected ~ species*p_disper, data=all_p)
eBD_disper_inter1 <- anova(lm_BD_pdist)
eBD_disper_inter1
eBD_disper_main1 <- Anova(lm_BD_pdist, type=2)
eBD_disper_main1
all_p %>%
    ggplot(aes(x=p_disper, y=eBD_log_infected)) +
    geom_point(aes(col=species), cex=3) 

#' Try a version where we keep zeros
lm_BD_pdist_nozeros <- lm(eBD_log ~ species*p_disper, data=all_p)
eBD_disper_inter <- anova(lm_BD_pdist_nozeros)
eBD_disper_inter
eBD_disper_main <- Anova(lm_BD_pdist_nozeros, type=2)
eBD_disper_main
all_p %>%
    ggplot(aes(x=p_disper, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 


all_stats["Dispersion",1:4] <- c("Dispersion"
                                 , paste0("p=",signif(PABD_disper_main$`Pr(>Chisq)`[2],3), " LR Chisq(",PABD_disper_main$Df[2],",N=",nrow(all_p),")=",signif(PABD_disper_main$`LR Chisq`[2],3))
                                 , paste0("p=",signif(eBD_disper_main$`Pr(>F)`[2],3), ", F(",eBD_disper_main$Df[2],",",eBD_disper_main$Df[4],")=",signif(eBD_disper_main$`F value`[2],3))
                                 , paste0("p=",signif(eBD_disper_main1$`Pr(>F)`[2],3), ", F(",eBD_disper_main1$Df[2],",",eBD_disper_main1$Df[4],")=",signif(eBD_disper_main1$`F value`[2],3))
)

#### PABD and Inhibitory ####

#' 
#' \
#' (3a) Does composition of microbiome influence BD infection risk?\
#' Now, we ask if composition-- specitically, the richness and percent of BD inhibitory bacteria-- 
#' influences infection risk in individuals. First, below, we use just a regular correlation
#' between richness and infection risk

glm_PABD_pinhibRich <- glm(PABD ~ species*p_inhibRich, data=all_p, family=binomial)
PABD_inhibRich_inter <- anova(glm_PABD_pinhibRich, test="Chisq")
PABD_inhibRich_inter
PABD_inhibRich_main <- Anova(glm_PABD_pinhibRich, type=2) #### SIG
PABD_inhibRich_main
Anova(glm(PABD ~ species + p_inhibRich, data=all_p, family=binomial), type=2)
all_p %>%
    ggplot(aes(x=p_inhibRich, y=PABD)) +
    geom_point(aes(col=species), cex=3)

#' Now let's do percent inhibitory of standardized values
glm_PABD_ppinhib <- glm(PABD ~ species*p_percInhib, data=all_p, family=binomial)
PABD_percInhib_inter <- anova(glm_PABD_ppinhib, test="Chisq")
PABD_percInhib_inter
PABD_percInhib_main <- Anova(glm_PABD_ppinhib, type=2)
PABD_percInhib_main
all_p %>%
    ggplot(aes(x=p_percInhib, y=PABD)) +
    geom_point(aes(col=species), cex=3)


#### eBD and Inhibitory ####

#' (3b) Does composition of microbiome influence BD infection intensity?\

#' Version without zeros
lm_eBD_pinhibRich <- lm(eBD_log_infected ~ species*p_inhibRich, data=all_p)
eBD_inhibRich_inter1 <- anova(lm_eBD_pinhibRich)
eBD_inhibRich_inter1
eBD_inhibRich_main1 <- Anova(lm_eBD_pinhibRich, type=2)
eBD_inhibRich_main1
all_p %>%
    ggplot(aes(x=p_inhibRich, y=eBD_log_infected)) +
    geom_point(aes(col=species), cex=3)

#' Version with zeros
lm_eBD_pinhibRich_nozeros <- lm(eBD_log ~ species*p_inhibRich, data=all_p)
eBD_inhibRich_inter <- anova(lm_eBD_pinhibRich_nozeros)
eBD_inhibRich_inter
eBD_inhibRich_main <- Anova(lm_eBD_pinhibRich_nozeros, type=2)
eBD_inhibRich_main
all_p %>%
    ggplot(aes(x=p_inhibRich, y=eBD_log)) +
    geom_point(aes(col=species), cex=3)

#' Now let's do percent inhibitory of standardized values; no zeros
lm_eBD_ppinhib <- lm(eBD_log_infected ~  species*p_percInhib, data=all_p)
eBD_percInhib_inter1 <- anova(lm_eBD_ppinhib)
eBD_percInhib_inter1
eBD_percInhib_main1 <- Anova(lm_eBD_ppinhib, type=2)
eBD_percInhib_inter1
all_p %>%
    ggplot(aes(x=p_percInhib, y=eBD_log_infected)) +
    geom_point(aes(col=species), cex=3) 

# Remove non-infected individuals and re-run
lm_eBD_ppinhib_nozeros <- lm(eBD_log ~  species*p_percInhib, data=all_p)
eBD_percInhib_inter<-anova(lm_eBD_ppinhib_nozeros)
eBD_percInhib_inter
eBD_percInhib_main <- Anova(lm_eBD_ppinhib_nozeros, type=2)
eBD_percInhib_main
all_p %>%
    ggplot(aes(x=p_percInhib, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 


all_stats["Inhibitory ASV Richness",1:4] <- c("Inhibitory ASV Richness"
                                              , paste0("p=",signif(PABD_inhibRich_main$`Pr(>Chisq)`[2],3), " LR Chisq(",PABD_inhibRich_main$Df[2],",N=",nrow(all_p),")=",signif(PABD_inhibRich_main$`LR Chisq`[2],3))
                                              , paste0("p=",signif(eBD_inhibRich_main$`Pr(>F)`[2],3), ", F(",eBD_inhibRich_main$Df[2],",",eBD_inhibRich_main$Df[4],")=",signif(eBD_inhibRich_main$`F value`[2],3))
                                              , paste0("p=",signif(eBD_inhibRich_main1$`Pr(>F)`[2],3), ", F(",eBD_inhibRich_main1$Df[2],",",eBD_inhibRich_main1$Df[4],")=",signif(eBD_inhibRich_main1$`F value`[2],3))
)
all_stats["Proportion Inhibitory",1:4] <- c("Proportion Inhibitory"
                                            , paste0("p=",signif(PABD_percInhib_main$`Pr(>Chisq)`[2],3), " LR Chisq(",PABD_percInhib_main$Df[2],",N=",nrow(all_p),")=",signif(PABD_percInhib_main$`LR Chisq`[2],3))
                                            , paste0("p=",signif(eBD_percInhib_main$`Pr(>F)`[2],3), ", F(",eBD_percInhib_main$Df[2],",",eBD_percInhib_main$Df[4],")=",signif(eBD_percInhib_main$`F value`[2],3))
                                            , paste0("p=",signif(eBD_percInhib_main1$`Pr(>F)`[2],3), ", F(",eBD_percInhib_main1$Df[2],",",eBD_percInhib_main1$Df[4],")=",signif(eBD_percInhib_main1$`F value`[2],3))
)


####Part II: Affect of BD infection on microbiome state####
#' 
#' Part II: Affect of BD infection on microbiome state
#' 
#' 
#### OTU Richness and PABD ####

#' (1a) Does BD infection state affect microbiome diversity?

lm_prich_PABD <- glm(p_rich ~ species*PABD, data=all_p_infected)
rich_PABD_inter <- anova(lm_prich_PABD)
rich_PABD_inter
rich_PABD_main <- Anova(lm_prich_PABD, type=2)
rich_PABD_main
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=p_rich)) +
    geom_boxplot() +
    geom_point(aes(col=species), cex=3, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)

#### OTU Richness and eBD ####
#' 
#' (1b) Does BD infection intensity affect microbiome diversity?

lm_prich_eBD <- lm(p_rich ~ species*eBD_log, data=all_p_infected)
rich_eBD_inter <- anova(lm_prich_eBD)
rich_eBD_inter
rich_eBD_main <-Anova(lm_prich_eBD, type=2)
rich_eBD_main
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_rich)) +
    geom_point(aes(col=species), cex=3)

lm_prich_eBD_nozeros <- lm(p_rich ~ species*eBD_log_infected, data=all_p_infected)
rich_eBD_inter1 <- anova(lm_prich_eBD_nozeros)
rich_eBD_inter1
rich_eBD_main1 <-Anova(lm_prich_eBD_nozeros, type=2)
rich_eBD_main1
all_p_infected %>%
    ggplot(aes(x=eBD_log_infected, y=p_rich)) +
    geom_point(aes(col=species), cex=3)


all_stats["ASV Richness",5:7] <- c( paste0("p=",signif(rich_PABD_main$`Pr(>Chisq)`[2],3), ", LR Chisq(",rich_PABD_main$Df[2],",N=",sum(!is.na(all_p_infected$p_rich)),")=",signif(rich_PABD_main$`LR Chisq`[2],3))
                                             , paste0("p=",signif(rich_eBD_main$`Pr(>F)`[2],3), ", F(",rich_eBD_main$Df[2],",",rich_eBD_main$Df[4],")=",signif(rich_eBD_main$`F value`[2],3))
                                    , paste0("p=",signif(rich_eBD_main1$`Pr(>F)`[2],3), ", F(",rich_eBD_main1$Df[2],",",rich_eBD_main1$Df[4],")=",signif(rich_eBD_main1$`F value`[2],3))
)


#### Instability and PABD ####

#' (2a) Does BD infection state affect microbiome instability?

lm_pdist_PABD <- lm(p_dist ~ species*PABD, data=all_p_infected)
dist_PABD_inter <- anova(lm_pdist_PABD)
dist_PABD_inter
dist_PABD_main <- Anova(lm_pdist_PABD)
dist_PABD_main
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=p_dist)) +
    geom_boxplot() +
    geom_point(aes(color=species), cex=4, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)

#### Dispersion and eBD ####

#' (2b) Does BD infection intensity affect microbiome instability?
#' 

lm_pdist_eBD <- lm(p_dist ~ species*eBD_log, data=all_p_infected)
dist_eBD_inter <- anova(lm_pdist_eBD)
dist_eBD_inter
dist_eBD_main <-Anova(lm_pdist_eBD)
dist_eBD_main
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_dist)) +
    geom_point(aes(color=species), cex=4)

# Remove non-infected individuals and re-run
lm_pdist_eBD_nozeros <- lm(p_dist ~ species*eBD_log_infected, data=all_p_infected)
dist_eBD_inter1 <- anova(lm_pdist_eBD_nozeros)
dist_eBD_inter1
dist_eBD_main1 <- Anova(lm_pdist_eBD_nozeros)
dist_eBD_main1
all_p_infected %>%
    ggplot(aes(x=eBD_log_infected, y=p_dist)) +
    geom_point(aes(color=species), cex=4)

all_stats["Dispersion",5:7] <- c( paste0("p=",signif(dist_PABD_main$`Pr(>F)`[2],3), ", F(",dist_PABD_main$Df[2],",",dist_PABD_main$Df[4],")=",signif(dist_PABD_main$`F value`[2],3))
                                   , paste0("p=",signif(dist_eBD_main$`Pr(>F)`[2],3), ", F(",dist_eBD_main$Df[2],",",dist_eBD_main$Df[4],")=",signif(dist_eBD_main$`F value`[2],3))
                                  , paste0("p=",signif(dist_eBD_main1$`Pr(>F)`[2],3), ", F(",dist_eBD_main1$Df[2],",",dist_eBD_main1$Df[4],")=",signif(dist_eBD_main1$`F value`[2],3))
)


#### Dispersion and PABD ####

#' (2a) Does BD infection state affect microbiome dispersion?
# 

lm_pdisper_PABD <- lm(p_disper ~ species*PABD, data=all_p_infected)
disper_PABD_inter <- anova(lm_pdisper_PABD) ## SIG interaction?
disper_PABD_inter
disper_PABD_main <- Anova(lm_pdisper_PABD)
disper_PABD_main
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=p_disper)) +
    geom_boxplot() +
    geom_point(aes(color=species), cex=4, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)

#### Dispersion and eBD ####

#' (2b) Does BD infection intensity affect microbiome dispersion?
#' 

lm_pdisper_eBD <- lm(p_disper ~ species*eBD_log, data=all_p_infected)
disper_eBD_inter <- anova(lm_pdisper_eBD)
disper_eBD_inter
disper_eBD_main <- Anova(lm_pdisper_eBD)
disper_eBD_main
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_disper)) +
    geom_point(aes(color=species), cex=4)

# Remove non-infected individuals and re-run

lm_pdisper_eBD_nozros <- lm(p_disper ~ species*eBD_log_infected, data=all_p_infected)
disper_eBD_inter1 <- anova(lm_pdisper_eBD_nozros)
disper_eBD_inter1
disper_eBD_main1 <- Anova(lm_pdisper_eBD_nozros)
disper_eBD_main1
all_p_infected %>%
    ggplot(aes(x=eBD_log_infected, y=p_disper)) +
    geom_point(aes(color=species), cex=4)

all_stats["Instability",5:7] <- c( paste0("p=",signif(disper_PABD_main$`Pr(>F)`[2],3), ", F(",disper_PABD_main$Df[2],",",disper_PABD_main$Df[4],")=",signif(disper_PABD_main$`F value`[2],3))
                                   , paste0("p=",signif(disper_eBD_main$`Pr(>F)`[2],3), ", F(",disper_eBD_main$Df[2],",",disper_eBD_main$Df[4],")=",signif(disper_eBD_main$`F value`[2],3))
                                   , paste0("p=",signif(disper_eBD_main1$`Pr(>F)`[2],3), ", F(",disper_eBD_main1$Df[2],",",disper_eBD_main1$Df[4],")=",signif(disper_eBD_main1$`F value`[2],3))
)

#### Inhibitory and PABD ####

#' (3a) Does BD infection state affect microbiome composition?

lm_pinhibRich_PABD <- lm(p_inhibRich ~ species*PABD, data=all_p_infected)
inhibRich_PABD_inter <- anova(lm_pinhibRich_PABD)
inhibRich_PABD_inter
inhibRich_PABD_main <-Anova(lm_pinhibRich_PABD)
inhibRich_PABD_main
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=p_inhibRich)) +
    geom_boxplot() +
    geom_point(aes(col=species), cex=3, position=position_jitter(width=0.15, height=0)) +
    facet_wrap(~species, nrow=1)

lm_ppercInhib_PABD <- lm(p_percInhib ~ species*PABD, data=all_p_infected)
percInhib_PABD_inter <- anova(lm_ppercInhib_PABD)
percInhib_PABD_inter
percInhib_PABD_main <- Anova(lm_ppercInhib_PABD)
percInhib_PABD_main
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=p_percInhib)) +
    geom_boxplot() +
    geom_point(aes(col=species), cex=3, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)


#### Inhibitory and eBD ####

#' 
#' (3b) Does BD infection intensity affect microbiome composition?

lm_pinhibRich_eBD <- lm(p_inhibRich ~ species*eBD_log, data=all_p_infected)
inhibRich_eBD_inter <-anova(lm_pinhibRich_eBD)
inhibRich_eBD_inter
inhibRich_eBD_main <- Anova(lm_pinhibRich_eBD, type="II")
inhibRich_eBD_main
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_inhibRich)) +
    geom_point(aes(col=species), cex=3) 

# Re-run with non-zeros
lm_pinhibRich_eBD_nozeros <- lm(p_inhibRich ~ species*eBD_log_infected, data=all_p_infected)
inhibRich_eBD_inter1 <- anova(lm_pinhibRich_eBD_nozeros)
inhibRich_eBD_inter1
inhibRich_eBD_main1 <- Anova(lm_pinhibRich_eBD_nozeros, type=2)
inhibRich_eBD_main1
all_p_infected %>%
    ggplot(aes(x=eBD_log_infected, y=p_inhibRich)) +
    geom_point(aes(col=species), cex=3) 

# Now do percent inhibitory
lm_ppercInhib_eBD <- lm(p_percInhib ~  species*eBD_log, data=all_p_infected)
percInhib_eBD_inter <- anova(lm_ppercInhib_eBD)
percInhib_eBD_inter
percInhib_eBD_main <- Anova(lm_ppercInhib_eBD, type = 2)
percInhib_eBD_main
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_percInhib)) +
    geom_point(aes(col=species), cex=3) +
    facet_wrap(~species, nrow=1)

# Re-run with non-zeros
lm_ppercInhib_eBD_nozeros <- lm(p_percInhib ~  species*eBD_log_infected, data=all_p_infected)
percInhib_eBD_inter1 <- anova(lm_ppercInhib_eBD_nozeros)
percInhib_eBD_inter1
percInhib_eBD_main1 <- Anova(lm_ppercInhib_eBD_nozeros, type = 2)
percInhib_eBD_main1
all_p_infected %>%
    ggplot(aes(x=eBD_log_infected, y=p_percInhib)) +
    geom_point(aes(col=species), cex=3) +
    facet_wrap(~species, nrow=1)



all_stats["Inhibitory ASV Richness",5:7] <- c( paste0("p=",signif(inhibRich_PABD_main$`Pr(>F)`[2],3), ", F(",inhibRich_PABD_main$Df[2],",",inhibRich_PABD_main$Df[4],")=",signif(inhibRich_PABD_main$`F value`[2],3))
                                               , paste0("p=",signif(inhibRich_eBD_main$`Pr(>F)`[2],3), ", F(",inhibRich_eBD_main$Df[2],",",inhibRich_eBD_main$Df[4],")=",signif(inhibRich_eBD_main$`F value`[2],3))
                                               , paste0("p=",signif(inhibRich_eBD_main1$`Pr(>F)`[2],3), ", F(",inhibRich_eBD_main1$Df[2],",",inhibRich_eBD_main1$Df[4],")=",signif(inhibRich_eBD_main1$`F value`[2],3))
)
all_stats["Proportion Inhibitory",5:7] <- c( paste0("p=",signif(percInhib_PABD_main$`Pr(>F)`[2],3), ", F(",percInhib_PABD_main$Df[2],",",percInhib_PABD_main$Df[4],")=",signif(percInhib_PABD_main$`F value`[2],3))
                                             , paste0("p=",signif(percInhib_eBD_main$`Pr(>F)`[2],3), ", F(",percInhib_eBD_main$Df[2],",",percInhib_eBD_main$Df[4],")=",signif(percInhib_eBD_main$`F value`[2],3))
                                             , paste0("p=",signif(percInhib_eBD_main1$`Pr(>F)`[2],3), ", F(",percInhib_eBD_main1$Df[2],",",percInhib_eBD_main1$Df[4],")=",signif(percInhib_eBD_main1$`F value`[2],3))
)


all_stats
write.csv(all_stats, "all_stats_driversofinfection.csv")

#### FOLLOW UP ####
# get all mf values
mf_all_without_init_infect <- mf.rare %>% 
    filter(SampleID %in%c(mf_con_without_init_infect$SampleID, mf_treat_without_init_infect$SampleID) )

# combine all_p?
# Is OTU richness and inhibitory bacterial richness related?
all_p_withcon <- all_p_withcon %>%
    separate(toadID, into=c("species","indiv"), remove = FALSE)
all_p_withcon %>%
    ggplot(aes(x=p_logRich, y=p_inhibRich)) +
    geom_point(aes(col=species), cex=3) +
    geom_smooth(method="lm")
cor.test(all_p$p_logRich, all_p$p_inhibRich, method="pearson")

anova(lm(p_inhibRich ~ p_logRich*species, data=all_p_withcon))
Anova(lm(p_inhibRich ~ p_logRich*species, data=all_p_withcon), type="II")

# No, it's not-- it means it's decoupled

# Is percent inhibitory and inhibitory bacterial richness related?

all_p_withcon %>%
    ggplot(aes(x=p_percInhib, y=p_inhibRich)) +
    geom_point(aes(col=species), cex=3) +
    geom_smooth(method="lm")
cor.test(all_p$p_percInhib, all_p$p_inhibRich, method="pearson")

anova(lm(p_inhibRich ~ p_percInhib*species, data=all_p))
Anova(lm(p_inhibRich ~ p_percInhib*species, data=all_p), type="II")

# Is inhibitory richness correlated with dispersion, since those are two driving factors of bd risk?
all_p_withcon %>%
    ggplot(aes(x=p_disper, y=p_inhibRich)) +
    geom_point(aes(col=species), cex=3) +
    geom_smooth(method="lm")
cor.test(all_p$p_disper, all_p$p_inhibRich, method="pearson")
anova(lm(p_disper ~ p_inhibRich*species, data=all_p))
Anova(lm(p_disper ~ p_inhibRich*species, data=all_p), type="II")

anova(lm(p_disper ~ p_inhibRich*species, data=all_p_withcon))
Anova(lm(p_disper ~ p_inhibRich*species, data=all_p_withcon), type="II")

# Is diversity/richness correlated with stability?
all_p_withcon %>%
    ggplot(aes(x=p_logRich, y=p_dist)) +
    geom_point(aes(col=species), cex=3) +
    geom_smooth(method="lm")
anova(lm(p_dist ~ p_logRich*species, data=all_p))
Anova(lm(p_dist ~ p_logRich*species, data=all_p), type="II")

#### NMDS of inhibitory OTUs ####

## CONTROLS
rownames(otu.inhibOnly.con) <- otu.inhibOnly.con$OTUID
otu.inhibOnly.con_fornmds <- otu.inhibOnly.con[,-match("OTUID",colnames(otu.inhibOnly.con))]
dm_bray_con_inhib <- vegdist(t(otu.inhibOnly.con_fornmds), method = "bray")
set.seed(01238)
nmds_con_inhib <- isoMDS(dm_bray_con_inhib, k=2)

mf_con_without_init_infect$NMDS1_inhib <- NA
mf_con_without_init_infect$NMDS2_inhib <- NA
mf_con_without_init_infect[match(rownames(nmds_con_inhib$points), mf_con_without_init_infect$SampleID),c("NMDS1_inhib","NMDS2_inhib")] <- nmds_con_inhib$points

mf_con_without_init_infect %>%
    # filter(species=="Osse") %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=species,alpha=time))
adonis(dm_bray_con_inhib ~ species*time, data=mf_con_without_init_infect)


## TREAT
rownames(otu.inhibOnly.treat) <- otu.inhibOnly.treat$OTUID
otu.inhibOnly.treat_fornmds <- otu.inhibOnly.treat[,-match("OTUID",colnames(otu.inhibOnly.treat))]
dm_bray_treat_inhib <- vegdist(t(otu.inhibOnly.treat_fornmds), method = "bray")
set.seed(01238)
nmds_treat_inhib <- isoMDS(dm_bray_treat_inhib, k=2)

mf_treat_without_init_infect$NMDS1_inhib <- NA
mf_treat_without_init_infect$NMDS2_inhib <- NA
mf_treat_without_init_infect[match(rownames(nmds_treat_inhib$points), mf_treat_without_init_infect$SampleID),c("NMDS1_inhib","NMDS2_inhib")] <- nmds_treat_inhib$points

mf_treat_without_init_infect %>%
    # filter(species=="Osse") %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=species,alpha=time))
adonis(dm_bray_treat_inhib ~ species*time, data=mf_treat_without_init_infect)

## TREAT-- POST exposure only
mf_treat_postexp <- mf_treat_without_init_infect %>%
    filter(prepost=="Pos")
otu.inhibOnly.treat_POS <- otu.inhibOnly.treat[,match(mf_treat_postexp$SampleID, colnames(otu.inhibOnly.treat))]
# rownames(otu.inhibOnly.treat_POS) <- otu.inhibOnly.treat_POS$OTUID
# otu.inhibOnly.treat_POS_fornmds <- otu.inhibOnly.treat_POS[,-match("OTUID",colnames(otu.inhibOnly.treat_POS))]
dm_bray_treat_inhib_POS <- vegdist(t(otu.inhibOnly.treat_POS), method = "bray")
set.seed(01238)
nmds_treat_inhib_POS <- isoMDS(dm_bray_treat_inhib_POS, k=2)

mf_treat_postexp$NMDS1_inhib <- NA
mf_treat_postexp$NMDS2_inhib <- NA
mf_treat_postexp[match(rownames(nmds_treat_inhib_POS$points), mf_treat_postexp$SampleID),c("NMDS1_inhib","NMDS2_inhib")] <- nmds_treat_inhib_POS$points

mf_treat_postexp %>%
    # filter(species=="Osse") %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=species,alpha=time))

## EFFECT OF INFECTION ON MICROBIOME??
adonis(dm_bray_treat_inhib_POS ~ species*time*PABD, data=mf_treat_postexp)
adonis(dm_bray_treat_inhib_POS ~ species*time*eBD_log, data=mf_treat_postexp)


## COMBINED
otu_inhibOnly_combo <- otu.inhibOnly.con %>%
    left_join(otu.inhibOnly.treat,by="OTUID")
mf_combined <- rbind(mf_con_without_init_infect, mf_treat_without_init_infect)
#make sure in same order
mf_combined <- mf_combined[match(colnames(otu_inhibOnly_combo), mf_combined$SampleID),]

rownames(otu_inhibOnly_combo) <- otu_inhibOnly_combo$OTUID
otu_inhibOnly_combo_fornmds <- otu_inhibOnly_combo[,-match("OTUID",colnames(otu_inhibOnly_combo))]
dm_bray_combo_inhib <- vegdist(t(otu_inhibOnly_combo_fornmds), method = "bray")
set.seed(01238)
nmds_combo_inhib <- isoMDS(dm_bray_combo_inhib, k=2)

mf_combined$NMDS1_inhib <- NA
mf_combined$NMDS2_inhib <- NA
mf_combined[match(rownames(nmds_combo_inhib$points), mf_combined$SampleID),c("NMDS1_inhib","NMDS2_inhib")] <- nmds_combo_inhib$points

mf_combined %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=species,alpha=time, pch=BD_infected))+
    scale_shape_manual(values=c(`y`=16, `n`=21))

## EFFECT OF EXPOSURE ON INHIB?
adonis(dm_bray_combo_inhib ~ species*BD_infected*prepost, data=mf_combined)

mf_combined %>%
    filter(species=="Anbo") %>%
    mutate(Treatment=BD_infected, Infected=PABD) %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=Treatment,size=(time), pch=Infected))+
    scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21))
mf_combined %>%
    filter(species=="Anma") %>%
    mutate(Treatment=BD_infected, Infected=PABD) %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=Treatment,size=(time), pch=Infected))+
    scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21))
mf_combined %>%
    filter(species=="Lica") %>%
    mutate(Treatment=BD_infected, Infected=PABD) %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=Treatment,size=(time), pch=Infected))+
    scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21))
mf_combined %>%
    filter(species=="Lipi") %>%
    mutate(Treatment=BD_infected, Infected=PABD) %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=Treatment,size=(time), pch=Infected))+
    scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21))
mf_combined %>%
    filter(species=="Osse") %>%
    mutate(Treatment=BD_infected, Infected=PABD) %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=Treatment,size=(time), pch=Infected))+
    scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21))



#### Inhibitory ASVs: who is important? ####

# remove things that are absent in both
zero_genera <- rbind(mf_con_with_inhibOTUs, mf_treat_with_inhibOTUs) %>%
    group_by(G) %>%
    summarize(total=sum(reads)) %>%
    filter(total==0) %>%
    pull(G)

# Collapse by genus
# Combined version
mf_treat_with_inhibOTUs$TreatmentGroup <- "Treatment"
mf_con_with_inhibOTUs$TreatmentGroup <- "Control"

mf_combo_with_inhibOTUs_G <- rbind(mf_treat_with_inhibOTUs,mf_con_with_inhibOTUs) %>%
    filter(!(G %in% zero_genera)) %>%
    group_by(SampleID, species, time, toadID, eBD_log, PABD, prepost,TreatmentGroup, OTUID, K, P, C, O, `F`, G) %>%
    summarize(genus_proportion = sum(reads))   %>%
    group_by(G) %>%
    mutate(max_g_prop = max(genus_proportion)) %>%
    ungroup() %>%
    mutate(standardized_g = ifelse(max_g_prop>0,genus_proportion/max_g_prop,0)) 

mf_combo_with_inhibOTUs_G %>%
    ggplot(aes(x=time, y=genus_proportion, col=species), cex=0.1) +
    geom_point(alpha=0.3) +
    geom_vline(aes(xintercept=5.5), col="grey", lty=2) +
    facet_grid(G~TreatmentGroup, scales = "free")

#### Effect of time, exposure, and BD infection ####

#' ### Effect of time, exposure, and BD infection ###
#' A problem with our dataset is that there are several levels of factors that can confound
#' the effect of BD infection. First, there can be an effect of TIME, which is weakly correlated with BD infection
cor.test(mf_treat_without_init_infect$time, mf_treat_without_init_infect$eBD_log, method="kendall")
#' Second, time is nested within BD exposure. This is obvious because ONLY latter time points have 
#' BD exposure. 
#' \
#' Therefore, in order to correctly determine wheter BD infection itself has an effect on the microbiome, we must control for time and
#' BD exposure.
#' To control for time, we can compare the effect of time in the controls against the effect of time in the treatments (which we do above), and
#' then assess whether there is a significant interaction of treatment*time or treatment*exposure. Here, we assume no interact effects between species and other factors.
#' The interaction between treatment and time will tell us whether there differences across time between treatments.

Anova(lm(logRich ~ species+time*BD_infected, data=mf_all_without_init_infect), type=3)
Anova(glm(distance_bray_curtis ~ species+time*BD_infected, data=mf_all_without_init_infect, family="quasibinomial"), type=3)
Anova(lm(log(disper_bray_curtis) ~ species+time*BD_infected, data=mf_all_without_init_infect), type=3)
Anova(glm(inhibRich ~ species+time*BD_infected, data=mf_all_without_init_infect, family="poisson"), type=3)
Anova(lm(percInhib ~ species+time*BD_infected, data=mf_all_without_init_infect), type=3)

Anova(lm(logRich ~ species+prepost*BD_infected, data=mf_all_without_init_infect), type=3)
Anova(glm(distance_bray_curtis ~ species+prepost*BD_infected, data=mf_all_without_init_infect, family="quasibinomial"), type=3)
Anova(lm(log(disper_bray_curtis) ~ species+prepost*BD_infected, data=mf_all_without_init_infect), type=3)
Anova(glm(inhibRich ~ species+prepost*BD_infected, data=mf_all_without_init_infect, family="poisson"), type=3)
Anova(lm(percInhib ~ species+prepost*BD_infected, data=mf_all_without_init_infect), type=3)


### BONUS: affect of interaction treatment and PABD
# For some reason, I can't get betareg to work so I'm using a logistic transformation instead.
anova(lm(qlogis(y.transf.betareg(percInhib)) ~ species +BD_infected/PABD, data=mf_all_without_init_infect))

pinhib_treat_main_glm_sponly<- betareg(y.transf.betareg(percInhib) ~ species, data=mf_treat_without_init_infect)
pinhib_treat_main_glm_exposure<- betareg(y.transf.betareg(percInhib) ~ species + prepost, data=mf_treat_without_init_infect)

pinhib_con_main_glm_sponly<- betareg(y.transf.betareg(percInhib) ~ species, data=mf_con_without_init_infect)
pinhib_con_main_glm_exposure<- betareg(y.transf.betareg(percInhib) ~ species + prepost, data=mf_con_without_init_infect)

pinhib_con_main_exposure <- lrtest(pinhib_con_main_glm_sponly,pinhib_con_main_glm_exposure) # plus prepost
pinhib_con_main_exposure
pinhib_treat_main_exposure <- lrtest(pinhib_treat_main_glm_sponly,pinhib_treat_main_glm_exposure) # plus prepost
pinhib_treat_main_exposure

### BONUS: affect of interaction treatment and BD exposure on dispersion
# For some reason, I can't get betareg to work so I'm using a logistic transformation instead.
# hist(log(mf_all_without_init_infect$disper_bray_curtis))
anova(lm(log((disper_bray_curtis)) ~ species*BD_infected*prepost, data=mf_all_without_init_infect))

#' It looks like only percent inhibitory has a significant interaction, which suggests that the treatment is behaving differently
#' over time than the control in terms of proportion of inhibitory bacteria. However, is it BD infection that is driving this difference
#' over time? Or BD exposure? To test this, we can assess the interaction between BD infection and exposure in a model for only treatment individuals.

anova(lm(logRich ~ species+prepost/PABD, data=mf_treat_without_init_infect))
Anova(glm(distance_bray_curtis ~ species+prepost/PABD, data=mf_treat_without_init_infect, family="quasibinomial"), type=2)
anova(lm(log(disper_bray_curtis) ~ species+prepost/PABD, data=mf_treat_without_init_infect))
Anova(glm(inhibRich ~ species+prepost/PABD, data=mf_treat_without_init_infect, family="poisson"), type=2)
anova(lm(percInhib ~ species+prepost/PABD, data=mf_treat_without_init_infect))

anova(lm(logRich ~ species+prepost/eBD_log, data=mf_treat_without_init_infect))
Anova(glm(distance_bray_curtis ~ species+prepost/eBD_log, data=mf_treat_without_init_infect, family="quasibinomial"), type=2)
anova(lm(log(disper_bray_curtis) ~ species+prepost/eBD_log, data=mf_treat_without_init_infect))
Anova(glm(inhibRich ~ species+prepost/eBD_log, data=mf_treat_without_init_infect, family="poisson"), type=2)
anova(lm(qlogis(y.transf.betareg(percInhib)) ~ species+prepost/eBD_log, data=mf_treat_without_init_infect))

#' What we find is that in dispersion and percent inhibitory, BD infection seems to have an effect on TOP of BD exposure (in terms
#' of both presence/absence infection, and infection load). Thus, we can infer that infection really does have an effect on percent inhibitory.
#' \ Side note: curiously, dispersion apparently differs here, even though it is not detected in our other models. 
#' 
#' Now, our next question is: is it the presence/absence of BD or the load that causes percent inhibitory to increase? We can test this by seeing if
#' intensity explains an additional portion of the model when compared to just presence/absence

Anova(lm(logRich ~ species+PABD/eBD_log, data=all_p_infected), type=3)
Anova(lm(distance_bray_curtis ~ species+PABD/eBD_log, data=all_p_infected), type=3)
Anova(lm(disper_bray_curtis ~ species+PABD/eBD_log, data=all_p_infected), type=3)
Anova(lm(inhibRich ~ species+PABD/eBD_log, data=all_p_infected), type=3)
Anova(lm(percInhib ~ species+PABD/eBD_log, data=all_p_infected), type=3)


Anova(lm(logRich ~ species*prepost, data=mf_con_without_init_infect), type=3)
Anova(lm(distance_bray_curtis ~ species*prepost, data=mf_con_without_init_infect), type=2)
Anova(lm(disper_bray_curtis ~ species*prepost, data=mf_con_without_init_infect), type=2)
Anova(lm(inhibRich ~ species*prepost, data=mf_con_without_init_infect), type=3)
Anova(lm(percInhib ~ species*prepost, data=mf_con_without_init_infect), type=3)


