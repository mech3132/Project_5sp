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
library(vegan) # for permanova
library(gridExtra)
library(betareg) # for beta distr
library(lmtest) # for beta analysis
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

# add a species column and PABD column
all_p <- all_p %>%
    mutate(PABD=ifelse(infect>0,1,0), infect = log(infect+1)) %>%
    rename(eBD_log=infect) %>%
    separate(toadID, into=c("species","indiv"), remove=FALSE)
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
    mutate(metric="logRich") %>%
    rename(value=logRich)
temp2 <- mf_con_without_init_infect %>%
    dplyr::select(species, inhibRich) %>%
    mutate(metric="inhibRich")%>%
    rename(value=inhibRich)
temp3 <- mf_con_without_init_infect %>%
    dplyr::select(species, percInhib) %>%
    mutate(metric="percInhib")%>%
    rename(value=percInhib)
temp4 <- mf_con_without_init_infect %>%
    dplyr::select(species, distance_bray_curtis) %>%
    mutate(metric="distance_bray_curtis")%>%
    rename(value=distance_bray_curtis)

gg_all <- rbind(temp1,temp2,temp3,temp4) %>%
    ggplot(aes(x=species, y=value)) +
    geom_boxplot() +
    geom_point(aes(col=species), position = position_jitter(width=0.1, height=0))+
    facet_grid(metric~., scales = "free") 
lay <- rbind(c(1,2),
             c(3,2))
#+ fig.height=10, fig.width=10
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
centroid_con_interaction_lm <- lm(log(distance.to.centroid) ~ species*time, data=mf_con_without_init_infect)
centroid_con_interaction <- anova(centroid_con_interaction_lm)
centroid_con_interaction
# Type II ANOVA with no interaction
centroid_con_main_lm <- lm(log(distance.to.centroid) ~ species + time, data=mf_con_without_init_infect)
centroid_con_main <- Anova(centroid_con_main_lm, type = 2)
centroid_con_main

#' Is there an effect of species and time on treatment??
# Type I ANOVA (to check for interaction) (AB | A,B)
centroid_treat_interaction_lm <- lm(log(distance.to.centroid) ~ species*time, data=mf_treat_without_init_infect)
centroid_treat_interaction <- anova(centroid_treat_interaction_lm)
centroid_treat_interaction
# Type II ANOVA with no interaction
centroid_treat_main_lm <- lm(log(distance.to.centroid) ~ species + time, data=mf_treat_without_init_infect)
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

#### PABD and OTU Richness ####

#' ### (1a) Does overall diversity of microbiome influence BD infection rate?\
#' The first thing we would like to know is whether microbiome richness of an individual influences 
#' its risk of becoming infected by BD. The most simple way to look at this would be
#' to plot OTU richness VS presence/absence of BD
#' Below, we fitted normal and lognormal distributions, respectively, to diversity (shannon) and
#' otu richness to individuals prior to BD infection. Now, we fit a binomial general
#' linearized model to see if there is a relationship between diversity and infection rate.
#' \

glm_PABD_prich <- glm(PABD ~ species*p_rich, data=all_p, family=binomial(link="logit"))
anova(glm_PABD_prich, test="Chisq") # test for interaction
Anova(glm_PABD_prich, type=2) # test for main effects
all_p %>%
    ggplot(aes(x=p_rich, y=PABD)) +
    geom_point(aes(col=species), cex=3)  

#' If anything, it looks like increased diversity and richness might increase infection risk

#'\

#### eBD and OTU Richness ####

#' (1b) Does overall diversity of microbiome influence BD infection intensity?\
#' \
#' The next thing we would like to know is if richness of the microbiome
#' influences infection intensity.
#' \

#' Now let's do richness
lm_eBD_prich <- lm(eBD_log ~ species*p_rich, data=all_p)
anova(lm_eBD_prich)
Anova(lm_eBD_prich, type=2)
all_p %>%
    ggplot(aes(x=p_rich, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 

#' Try a version where we remove zeros so that we do not have zero-inflated data
lm_eBD_prich_nozeros <- lm(eBD_log_infected ~ species*p_rich, data=all_p)
anova(lm_eBD_prich_nozeros)
Anova(lm_eBD_prich_nozeros, type=2)
all_p %>%
    ggplot(aes(x=p_rich, y=eBD_log_infected)) +
    geom_point(aes(col=species), cex=3) 

#### PABD and Instability  ####

#' (2a) Does instability of microbiome influence BD infection rate?\
#' Here we look at average distance travelled (bray-curtis) between samples
#' prior to being infected. We see if it is correlated to infection risk.\

glm_PABD_pbc <- glm(PABD ~ species*p_mu, data=all_p, family=binomial)
anova(glm_PABD_pbc, test="Chisq")
Anova(glm_PABD_pbc, type=2)
all_p %>%
    ggplot(aes(x=p_mu, y=PABD)) +
    geom_point(aes(col=species), cex=3) 


#### eBD and Instability ####
 
#' \
#' (2b) Does instability of microbiome influence BD infection intensity?\

lm_BD_pbc <- lm(eBD_log ~ species*p_mu, data=all_p)
anova(lm_BD_pbc)
Anova(lm_BD_pbc, type=2)
all_p %>%
    ggplot(aes(x=p_mu, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 

#' What if we remove those not infected?
lm_BD_pbc_nozeros <- lm(eBD_log_infected ~ species*p_mu, data=all_p)
anova(lm_BD_pbc_nozeros)
Anova(lm_BD_pbc_nozeros, type=2)
all_p %>%
    ggplot(aes(x=p_mu, y=eBD_log_infected)) +
    geom_point(aes(col=species), cex=3) 

#### PABD and Dispersion ####

#' (2a) Does dispersion of microbiome influence BD infection rate?\
#' Here we look at average distance to centroid (bray-curtis) between samples
#' prior to being infected at same time point. We see if it is correlated to infection risk.

glm_PABD_pdist <- glm(PABD ~ species*p_distmu, data=all_p, family=binomial)
anova(glm_PABD_pdist, test="Chisq")
Anova(glm_PABD_pdist, type=2)
all_p %>%
    ggplot(aes(x=p_distmu, y=PABD)) +
    geom_point(aes(col=species), cex=3) 


#### eBD and Dispersion ####

#' 
#' \
#' (2b) Does dispersion of microbiome influence BD infection intensity?\

lm_BD_pdist <- lm(eBD_log ~ species*p_distmu, data=all_p)
anova(lm_BD_pdist)
Anova(lm_BD_pdist, type=2) ## SIG
all_p %>%
    ggplot(aes(x=p_distmu, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 

#' Let's try removing zeros
lm_BD_pdist_nozeros <- lm(eBD_log_infected ~ species*p_distmu, data=all_p)
anova(lm_BD_pdist_nozeros)
Anova(lm_BD_pdist_nozeros, type=2) ## SIG
all_p %>%
    ggplot(aes(x=p_distmu, y=eBD_log_infected)) +
    geom_point(aes(col=species), cex=3) 

#### PABD and Inhibitory ####

#' 
#' \
#' (3a) Does composition of microbiome influence BD infection risk?\
#' Now, we ask if composition-- specitically, the richness and percent of BD inhibitory bacteria-- 
#' influences infection risk in individuals. First, below, we use just a regular correlation
#' between richness and infection risk

glm_PABD_pinhibRich <- glm(PABD ~ species*p_inhibRich, data=all_p, family=binomial)
anova(glm_PABD_pinhibRich, test="Chisq")
Anova(glm_PABD_pinhibRich, type=2) #### SIG
Anova(glm(PABD ~ species + p_inhibRich, data=all_p, family=binomial), type=2)
all_p %>%
    ggplot(aes(x=p_inhibRich, y=PABD)) +
    geom_point(aes(col=species), cex=3)

#' Now let's do percent inhibitory of standardized values
glm_PABD_ppinhib <- glm(PABD ~ species*p_pinhib, data=all_p, family=binomial)
anova(glm_PABD_ppinhib, test="Chisq")
Anova(glm_PABD_ppinhib, type=2)
all_p %>%
    ggplot(aes(x=p_pinhib, y=PABD)) +
    geom_point(aes(col=species), cex=3)

#### eBD and Inhibitory ####

#' (3b) Does composition of microbiome influence BD infection intensity?\


lm_eBD_pinhibRich <- lm(eBD_log ~ species*p_inhibRich, data=all_p)
anova(lm_eBD_pinhibRich)
Anova(lm_eBD_pinhibRich, type=2)
all_p %>%
    ggplot(aes(x=p_inhibRich, y=eBD_log)) +
    geom_point(aes(col=species), cex=3)

# Remove non-infected individuals and re-run
lm_eBD_pinhibRich_nozeros <- lm(eBD_log_infected ~ species*p_inhibRich, data=all_p)
anova(lm_eBD_pinhibRich_nozeros)
Anova(lm_eBD_pinhibRich_nozeros, type=2)
all_p %>%
    ggplot(aes(x=p_inhibRich, y=eBD_log_infected)) +
    geom_point(aes(col=species), cex=3)

#' Now let's do percent inhibitory of standardized values
lm_eBD_ppinhib <- lm(eBD_log ~  species*p_pinhib, data=all_p)
anova(lm_eBD_ppinhib)
Anova(lm_eBD_ppinhib, type=2)
all_p %>%
    ggplot(aes(x=p_pinhib, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 

# Remove non-infected individuals and re-run
lm_eBD_ppinhib_nozeros <- lm(eBD_log_infected ~  species*p_pinhib, data=all_p)
anova(lm_eBD_ppinhib_nozeros)
Anova(lm_eBD_ppinhib_nozeros, type=2)
all_p %>%
    ggplot(aes(x=p_pinhib, y=eBD_log_infected)) +
    geom_point(aes(col=species), cex=3) 

####Part II: Affect of BD infection on microbiome state####
#' 
#' Part II: Affect of BD infection on microbiome state
#' 
#' 
#### OTU Richness and PABD ####

#' (1a) Does BD infection state affect microbiome diversity?

lm_prich_PABD <- lm(p_rich ~ species*PABD, data=all_p_infected)
anova(lm_prich_PABD)
Anova(lm_prich_PABD, type=2)
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
anova(lm_prich_eBD)
Anova(lm_prich_eBD, type=2)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_rich)) +
    geom_point(aes(col=species), cex=3)

lm_prich_eBD_nozeros <- lm(p_rich ~ species*eBD_log_infected, data=all_p_infected)
anova(lm_prich_eBD_nozeros)
Anova(lm_prich_eBD_nozeros, type=2)
all_p_infected %>%
    ggplot(aes(x=eBD_log_infected, y=p_rich)) +
    geom_point(aes(col=species), cex=3)

#### Instability and PABD ####

#' (2a) Does BD infection state affect microbiome instability?

lm_pbc_PABD <- lm(p_BC ~ species*PABD, data=all_p_infected)
anova(lm_pbc_PABD)
Anova(lm_pbc_PABD)
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=p_BC)) +
    geom_boxplot() +
    geom_point(aes(color=species), cex=4, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)

#### Dispersion and eBD ####

#' (2b) Does BD infection intensity affect microbiome instability?
#' 

lm_pbc_eBD <- lm(p_BC ~ species*eBD_log, data=all_p_infected)
anova(lm_pbc_eBD)
Anova(lm_pbc_eBD)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_BC)) +
    geom_point(aes(color=species), cex=4)

# Remove non-infected individuals and re-run
lm_pbc_eBD_nozeros <- lm(p_BC ~ species*eBD_log_infected, data=all_p_infected)
anova(lm_pbc_eBD_nozeros)
Anova(lm_pbc_eBD_nozeros)
all_p_infected %>%
    ggplot(aes(x=eBD_log_infected, y=p_BC)) +
    geom_point(aes(color=species), cex=4)

#### Dispersion and PABD ####

#' (2a) Does BD infection state affect microbiome dispersion?
# 

lm_pdist_PABD <- lm(p_BCdist ~ species*PABD, data=all_p_infected)
anova(lm_pdist_PABD) ## SIG
Anova(lm_pdist_PABD)
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=p_BCdist)) +
    geom_boxplot() +
    geom_point(aes(color=species), cex=4, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)

#### Dispersion and eBD ####

#' (2b) Does BD infection intensity affect microbiome dispersion?
#' 

lm_pbcdist_eBD <- lm(p_BCdist ~ species*eBD_log, data=all_p_infected)
anova(lm_pbcdist_eBD)
Anova(lm_pbcdist_eBD)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_BCdist)) +
    geom_point(aes(color=species), cex=4)

# Remove non-infected individuals and re-run

lm_pbcdist_eBD_nozros <- lm(p_BCdist ~ species*eBD_log_infected, data=all_p_infected)
anova(lm_pbcdist_eBD_nozros)
Anova(lm_pbcdist_eBD_nozros)
all_p_infected %>%
    ggplot(aes(x=eBD_log_infected, y=p_BCdist)) +
    geom_point(aes(color=species), cex=4)

#### Inhibitory and PABD ####

#' (3a) Does BD infection state affect microbiome composition?

lm_pinhibRich_PABD <- lm(p_inhibRich ~ species*PABD, data=all_p_infected)
anova(lm_pinhibRich_PABD)
Anova(lm_pinhibRich_PABD)
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=p_inhibRich)) +
    geom_boxplot() +
    geom_point(aes(col=species), cex=3, position=position_jitter(width=0.15, height=0)) +
    facet_wrap(~species, nrow=1)

lm_ppercInhib_PABD <- lm(p_percInhib ~ species*PABD, data=all_p_infected)
anova(lm_ppercInhib_PABD)
Anova(lm_ppercInhib_PABD)
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
Anova(lm_pinhibRich_eBD)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_inhibRich)) +
    geom_point(aes(col=species), cex=3) 

# Re-run with non-zeros
lm_pinhibRich_eBD_nozeros <- lm(p_inhibRich ~ species*eBD_log_infected, data=all_p_infected)
Anova(lm_pinhibRich_eBD_nozeros)
Anova(lm_pinhibRich_eBD_nozeros, type=2)
all_p_infected %>%
    ggplot(aes(x=eBD_log_infected, y=p_inhibRich)) +
    geom_point(aes(col=species), cex=3) 

# No do percent inhibitory
lm_ppercInhib_eBD <- lm(p_percInhib ~  species*eBD_log, data=all_p_infected)
anova(lm_ppercInhib_eBD)
Anova(lm_ppercInhib_eBD, type = 2)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_percInhib)) +
    geom_point(aes(col=species), cex=3) +
    facet_wrap(~species, nrow=1)

# Re-run with non-zeros
lm_ppercInhib_eBD_nozeros <- lm(p_percInhib ~  species*eBD_log_infected, data=all_p_infected)
anova(lm_ppercInhib_eBD_nozeros)
Anova(lm_ppercInhib_eBD_nozeros, type = 2)
all_p_infected %>%
    ggplot(aes(x=eBD_log_infected, y=p_percInhib)) +
    geom_point(aes(col=species), cex=3) +
    facet_wrap(~species, nrow=1)

#### FOLLOW UP ####


# Is OTU richness and inhibitory bacterial richness related?

all_p %>%
    ggplot(aes(x=p_rich, y=p_inhibRich)) +
    geom_point(aes(col=species), cex=3) +
    geom_smooth(method="lm")
anova(lm(p_inhibRich ~ p_rich, data=all_p))

# No, it's not-- it means it's decoupled


# Is percent inhibitory and inhibitory bacterial richness related?

all_p %>%
    ggplot(aes(x=p_pinhib, y=p_inhibRich)) +
    geom_point(aes(col=species), cex=3) +
    geom_smooth(method="lm")
anova(lm(p_inhibRich ~ p_pinhib, data=all_p))

# No, it's not-- it means it's decoupled


