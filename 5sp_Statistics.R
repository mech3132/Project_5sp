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
# Mapping files
load("mf_con_without_init_infect.RData")
load("mf_treat_without_init_infect.RData")
load("mf.rare.RData")

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

#### PART I ####

#' Part I: Microbiome state and effect on infection risk and intensity\
#' \
#' (1a) Does overall diversity of microbiome influence BD infection rate?\
#' The first thing we would like to know is whether microbiome richness of an individual influences 
#' its risk of becoming infected by BD. The most simple way to look at this would be
#' to plot diversity VS presence/absence of BD and OTU richness VS presence/absence of BD
#' Below, we fitted normal and lognormal distributions, respectively, to diversity (shannon) and
#' otu richness to individuals prior to BD infection. Now, we fit a binomial general
#' linearized model to see if there is a relationship between diversity and infection rate.
#' \

glm_PABD_shan <- glm(PABD ~ species*exp_shan, data=all_p, family=binomial(link="logit"))
Anova(glm_PABD_shan)
all_p %>%
    ggplot(aes(x=exp_shan, y=PABD)) +
    geom_point(aes(col=species), cex=3)  

#' Now for observed OTUs
#' \

glm_PABD_rich <- glm(PABD ~ species*exp_rich, data=all_p, family=binomial(link="logit"))
Anova(glm_PABD_rich)
all_p %>%
    ggplot(aes(x=exp_rich, y=PABD)) +
    geom_point(aes(col=species), cex=3)

#' Clearly, we see that neither diversity nor richness seems to signficiatly affect
#' infection risk. 
#'\
#'However, let's try normalizing it by species; using each species to determine whether it is "out of the
#'ordinary"

glm_PABD_pshan <- glm(PABD ~ species*p_shan, data=all_p, family=binomial(link="logit"))
Anova(glm_PABD_pshan)
all_p %>%
    ggplot(aes(x=p_shan, y=PABD)) +
    geom_point(aes(col=species), cex=3) 

#' Now for observed OTUs
#' \

glm_PABD_prich <- glm(PABD ~ species*p_rich, data=all_p, family=binomial(link="logit"))
Anova(glm_PABD_prich)
all_p %>%
    ggplot(aes(x=p_rich, y=PABD)) +
    geom_point(aes(col=species), cex=3)  

#### TESTING SOMETHING
Anova(glm(PABD ~ p_shan, data=all_p, family=binomial(link="logit")))
Anova(glm(PABD ~ p_rich, data=all_p, family=binomial(link="logit")))

#'\
#' If anything, it looks like increased diversity and richness might increase infection risk


#### eBD and diversity ####

#'\
#' (1b) Does overall diversity of microbiome influence BD infection intensity?\
#' \
#' The next thing we would like to know is if diversity of richness of the microbiome
#' influences infection intensity.
#' \

lm_eBD_shan <- lm(eBD_log ~ species*exp_shan, data=all_p)
Anova(lm_eBD_shan)
all_p %>%
    ggplot(aes(x=exp_shan, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 

#' Now let's do richness
lm_eBD_rich <- lm(eBD_log ~ species*exp_rich, data=all_p)
Anova(lm_eBD_rich)
all_p %>%
    ggplot(aes(x=exp_rich, y=eBD_log)) +
    geom_point(aes(col=species), cex=3)  

#' \ 
#' It looks like diversity might slightly increasing risk of infection, which is opposite what we might expect
#' Let's try using the normalized values.

lm_eBD_pshan <- lm(eBD_log ~ species*p_shan, data=all_p)
Anova(lm_eBD_pshan)
all_p %>%
    ggplot(aes(x=p_shan, y=eBD_log)) +
    geom_point(aes(col=species), cex=3)

#' Now let's do richness
lm_eBD_prich <- lm(eBD_log ~ species*p_rich, data=all_p)
Anova(lm_eBD_prich)
all_p %>%
    ggplot(aes(x=p_shan, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 

#### TESTING SOMETHING
Anova(lm(eBD_log ~ p_shan, data=all_p))
Anova(lm(eBD_log ~ p_rich, data=all_p))

#### PABD and instability ####

#' (2a) Does instability of microbiome influence BD infection rate?\
#' Here we look at average distance travelled (bray-curtis) between samples
#' prior to being infected. We see if it is correlated to infection risk.

glm_PABD_bc <- glm(PABD ~ species*exp_mu, data=all_p, family=binomial)
Anova(glm_PABD_bc)
all_p %>%
    ggplot(aes(x=exp_mu, y=PABD)) +
    geom_point(aes(col=species), cex=3) 

glm_PABD_pbc <- glm(PABD ~ species*p_mu, data=all_p, family=binomial)
Anova(glm_PABD_pbc)
all_p %>%
    ggplot(aes(x=p_mu, y=PABD)) +
    geom_point(aes(col=species), cex=3) 

### TESTING SOMETHING
Anova( glm(PABD ~ p_mu, data=all_p, family=binomial))
#### eBD and instability ####

#' 
#' \
#' (2b) Does instability of microbiome influence BD infection intensity?\

lm_eBD_bc <- lm(eBD_log ~ species*exp_mu, data=all_p)
Anova(lm_eBD_bc)
all_p %>%
    ggplot(aes(x=exp_mu, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 

lm_BD_pbc <- lm(eBD_log ~ species*p_mu, data=all_p)
Anova(lm_BD_pbc)
all_p %>%
    ggplot(aes(x=p_mu, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 

#### TESTING SOMETHING
Anova(lm(eBD_log ~ exp_mu, data=all_p))
Anova(lm(eBD_log ~ p_mu, data=all_p))


#### PABD and inhibitory ####

#' 
#' \
#' (3a) Does composition of microbiome influence BD infection risk?\
#' Now, we ask if composition-- specitically, the richness and percent of BD inhibitory bacteria-- 
#' influences infection risk in individuals. First, below, we use just a regular correlation
#' between richness and infection risk

glm_PABD_inhibRich <- glm(PABD ~ species*exp_inhibRich, data=all_p, family=binomial)
Anova(glm_PABD_inhibRich)
all_p %>%
    ggplot(aes(x=exp_inhibRich, y=PABD)) +
    geom_point(aes(col=species), cex=3)

#' Now let's do percent inhibitory
glm_PABD_percInhib <- glm(PABD ~ species*exp_pinhib, data=all_p, family=binomial)
Anova(glm_PABD_percInhib)
summary(glm_PABD_percInhib)
all_p %>%
    ggplot(aes(x=exp_pinhib, y=PABD)) +
    geom_point(aes(col=species), cex=3) 

#' Now let's look at the standardized values

glm_PABD_pinhibRich <- glm(PABD ~ species*p_inhibRich, data=all_p, family=binomial)
Anova(glm_PABD_pinhibRich) #### SIG
all_p %>%
    ggplot(aes(x=p_inhibRich, y=PABD)) +
    geom_point(aes(col=species), cex=3)

#' Now let's do percent inhibitory of standardized values
glm_PABD_ppinhib <- glm(PABD ~ species*p_pinhib, data=all_p, family=binomial)
Anova(glm_PABD_ppinhib)
summary(glm_PABD_ppinhib)
all_p %>%
    ggplot(aes(x=p_pinhib, y=PABD)) +
    geom_point(aes(col=species), cex=3)

#### TESTING SOMETHING
Anova(glm(PABD ~ p_inhibRich, data=all_p, family=binomial))
Anova(glm(PABD ~ p_pinhib, data=all_p, family=binomial))

#### eBD and inhibitory ####

#' (3b) Does composition of microbiome influence BD infection intensity?\

lm_eBD_inhibRich <- lm(eBD_log ~ species*exp_inhibRich, data=all_p)
Anova(lm_eBD_inhibRich)
all_p %>%
    ggplot(aes(x=exp_inhibRich, y=eBD_log)) +
    geom_point(aes(col=species), cex=3)

#' Now let's do percent inhibitory
lm_eBD_percInhib <- lm(eBD_log ~  species*exp_pinhib, data=all_p)
Anova(lm_eBD_percInhib)
all_p %>%
    ggplot(aes(x=exp_pinhib, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 

#' Now let's look at the standardized values

lm_eBD_pinhibRich <- lm(eBD_log ~ species*p_inhibRich, data=all_p)
Anova(lm_eBD_pinhibRich)
all_p %>%
    ggplot(aes(x=p_inhibRich, y=eBD_log)) +
    geom_point(aes(col=species), cex=3)

#' Now let's do percent inhibitory of standardized values
lm_eBD_ppinhib <- lm(eBD_log ~  species*p_pinhib, data=all_p)
Anova(lm_eBD_ppinhib)
all_p %>%
    ggplot(aes(x=p_pinhib, y=eBD_log)) +
    geom_point(aes(col=species), cex=3) 

### Test something 
Anova(lm(eBD_log ~ p_inhibRich, data=all_p))
Anova(lm(eBD_log ~  p_pinhib, data=all_p))

#' 
#' Part II: Affect of BD infection on microbiome state
#' 
#' (1a) Does BD infection state affect microbiome diversity?
#' - OTU richness vs BD infection
#' - Chao1 richness vs BD infection
#' - Shannon richness vs BD infection
#' - PD vs BD infection

lm_shan_PABD <- lm(shannon ~ species*PABD, data=all_p_infected)
Anova(lm_shan_PABD)
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=shannon)) +
    geom_violin() +
    geom_point(aes(col=species), cex=3, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)



lm_rich_PABD <- lm(logRich ~ species*PABD, data=all_p_infected)
Anova(lm_rich_PABD)
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=logRich)) +
    geom_violin() +
    geom_point(aes(col=species), cex=3, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)



#' Now let's try the standardized values
lm_pshan_PABD <- lm(p_shan ~ species*PABD, data=all_p_infected)
Anova(lm_pshan_PABD)
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=p_shan)) +
    geom_violin() +
    geom_point(aes(col=species), cex=3, position=position_jitter(width=0.15, height=0)) +
    facet_wrap(~species, nrow=1)

lm_prich_PABD <- lm(p_rich ~ species*PABD, data=all_p_infected)
Anova(lm_prich_PABD)
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=p_rich)) +
    geom_violin() +
    geom_point(aes(col=species), cex=3, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)


#' 
#' (1b) Does BD infection intensity affect microbiome diversity?

lm_shan_eBD <- lm(shannon ~ species*eBD_log, data=all_p_infected)
Anova(lm_shan_eBD)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=shannon)) +
    geom_point(aes(col=species), cex=3)

lm_rich_eBD <- lm(logRich ~ species*eBD_log, data=all_p_infected)
Anova(lm_rich_eBD)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=logRich)) +
    geom_point(aes(col=species), cex=3)

#' Now let's try the standardized values
lm_pshan_eBD <- lm(p_shan ~ species*eBD_log, data=all_p_infected)
Anova(lm_pshan_eBD)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_shan)) +
    geom_point(aes(col=species), cex=3) 

lm_prich_eBD <- lm(p_rich ~ species*eBD_log, data=all_p_infected)
Anova(lm_prich_eBD)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_rich)) +
    geom_point(aes(col=species), cex=3)

#' (2a) Does BD infection state affect microbiome instability?

lm_bc_PABD <- lm(distance_bray_curtis ~ species*PABD, data=all_p_infected)
Anova(lm_bc_PABD)
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=distance_bray_curtis)) +
    geom_violin() +
    geom_point(aes(color=species), cex=4, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)

# try standardized
lm_pbc_PABD <- lm(p_BC ~ species*PABD, data=all_p_infected)
Anova(lm_pbc_PABD)
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=p_BC)) +
    geom_violin() +
    geom_point(aes(color=species), cex=4, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)

#' (2b) Does BD infection intensity affect microbiome instability?
#' 
lm_bc_eBD <- lm(distance_bray_curtis ~ species*eBD_log, data=all_p_infected)
Anova(lm_bc_eBD)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=distance_bray_curtis)) +
    geom_point(aes(color=species), cex=4)

# try standardized
lm_pbc_eBD <- lm(p_BC ~ species*eBD_log, data=all_p_infected)
Anova(lm_pbc_eBD)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_BC)) +
    geom_point(aes(color=species), cex=4)

#' (3a) Does BD infection state affect microbiome composition?

lm_inhibRich_PABD <- lm(inhibRich ~ species*PABD, data=all_p_infected)
Anova(lm_inhibRich_PABD)
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=inhibRich)) +
    geom_violin() +
    geom_point(aes(col=species), cex=3, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)


lm_percInhib_PABD <- lm(percInhib ~ species*PABD, data=all_p_infected)
Anova(lm_percInhib_PABD)
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=percInhib)) +
    geom_violin() +
    geom_point(aes(col=species), cex=3, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)


#' Now let's try the standardized values
lm_pinhibRich_PABD <- lm(p_inhibRich ~ species*PABD, data=all_p_infected)
Anova(lm_pinhibRich_PABD)
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=p_inhibRich)) +
    geom_violin() +
    geom_point(aes(col=species), cex=3, position=position_jitter(width=0.15, height=0)) +
    facet_wrap(~species, nrow=1)

lm_ppercInhib_PABD <- lm(p_percInhib ~ species*PABD, data=all_p_infected)
Anova(lm_ppercInhib_PABD)
all_p_infected %>%
    mutate(PABD = factor(PABD)) %>%
    ggplot(aes(x=PABD, y=p_percInhib)) +
    geom_violin() +
    geom_point(aes(col=species), cex=3, position=position_jitter(width=0.15, height=0))+
    facet_wrap(~species, nrow=1)


#' 
#' (3b) Does BD infection intensity affect microbiome composition?

lm_inhibRich_eBD <- lm(inhibRich ~ species*eBD_log, data=all_p_infected)
Anova(lm_inhibRich_eBD)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=inhibRich)) +
    geom_point(aes(col=species), cex=3)

lm_percInhib_eBD <- lm(percInhib ~ species*eBD_log, data=all_p_infected)
Anova(lm_percInhib_eBD)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=percInhib)) +
    geom_point(aes(col=species), cex=3)



#' Now let's try the standardized values
lm_pinhibRich_eBD <- lm(p_inhibRich ~ species*eBD_log, data=all_p_infected)
Anova(lm_pinhibRich_eBD)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_inhibRich)) +
    geom_point(aes(col=species), cex=3) 

lm_ppercInhib_eBD <- lm(p_percInhib ~  species*eBD_log, data=all_p_infected)
Anova(lm_ppercInhib_eBD)
all_p_infected %>%
    ggplot(aes(x=eBD_log, y=p_percInhib)) +
    geom_point(aes(col=species), cex=3)
