### Creating figures for 5 species manuscript ####

library(tidyverse)
library(ggplot2)
library(gridExtra)

# Make directory for figures
dir.create("FIGURES")

load("mf_con_without_init_infect.RData")
load("mf_treat_without_init_infect.RData")
load("mf.rare.RData")
# OTU table of inhibitory bacteria
load("otu.inhibOnly.treat.RData")
load("otu.inhibOnly.con.RData")

# Previous analyses summaries
load("all_p.RData")
load("all_p_infected.RData")
load("con_exp_indiv.RData")

# add a species column and PABD column
all_p <- all_p %>%
    mutate(PABD=ifelse(infect>0,1,0), infect = log(infect+1)) %>%
    rename(eBD_log=infect) %>%
    separate(toadID, into=c("species","indiv"), remove=FALSE)
all_p_infected <- all_p_infected %>%
    mutate(PABD=ifelse(eBD_log>0,1,0)) %>%
    separate(toadID, into=c("species","indiv"), remove=FALSE)

#### experimental design ####
pdf(file = "FIGURES/experimental_design.pdf", width = 5,height = 12)
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
dev.off()



#### Control data: how does it vary? ####
gg_NMDS <- mf_con_without_init_infect %>%
    ggplot(aes(x=NMDS1, y=NMDS2)) +
    geom_point(aes(col=species), cex=3, show.legend = FALSE)
gg_infect <- mf_treat_without_init_infect  %>%
    ggplot(aes(x=species, y=eBD_log)) +
    geom_point(aes(col=species), cex=2.5, position = position_jitter(width=0.1, height=0.05), show.legend = FALSE)+
    xlab("Species") +
    ylab("log Bd Load")

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
    dplyr::select(species, distance_bray_curtis) %>%
    mutate(metric="Bray_Curtis_Distance")%>%
    rename(value=distance_bray_curtis)

gg_all <- rbind(temp1,temp2,temp3,temp4) %>%
    ggplot(aes(x=species, y=value)) +
    geom_boxplot() +
    geom_point(aes(col=species), position = position_jitter(width=0.1, height=0))+
    facet_grid(metric~., scales = "free") +
    ylab("")+
    xlab("Species")
lay <- rbind(c(1,2),
             c(3,2))

pdf("FIGURES/data_summary_controls.pdf", height = 10, width = 10)
grid.arrange(gg_NMDS, gg_all, gg_infect, layout_matrix = lay)
dev.off()


#### Inihibitory Richness affects chances of infection ####

# Note; should incorporate species
Anova(glm(PABD ~ species + p_inhibRich, family = binomial(link="logit"), data=all_p))
glm_PABD_inhibRich <- glm(PABD ~ p_inhibRich, family = binomial(link="logit"), data=all_p)

# x.predict <- data.frame(p_inhibRich=rep(seq(0,1, length.out = 100),times=5), species=rep(c("Anbo","Anma","Lica","Lipi","Osse"), each=100))
x.predict <- data.frame(p_inhibRich=seq(0,1, length.out = 100))
y.predict <- predict(glm_PABD_inhibRich, newdata = x.predict, type = "link", se.fit = TRUE)

# Combine the hypothetical data and predicted values
new.data <- cbind(x.predict, y.predict)

# Calculate confidence intervals
std <- qnorm(0.95 / 2 + 0.5)
new.data$ymin <- glm_PABD_inhibRich$family$linkinv(new.data$fit - std * new.data$se)
new.data$ymax <- glm_PABD_inhibRich$family$linkinv(new.data$fit + std * new.data$se)
new.data$fit <- glm_PABD_inhibRich$family$linkinv(new.data$fit)  # Rescale to 0-1

pdf("FIGURES/inhibRich_PABD.pdf",height=5, width=7)
all_p %>%
    rename(Species=species) %>%
    ggplot() +
    geom_point(aes(x=p_inhibRich, y=PABD,col=Species), position=position_jitter(height=0.05, width=0)) +
    geom_ribbon(data=new.data, aes(x=p_inhibRich, ymin=ymin, ymax=ymax), alpha=0.2, col="lightgrey") +
    geom_line(data=new.data, aes(x=p_inhibRich, y=fit)) +
    ylab("Probability of infection") + 
    xlab("Inhibitory OTU Richness (Percentile of species)")
dev.off()

#### Percent inhibitory changes after infection ####






#### Correlation between Richness and Inhibitory Richness ####




all_p %>%
    dplyr::select(toadID, exp_rich, p_rich, exp_inhibRich, p_inhibRich) %>%
    # full_join(con_exp_indiv, by = "toadID") %>%
    rbind(con_exp_indiv) %>%
    dplyr::select(p_inhibRich, p_rich, toadID) %>%
    separate(toadID, into=c("species","n"), remove=FALSE) %>%
    ggplot(aes(x=p_rich, y=p_inhibRich)) +
    geom_point(aes(col=species), cex=3)