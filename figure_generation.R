### Creating figures for 5 species manuscript ####

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(car) #Anova
library(RColorBrewer) # colors for barplots
library(MASS) #fitdistr

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
    mutate(Contaminated = factor(ifelse(orig_contam ==1, "!Contam upon arrival",NA), levels=c("!Contam upon arrival"))) %>%
    ggplot(aes(x=time, y=indiv)) +
    geom_line(aes(group=toadID, col=Treatment)) +
    geom_point(aes(group=toadID,bg=LnBd_load), cex=4, pch=21)+
    scale_color_manual(values=c("black","blue","orange")) +
    scale_fill_gradient(low = "white", high = "red") +
    geom_vline(aes(xintercept=5.5), col="orange")+
    geom_point(aes(group=toadID, col=Contaminated), cex=1, pch=19)+ ## NEW LINE
    facet_wrap(~species, nrow=5) +
    xlab("Time") +
    ylab("Individual Toad")
dev.off()



#### Control data: how does it vary? ####
# Calculate the hulls for each group
# hull_toad <- mf_con_without_init_infect %>%
#     group_by(toadID) %>%
#     slice(chull(NMDS1, NMDS2))

gg_NMDS <- mf_con_without_init_infect %>%
    ggplot(aes(x=NMDS1, y=NMDS2)) +
    # geom_line(aes(group=toadID,x=NMDS1, y=NMDS2), col="lightgrey") +
    geom_point(aes(col=species, alpha=(time)),  cex=2.5, show.legend = FALSE) +
    # geom_polygon(data = hull_toad, aes(fill=species), alpha = 0.1) +
    theme_classic()
# gg_infect <- mf_treat_without_init_infect  %>%
#     ggplot(aes(x=species, y=eBD_log)) +
#     geom_point(aes(col=species), cex=2, position = position_jitter(width=0.1, height=0.05), show.legend = FALSE)+
#     xlab("Species") +
#     ylab("log Bd Load") +
#     theme_classic()
gg_infect <- mf_treat_without_init_infect  %>%
    group_by(species, toadID) %>%
    summarize(eBD_log = max(eBD_log)) %>%
    ggplot(aes(x=species, y=eBD_log)) +
    geom_point(aes(col=species), cex=2, position = position_jitter(width=0.3, height=0), show.legend = FALSE)+
    xlab("Species") +
    ylab("Max log Bd Load of each individual") +
    theme_classic()

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
    rename(Species=species) %>%
    mutate(metric = factor(metric, levels=c("log_OTU_Richness","Inhibitory_OTU_Richness","Percent_Inhibitory","Bray_Curtis_Distance"))) %>%
    mutate(Metric = gsub("_"," ",metric, fixed=TRUE)) %>%
    ggplot(aes(x=Species, y=value)) +
    geom_boxplot() +
    geom_point(aes(col=Species), position = position_jitter(width=0.1, height=0), alpha=1/3)+
    facet_grid(Metric~., scales = "free", switch="y") +
    ylab("")+
    xlab("Species") 
lay <- rbind(c(1,2),
             c(3,2))

pdf("FIGURES/data_summary_controls.pdf", height = 8, width = 8)
grid.arrange(gg_NMDS, gg_all, gg_infect, layout_matrix = lay)
dev.off()

# Note about data summary control:
# The NMDS is all time points for control (uninfected) only; fading is time.
# The infection is MAX infection
# the others are all time points for control

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

pdf("FIGURES/inhibRich_PABD.pdf",height=4, width=6)
all_p %>%
    rename(Species=species) %>%
    ggplot() +
    geom_point(aes(x=p_inhibRich, y=PABD,col=Species), position=position_jitter(height=0.05, width=0)) +
    geom_ribbon(data=new.data, aes(x=p_inhibRich, ymin=ymin, ymax=ymax), alpha=0.2, col="lightgrey") +
    geom_line(data=new.data, aes(x=p_inhibRich, y=fit)) +
    ylab("Probability of infection") + 
    xlab("Inhibitory OTU Richness (Percentile of species)") +
    theme_classic()
dev.off()

#### Percent inhibitory changes after infection ####

#Exact same?
otu.inhibOnly.con$OTUID == otu.inhibOnly.treat$OTUID
# Good

# Look at "must abundant" otus
inhibPWD <- "../ANTIFUNGAL/inhibitory_metadata_MANUAL.txt"
inhib <- read.delim(paste0(inhibPWD), header=FALSE, as.is=TRUE)

inhib.tb <- inhib %>%
    as_tibble() %>%
    rename(Name=V1, inhib=V2, num=V3, seq=V4,sampleID=V5, taxa=V6)
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


## Control
# Going to leave zeros in for now, bc only 3-5 in each control and treatment
temp <- otu.inhibOnly.con %>%
    dplyr::select(-OTUID) %>%
    unlist()
otu_long_con <- cbind(rep(colnames(otu.inhibOnly.con)[-ncol(otu.inhibOnly.con)], each=nrow(otu.inhibOnly.con))
                        ,temp
                        ,rep(otu.inhibOnly.con$OTUID,times=ncol(otu.inhibOnly.con)-1)) %>%
    as_tibble() %>%
    rename(SampleID=V1, reads=temp, OTUID=V2) %>%
    mutate(reads=as.numeric(reads))
mf_con_with_inhibOTUs <- mf_con_without_init_infect %>%
    dplyr::select(SampleID, species, time, toadID, eBD_log, PABD, prepost, BD_infected) %>%
    left_join(otu_long_con)%>%
    mutate(taxa = (inhib.tb[match(OTUID, inhib.tb$seq),"taxa"])$taxa ) %>%
    separate(taxa, into = c("K","P","C","O","F","G"), sep = ";", remove = FALSE, fill="left")


## PLot
mf_con_with_inhibOTUs_summarized <- mf_con_with_inhibOTUs %>%
    group_by(toadID, time, species, G) %>%
    summarise(AverageProportion = mean(reads)) %>%
    mutate(Group="Control")
mf_treat_with_inhibOTUs_summarized <- mf_treat_with_inhibOTUs %>%
    group_by(toadID, time, species, G) %>%
    summarise(AverageProportion = mean(reads)) %>%
    mutate(Group="Treatment")
mf_combined_with_inhibOTUS_summarized <- full_join(mf_con_with_inhibOTUs_summarized, mf_treat_with_inhibOTUs_summarized, by=c("time","species","G","AverageProportion","Group"))


# 
# gg_conpercinhib <- mf_con_with_inhibOTUs %>%
#     group_by(toadID, time, species, G) %>%
#     summarise(AverageProportion = mean(reads)) %>%
#     ggplot(aes(x=time, y=AverageProportion)) +
#     geom_bar(aes(fill=G), stat="identity", show.legend = FALSE) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.key.size =  unit(0.2, "cm") ) +
#     facet_wrap(~species, nrow=1) +
#     scale_fill_manual(values=set_col) +
#     geom_vline(aes(xintercept=5.5), col="grey", lty=2)
# gg_treatpercinhib <- mf_treat_with_inhibOTUs %>%
#     group_by(toadID, time, species, G) %>%
#     summarise(AverageProportion = mean(reads)) %>%
#     ggplot(aes(x=time, y=AverageProportion)) +
#     geom_bar(aes(fill=G), stat="identity", show.legend = FALSE) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.key.size =  unit(0.2, "cm") ) +
#     facet_wrap(~species, nrow=1) +
#     scale_fill_manual(values=set_col) +
#     geom_vline(aes(xintercept=5.5), col="grey", lty=2)
# 
g_legend <- function(a.gplot){ # from stack overflow
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
gg_legend <- g_legend(mf_treat_with_inhibOTUs %>%
                          group_by(toadID, time, species, G) %>%
                          summarise(AverageProportion = mean(reads)) %>%
                          rename(Genus=G) %>%
                          ggplot(aes(x=time, y=AverageProportion)) +
                          geom_bar(aes(fill=Genus), stat="identity") +
                          scale_fill_manual(values=set_col) )

gg_combined_inhibOTUs <- mf_combined_with_inhibOTUS_summarized %>%
    rename(Time=time, Species=species, Genus=G) %>%
    ggplot(aes(x=Time, y=AverageProportion)) +
    geom_bar(aes(fill=Genus), stat="identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.key.size =  unit(0.2, "cm") ) +
    facet_grid(Group~Species) +
    scale_fill_manual(values=set_col) +
    geom_vline(aes(xintercept=5.5), col="grey", lty=2) +
    theme_classic()
lay <- rbind(c(1,1,2))

pdf(file = "FIGURES/percent_inhib.pdf", width=10, height=5.5)
grid.arrange(gg_combined_inhibOTUs, gg_legend,layout_matrix = lay )
dev.off()

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

# mf_con_statdiff %>%
#     filter(!is.na(significant)) %>%
#     ggplot(aes(x=G)) +
#     geom_point(aes(y=temp, col=species, pch=significant), cex=2, position=position_jitter(width=0, height=0.25)) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     scale_shape_manual(values=c(21,19)) +
#     geom_hline(aes(yintercept=0))


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

# mf_treat_statdiff %>%
#     filter(!is.na(significant)) %>%
#     ggplot(aes(x=G)) +
#     geom_point(aes(y=temp, col=species, pch=significant), cex=2, position=position_jitter(width=0, height=0.25)) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     scale_shape_manual(values=c(21,19)) +
#     geom_hline(aes(yintercept=0))

####### Looking at specific OTUs #########

# CONTROL 
# mf_con_statdiff <- mf_con_with_inhibOTUs %>%
#     group_by(toadID, species, prepost, G) %>%
#     summarize(meanReads=mean(reads), meanBd = mean(eBD_log)) %>%
#     dplyr::select(-meanBd) %>%
#     spread(key=prepost, value=meanReads) %>%
#     mutate(fc=log((Pos-Pre)/Pre +1), significant = NA, p=NA) %>%
#     mutate(temp= ifelse(is.finite(fc), fc, ifelse(is.infinite(fc), 10, 0))) 
mf_con_statdiff <- mf_con_with_inhibOTUs %>%
    group_by(toadID, species, prepost, G) %>%
    summarize(meanReads=mean(reads), meanBd = mean(eBD_log)) %>%
    dplyr::select(-meanBd) %>%
    spread(key=prepost, value=meanReads) %>%
    mutate(fc=log((Pos-Pre)/Pre +1), significant = NA, p=NA) %>%
    mutate(temp= ifelse(is.finite(fc), fc, ifelse(fc>0, 10, -10)))

    
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


# TREAT

mf_treat_statdiff <- mf_treat_with_inhibOTUs %>%
    group_by(toadID, species, prepost, G) %>%
    summarize(meanReads=mean(reads), meanBd = mean(eBD_log)) %>%
    dplyr::select(-meanBd) %>%
    spread(key=prepost, value=meanReads) %>%
    # mutate(diff=Pos-Pre, significant = NA, p=NA) %>%
    mutate(fc=log((Pos-Pre)/Pre +1), significant = NA, p=NA) %>%
    mutate(temp= ifelse(is.finite(fc), fc, ifelse(fc>0, 10, -10)))

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

## Plot

mf_con_statdiff$Group <- "Control"
mf_treat_statdiff$Group <- "Treatment"
mf_combined_statdiff <- full_join(mf_con_statdiff, mf_treat_statdiff)

pdf("FIGURES/otu_change.pdf", width=5, height=5)
mf_combined_statdiff %>%
    rename(FoldChange=temp, Species=species, Significant=significant, Genus=G) %>%
    filter(!is.na(Significant)) %>%
    ggplot(aes(x=Genus)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)
          , panel.background = element_rect(fill = "white",
                                        colour = "black")
          , panel.grid.major.x = element_line(colour = "grey")
          ) +
    geom_point(aes(y=FoldChange, col=Species, pch=Significant), cex=2, position=position_jitter(width=0.05, height=0.5)) +
    scale_shape_manual(values=c(21,19)) +
    geom_hline(aes(yintercept=c(0)), lty=1, col="black") +
    geom_hline(aes(yintercept=c(-10)), lty=3, col="grey") +
    geom_hline(aes(yintercept=c(10)), lty=3, col="grey") +
    facet_grid(Group ~ .) 
dev.off()




#### Correlation between Richness and Inhibitory Richness ####

pdf("FIGURES/corr_rich_inhibrich.pdf", width=5, height=4)
all_p %>%
    dplyr::select(toadID, exp_rich, p_rich, exp_inhibRich, p_inhibRich) %>%
    # full_join(con_exp_indiv, by = "toadID") %>%
    rbind(con_exp_indiv) %>%
    dplyr::select(p_inhibRich, p_rich, toadID) %>%
    separate(toadID, into=c("species","n"), remove=FALSE) %>%
    rename(Species=species) %>%
    ggplot(aes(x=p_rich, y=p_inhibRich)) +
    geom_point(aes(col=Species), cex=3) +    
    ylab("OTU Richness (Percentile of species)") + 
    xlab("Inhibitory OTU Richness (Percentile of species)")
dev.off()


##### SUPPLEMENTARY FIGURES ######

#### Fit models for each trait ####

## Logrich

logrich.normfit <- fitdistr(x=mf_con_without_init_infect$logRich, densfun = "normal")
x.pred <- seq(min(mf_con_without_init_infect$logRich)-sd(mf_con_without_init_infect$logRich)
              , max(mf_con_without_init_infect$logRich)+sd(mf_con_without_init_infect$logRich)
              , length.out = 100)
y.expect <- dnorm(x=x.pred, mean=logrich.normfit$estimate[1], sd=logrich.normfit$estimate[2])

gg_logrich_fit <- mf_con_without_init_infect %>%
    ggplot(aes(x=logRich, y=..density..)) +
    geom_histogram(bins=25) +
    geom_line(data=data.frame(x=x.pred, y=y.expect), aes(x=x, y=y), col="red")


## Dispersion

disp.normfit <- fitdistr(x=mf_con_without_init_infect$distance.to.centroid, densfun = "normal")
disp.lognormfit <- fitdistr(x=mf_con_without_init_infect$distance.to.centroid, densfun = "lognormal")
x.pred <- seq(min(mf_con_without_init_infect$distance.to.centroid, na.rm = TRUE)-sd(mf_con_without_init_infect$distance.to.centroid, na.rm = TRUE)
              , max(mf_con_without_init_infect$distance.to.centroid, na.rm = TRUE)+sd(mf_con_without_init_infect$distance.to.centroid, na.rm = TRUE)
              , length.out = 100)
y.expect.norm <- dnorm(x=x.pred, mean=disp.normfit$estimate[1], sd=disp.normfit$estimate[2])
y.expect.lognorm <- dlnorm(x=x.pred, meanlog=disp.lognormfit$estimate[1], sdlog = disp.lognormfit$estimate[2])

gg_disp_fit <- mf_con_without_init_infect %>%
    ggplot(aes(x=distance.to.centroid, y=..density..)) +
    geom_histogram(bins=25) +
    geom_line(data=data.frame(x=x.pred, y=y.expect.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.pred, y=y.expect.lognorm), aes(x=x, y=y), col="blue") 


## Instability

instab.normfit <- fitdistr(x=mf_con_without_init_infect$distance_bray_curtis[!is.na(mf_con_without_init_infect$distance_bray_curtis)], densfun = "normal")
instab.lognormfit <- fitdistr(x=mf_con_without_init_infect$distance_bray_curtis[!is.na(mf_con_without_init_infect$distance_bray_curtis)], densfun = "beta", start = list(shape1=1,shape2=1))
x.pred <- seq(min(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)-sd(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)
              , max(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)+sd(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)
              , length.out = 100)
y.expect.norm <- dnorm(x=x.pred, mean=instab.normfit$estimate[1], sd=instab.normfit$estimate[2])
y.expect.beta <- dbeta(x=x.pred, shape1 = instab.lognormfit$estimate[1], shape2 = instab.lognormfit$estimate[2])

mf_con_without_init_infect %>%
    ggplot(aes(x=distance_bray_curtis, y=..density..)) +
    geom_histogram(bins=25) +
    geom_line(data=data.frame(x=x.pred, y=y.expect.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.pred, y=y.expect.beta), aes(x=x, y=y), col="blue") 


## Instability

instab.normfit <- fitdistr(x=mf_con_without_init_infect$distance_bray_curtis[!is.na(mf_con_without_init_infect$distance_bray_curtis)], densfun = "normal")
instab.lognormfit <- fitdistr(x=mf_con_without_init_infect$distance_bray_curtis[!is.na(mf_con_without_init_infect$distance_bray_curtis)], densfun = "beta", start = list(shape1=1,shape2=1))
x.pred <- seq(min(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)-sd(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)
              , max(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)+sd(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)
              , length.out = 100)
y.expect.norm <- dnorm(x=x.pred, mean=instab.normfit$estimate[1], sd=instab.normfit$estimate[2])
y.expect.beta <- dbeta(x=x.pred, shape1 = instab.lognormfit$estimate[1], shape2 = instab.lognormfit$estimate[2])

gg_instab_fit <- mf_con_without_init_infect %>%
    ggplot(aes(x=distance_bray_curtis, y=..density..)) +
    geom_histogram(bins=25) +
    geom_line(data=data.frame(x=x.pred, y=y.expect.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.pred, y=y.expect.beta), aes(x=x, y=y), col="blue") 

## InhibRich

?fitdistr
inhibRich.normfit <- fitdistr(x=mf_con_without_init_infect$inhibRich[!is.na(mf_con_without_init_infect$inhibRich)], densfun = "normal")
inhibRich.posfit <- fitdistr(x=mf_con_without_init_infect$inhibRich[!is.na(mf_con_without_init_infect$inhibRich)], densfun = "Poisson")
x.pred <- round(seq(min(mf_con_without_init_infect$inhibRich, na.rm = TRUE)-sd(mf_con_without_init_infect$inhibRich, na.rm = TRUE)
              , max(mf_con_without_init_infect$inhibRich, na.rm = TRUE)+sd(mf_con_without_init_infect$inhibRich, na.rm = TRUE)
              , length.out = 100))
y.expect.norm <- dnorm(x=x.pred, mean=inhibRich.normfit$estimate[1], sd=inhibRich.normfit$estimate[2])
y.expect.pois <- dpois(x=x.pred, lambda = inhibRich.posfit$estimate[1])

gg_inhibRich_fit <- mf_con_without_init_infect %>%
    ggplot(aes(x=inhibRich, y=..density..)) +
    geom_histogram(bins=25) +
    geom_line(data=data.frame(x=x.pred, y=y.expect.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.pred, y=y.expect.pois), aes(x=x, y=y), col="purple") 


