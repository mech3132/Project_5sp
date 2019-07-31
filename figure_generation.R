### Creating figures for 5 species manuscript ####

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(car) #Anova
library(RColorBrewer) # colors for barplots
library(MASS) #fitdistr
library(vegan) # braycurtis distance for nmds of inhib otus

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
load("all_p_withcon.RData")

# add a species column and PABD column
all_p <- all_p %>%
    mutate(PABD=ifelse(infect>0,1,0), infect = log(infect+1)) %>%
    rename(eBD_log=infect) %>%
    separate(toadID, into=c("species","indiv"), remove=FALSE)
all_p_infected <- all_p_infected %>%
    mutate(PABD=ifelse(eBD_log>0,1,0)) %>%
    separate(toadID, into=c("species","indiv"), remove=FALSE)

#### experimental design ####
pdf(file = "FIGURES/experimental_design.pdf", width = 6,height = 12)
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
    ylab("Individual Toad") + 
    theme_classic()
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
    mutate(metric="log_ASV_Richness") %>%
    rename(value=logRich)
temp2 <- mf_con_without_init_infect %>%
    dplyr::select(species, inhibRich) %>%
    mutate(metric="Inhib_ASV_Richness")%>%
    rename(value=inhibRich)
temp3 <- mf_con_without_init_infect %>%
    dplyr::select(species, percInhib) %>%
    mutate(metric="Percent_Inhibitory")%>%
    rename(value=percInhib)
temp4 <- mf_con_without_init_infect %>%
    dplyr::select(species, disper_bray_curtis) %>%
    mutate(metric="Dispersion")%>%
    rename(value=disper_bray_curtis)
temp5 <- mf_con_without_init_infect %>%
    dplyr::select(species, distance_bray_curtis) %>%
    mutate(metric="Instability")%>%
    rename(value=distance_bray_curtis)

gg_all <- rbind(temp1,temp2,temp3,temp4, temp5) %>%
    rename(Species=species) %>%
    mutate(Metric = gsub("_"," ",metric, fixed=TRUE)) %>%
    mutate(Metric = factor(Metric, levels=c("log ASV Richness","Inhib ASV Richness","Percent Inhibitory","Dispersion", "Instability"))) %>%
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
    xlab("Inhibitory ASV Richness (Percentile of species)") +
    theme_classic()
dev.off()


#### Dispersion affects infection intensity ####

# Note; should incorporate species
Anova(lm(eBD_log ~ species + p_disper, data=all_p))

# # Calculate confidence intervals
# std <- qnorm(0.95 / 2 + 0.5)
# new.data$ymin <- glm_eBDlog_disper$
# new.data$ymin <- glm_eBDlog_disper$family$linkinv(new.data$fit - std * new.data$se)
# new.data$ymax <- glm_eBDlog_disper$family$linkinv(new.data$fit + std * new.data$se)
# new.data$fit <- glm_eBDlog_disper$family$linkinv(new.data$fit)  # Rescale to 0-1

pdf("FIGURES/disper_eBDlog.pdf",height=4, width=6)
all_p %>%
    rename(Species=species) %>%
    ggplot(aes(x=p_disper, y=eBD_log))+
    geom_point(aes(col=Species), position=position_jitter(height=0, width=0)) +
    # geom_ribbon(data=new.data, aes(x=p_inhibRich, ymin=ymin, ymax=ymax), alpha=0.2, col="lightgrey") +
    # geom_line(data=new.data, aes(x=p_inhibRich, y=fit)) +
    geom_smooth(method="lm", col="black")+
    ylab("Infection intensity (log BD load)") + 
    xlab("Dispersion (Distance from centroid)") +
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
save(mf_treat_with_inhibOTUs, file="mf_treat_with_inhibOTUs.Rdata")

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
save(mf_con_with_inhibOTUs, file="mf_con_with_inhibOTUs.RData")

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

pdf(file = "FIGURES/percent_inhib.pdf", width=13, height=5)
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


# TREAT

mf_treat_statdiff <- mf_treat_with_inhibOTUs %>%
    group_by(toadID, species, prepost, G) %>%
    summarize(meanReads=mean(reads), meanBd = mean(eBD_log)) %>%
    dplyr::select(-meanBd) %>%
    spread(key=prepost, value=meanReads) %>%
    # mutate(diff=Pos-Pre, significant = NA, p=NA) %>%
    mutate(fc=log((Pos-Pre)/Pre +1), significant = NA, p=NA) %>%
    mutate(temp= ifelse(is.finite(fc), fc, ifelse(is.infinite(fc), 10, 0))) %>%
    mutate(maxBD=NA)

#+ warning=FALSE, message=FALSE
allInhibTaxa <- unique(mf_treat_with_inhibOTUs$G)
allIndiv <- unique(mf_treat_with_inhibOTUs$toadID)
for ( inhib in allInhibTaxa) {
    # inhib <- allInhibTaxa[1]
    # inhib <- "Delftia"
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
        mf_treat_statdiff[which(mf_treat_statdiff$toadID==indiv & mf_treat_statdiff$G == inhib),"maxBD"] <- mf_treat_with_inhibOTUs %>%
            filter(G==inhib, toadID==indiv) %>%
            group_by(toadID) %>%
            summarise(maxBD=max(eBD_log, na.rm=TRUE)) %>%
            pull(maxBD)
    }
}

# 
# ####### Looking at specific OTUs #########
# 
# mf_con_statdiff <- mf_con_with_inhibOTUs %>%
#     group_by(toadID, species, prepost, G) %>%
#     summarize(meanReads=mean(reads), meanBd = mean(eBD_log)) %>%
#     dplyr::select(-meanBd) %>%
#     spread(key=prepost, value=meanReads) %>%
#     mutate(fc=log((Pos-Pre)/Pre +1), significant = NA, p=NA) %>%
#     mutate(temp= ifelse(is.finite(fc), fc, ifelse(fc>0, 10, -10)))
# 
#     
# allInhibTaxa <- unique(mf_con_with_inhibOTUs$G)
# allIndiv <- unique(mf_con_with_inhibOTUs$toadID)
# for ( inhib in allInhibTaxa) {
#     # inhib <- allInhibTaxa[1]
#     # inhib <- "Arthrobacter"
#     for ( indiv in allIndiv) {
#         if (exists("stat_temp")) {
#             remove(stat_temp)
#         }
#         # indiv <- allIndiv[1]
#         # indiv <- "Anbo_2"
#         mf_temp <- mf_con_with_inhibOTUs %>%
#             filter(G==inhib, toadID==indiv) %>%
#             group_by(toadID, G, prepost, time) %>%
#             summarise(reads=mean(reads)) %>%
#             spread(key=prepost, value=reads)
#         # stat_temp <- anova(lm(reads ~ prepost, data=mf_temp))
#         #+ warning=FALSE, message=FALSE
#         stat_temp <- wilcox.test(mf_temp$Pre, mf_temp$Pos)
#         mf_con_statdiff[which(mf_con_statdiff$toadID==indiv & mf_con_statdiff$G == inhib),"p"] <- (stat_temp$p.value)
#         mf_con_statdiff[which(mf_con_statdiff$toadID==indiv & mf_con_statdiff$G == inhib),"significant"] <- (stat_temp$p.value<0.05)
#         
#     }
# }
# 
# 
# # TREAT
# 
# mf_treat_statdiff <- mf_treat_with_inhibOTUs %>%
#     group_by(toadID, species, prepost, G) %>%
#     summarize(meanReads=mean(reads), meanBd = mean(eBD_log)) %>%
#     dplyr::select(-meanBd) %>%
#     spread(key=prepost, value=meanReads) %>%
#     # mutate(diff=Pos-Pre, significant = NA, p=NA) %>%
#     mutate(fc=log((Pos-Pre)/Pre +1), significant = NA, p=NA) %>%
#     mutate(temp= ifelse(is.finite(fc), fc, ifelse(fc>0, 10, -10)))
# 
# allInhibTaxa <- unique(mf_treat_with_inhibOTUs$G)
# allIndiv <- unique(mf_treat_with_inhibOTUs$toadID)
# for ( inhib in allInhibTaxa) {
#     # inhib <- allInhibTaxa[1]
#     # inhib <- "Arthrobacter"
#     for ( indiv in allIndiv) {
#         if (exists("stat_temp")) {
#             remove(stat_temp)
#         }
#         # indiv <- allIndiv[1]
#         # indiv <- "Anbo_2"
#         mf_temp <- mf_treat_with_inhibOTUs %>%
#             filter(G==inhib, toadID==indiv) %>%
#             group_by(toadID, G, prepost, time) %>%
#             summarise(reads=mean(reads)) %>%
#             spread(key=prepost, value=reads)
#         # stat_temp <- anova(lm(reads ~ prepost, data=mf_temp))
#         #+ warning=FALSE, message=FALSE
#         stat_temp <- wilcox.test(mf_temp$Pre, mf_temp$Pos)
#         mf_treat_statdiff[which(mf_treat_statdiff$toadID==indiv & mf_treat_statdiff$G == inhib),"p"] <- (stat_temp$p.value)
#         mf_treat_statdiff[which(mf_treat_statdiff$toadID==indiv & mf_treat_statdiff$G == inhib),"significant"] <- (stat_temp$p.value<0.05)
#         
#     }
# }

## Plot

mf_con_statdiff$Group <- "Control"
mf_treat_statdiff$Group <- "Treatment"
mf_combined_statdiff <- full_join(mf_con_statdiff, mf_treat_statdiff)

# PROPORTION OF SAMPLES THAT ARE SIGNIFICANTLY DIFFERENT
pdf("FIGURES/fold_change_inhib_histograms.pdf")
par(mfrow=c(1,2))
hist(mf_con_statdiff$fc, main="Control Fold Change", xlab="Fold Change", freq=FALSE, xlim=c(-7,7))
hist(mf_treat_statdiff$fc, main="Treat Fold Change", ylab="Fold Change", freq=FALSE, xlim=c(-7,7))
dev.off()


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

#### Inhibitory taxa

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

pdf("FIGURES/inhib_taxa_all.pdf", height=50, width=5)
mf_combo_with_inhibOTUs_G %>%
    ggplot(aes(x=time, y=genus_proportion, col=species), cex=0.1) +
    geom_point(alpha=0.3) +
    geom_vline(aes(xintercept=5.5), col="grey", lty=2) +
    facet_grid(G~TreatmentGroup, scales = "free")
dev.off()


# summary of all individuals and species counts
summary_individuals <- mf_combo_with_inhibOTUs_G %>%
    group_by(toadID, TreatmentGroup) %>%
    summarize(species=unique(species))
# how many controls in each group?
con_counts <- table(summary_individuals %>%
    filter(TreatmentGroup=="Control") %>%
    pull(species))
# how many treatment in each group?
treat_counts <- table(summary_individuals %>%
                        filter(TreatmentGroup=="Treatment") %>%
                        pull(species))

## Making numbered controls and treatment individuals for facetting
numbered_indiv <- mf_combo_with_inhibOTUs_G %>%
    group_by(toadID) %>%
    summarize(TreatmentGroup=unique(TreatmentGroup)) %>%
    separate(col=toadID, into=c("species","indiv"), remove = FALSE) %>%
    mutate(indiv = as.numeric(indiv)) %>%
    arrange(species, TreatmentGroup, indiv)
runlength <- numbered_indiv %>%
    pull(TreatmentGroup) %>%
    rle()
numbered_indiv$new_indiv <- sequence(runlength$lengths)

# Making custom legend
unique(mf_combo_with_inhibOTUs_G$G)
color_taxa <- c(
    "Pseudomonas"="purple"
    , "Acinetobacter"="orange"
    ,"Brevundimonas"="lightblue"             
    , "Rhizobium" ="brown"
    ,"Stenotrophomonas"="darkseagreen1"
    ,"Lysobacter" ="orangered"                 
    ,"Flavobacterium"="goldenrod"
    ,"Pedobacter"="darkmagenta"
    ,"Sphingobacterium" ="olivedrab1"          
    ,"Novosphingobium"="magenta"
    ,"Serratia"="orangered2"
    ,"Dyella"="mediumvioletred"                     
    ,"Terrimonas"="darkolivegreen4"
    ,"Elizabethkingia"="slateblue1"
    ,"Empedobacter"  ="grey"                
    , "Chryseobacterium"="forestgreen"
    ,"Arthrobacter"="saddlebrown"
    ,"Microbacterium" ="pink"               
    ,"Brevibacterium"="blue"
    ,"Ralstonia"="plum"
    ,"Burkholderiales_incertae_sedis"="gold"
    ,"Variovorax"="cyan3"
    ,"Mitsuaria"="darkslateblue"
    ,"Duganella" ="yellow"                    
    , "Janthinobacterium"="purple3"
    ,"Delftia"="green"
    ,"Chitinimonas"    ="wheat1"              
    ,"Lactococcus"="deepskyblue"
    ,"Bacillus"    ="mediumpurple1"  
)

### Control plotting
mf_combo_with_inhibOTUs_G$Indiv <- pull(numbered_indiv[match(mf_combo_with_inhibOTUs_G$toadID, numbered_indiv$toadID),"new_indiv"])
totalAbund <- mf_combo_with_inhibOTUs_G %>%
    group_by(SampleID,toadID, time) %>%
    summarize(totalAbund=sum(genus_proportion)) %>%
    ungroup()
mf_combo_with_inhibOTUs_G$totalAbund <- pull(totalAbund[match(mf_combo_with_inhibOTUs_G$SampleID, totalAbund$SampleID),"totalAbund"])
pdf("FIGURES/inhib_taxa_con.pdf", width=10, height=5)
mf_combo_with_inhibOTUs_G %>%
    filter(TreatmentGroup=="Control") %>%
    mutate(toadID = factor(toadID), time=factor(time), `Proportional Abundance`=genus_proportion, Time=time) %>%
    ggplot() +
    geom_bar(aes(x=Time, y=`Proportional Abundance`, fill=G), stat="identity", na.rm = FALSE, show.legend = TRUE) +
    geom_vline(aes(xintercept=5.5), col="grey", lty=2) +
    facet_grid(species~Indiv, scales = "free_y", drop = FALSE) +
    theme(
        strip.text.x = element_blank()
    ) +
    scale_fill_manual(values=color_taxa) +
    theme(legend.key.size = unit(0.25, "cm"), axis.text.x = element_text(size=5))

dev.off()

pdf("FIGURES/inhib_taxa_treat.pdf", width=12, height=5)
mf_combo_with_inhibOTUs_G %>%
    filter(TreatmentGroup=="Treatment") %>%
    mutate(toadID = factor(toadID), time=factor(time), `Proportional Abundance`=genus_proportion, Time=time) %>%
    ggplot() +
    geom_bar(aes(x=Time, y=`Proportional Abundance`, fill=G), stat="identity", na.rm = FALSE, show.legend = TRUE) +
    geom_vline(aes(xintercept=5.5), col="grey", lty=2) +
    facet_grid(species~Indiv, scales = "free_y", drop = FALSE) +
    theme(
        strip.text.x = element_blank()
    ) +
    scale_fill_manual(values=color_taxa)+
    scale_color_continuous(low="white",high="darkred")+
    theme(legend.key.size = unit(0.25, "cm"), axis.text.x = element_text(size=5)) +
    geom_point(aes(x=Time,y=(totalAbund+0.05), col=eBD_log), pch=8, cex=0.2)
dev.off()


pdf("FIGURES/combined_inhib_taxa_indiv.pdf",width=12, height=10 )
gg1 <- mf_combo_with_inhibOTUs_G %>%
    filter(TreatmentGroup=="Control") %>%
    mutate(toadID = factor(toadID), time=factor(time), `Proportional Abundance`=genus_proportion, Time=time) %>%
    ggplot() +
    geom_bar(aes(x=Time, y=`Proportional Abundance`, fill=G), stat="identity", na.rm = FALSE, show.legend = TRUE) +
    geom_vline(aes(xintercept=5.5), col="grey", lty=2) +
    facet_grid(species~Indiv, scales = "free_y", drop = FALSE) +
    # theme(
    #     strip.text.x = element_blank()
    # ) +
    scale_fill_manual(values=color_taxa) +
    theme(legend.key.size = unit(0.25, "cm"), axis.text.x = element_text(size=8,angle = 90, hjust = 1),panel.grid = element_blank(), strip.text.x = element_blank()) +
    labs(title="CONTROL INDIVIDUALS")
    # theme_classic()
gg2 <-  mf_combo_with_inhibOTUs_G %>%
    filter(TreatmentGroup=="Treatment") %>%
    mutate(eBD_log_zorona = ifelse(eBD_log==0, as.numeric(NA), eBD_log)) %>%
    mutate(toadID = factor(toadID), time=factor(time), `Proportional Abundance`=genus_proportion, Time=time) %>%
    ggplot() +
    geom_bar(aes(x=Time, y=`Proportional Abundance`, fill=G), stat="identity", na.rm = FALSE, show.legend = FALSE) +
    geom_vline(aes(xintercept=5.5), col="grey", lty=2) +
    facet_grid(species~Indiv, scales = "free_y", drop = FALSE) +
    # theme(
    #     strip.text.x = element_blank()
    #     , axis.text=element_text(size=8)
    # ) +
    scale_fill_manual(values=color_taxa)+
    theme(legend.key.size = unit(0.25, "cm"), axis.text.x = element_text(size=5,angle = 90, hjust = 1), panel.grid=element_blank(), strip.text.x = element_blank()) +
    geom_point(aes(x=Time,y=(totalAbund+0.05), col=eBD_log_zorona), pch=8, cex=0.2)+
    scale_color_continuous(low="white",high="darkred", na.value = NA)+
    labs(title="TREATMENT INDIVIDUALS")
grid.arrange(gg1
             ,gg2
             # , layout_matrix <- rbind(c(1,2,2))
             )
dev.off()



### INHIB: trying to identify inhibitory OTUs that increase after exposure ####

# # Get controll increased
# con_increased <- mf_con_statdiff %>%
#     filter(significant==TRUE) %>%
#     mutate(fc_change= ifelse(fc>0,TRUE,FALSE)) %>%
#     unite(speciesGfc, species,G,fc_change, remove=FALSE)
# # Get treatment increased, but subtract those also increased in control
# treat_increased_controlled <- mf_treat_statdiff %>%
#     mutate(fc_change= ifelse(fc>0,TRUE,FALSE)) %>%
#     unite(speciesGfc, species, G, fc_change, remove=FALSE) %>%
#     unite(toadIDGfc, toadID, G, fc_change, remove=FALSE) %>%
#     filter(significant==TRUE, speciesGfc %in% con_increased$speciesGfc, fc_change==TRUE) 

### SET COLORS 
color_taxa <- c(
    "Pseudomonas"="purple"
    , "Acinetobacter"="orange"
    ,"Brevundimonas"="lightblue"             
    , "Rhizobium" ="brown"
    ,"Stenotrophomonas"="darkseagreen1"
    ,"Lysobacter" ="orangered"                 
    ,"Flavobacterium"="goldenrod"
    ,"Pedobacter"="darkmagenta"
    ,"Sphingobacterium" ="olivedrab1"          
    ,"Novosphingobium"="magenta"
    ,"Serratia"="orangered2"
    ,"Dyella"="mediumvioletred"                     
    ,"Terrimonas"="darkolivegreen4"
    ,"Elizabethkingia"="slateblue1"
    ,"Empedobacter"  ="grey"                
    , "Chryseobacterium"="forestgreen"
    ,"Arthrobacter"="saddlebrown"
    ,"Microbacterium" ="pink"               
    ,"Brevibacterium"="blue"
    ,"Ralstonia"="plum"
    ,"Burkholderiales_incertae_sedis"="gold"
    ,"Variovorax"="cyan3"
    ,"Mitsuaria"="darkslateblue"
    ,"Duganella" ="yellow"                    
    , "Janthinobacterium"="purple3"
    ,"Delftia"="green"
    ,"Chitinimonas"    ="wheat1"              
    ,"Lactococcus"="deepskyblue"
    ,"Bacillus"    ="mediumpurple1"  
)

#### DO SELECT MICROBES CORRELATE WITH BD EXPOSURE, BY SPECIES??
dir.create("FIGURES/by_EXPOSURE")

mf_treat_with_inhibOTUs <- mf_treat_with_inhibOTUs %>%
    unite(speciesG, species, G, remove=FALSE)

all_tests_inhib <- as.data.frame(matrix(ncol=10, nrow=length(unique(mf_treat_with_inhibOTUs$speciesG)), dimnames = list(unique(mf_treat_with_inhibOTUs$speciesG),c("Species","G","change","changetype","p_exposure","p_stat","change_con","change_con_type","p_exposure_control","p_stat_control"))))
for ( n in unique(mf_treat_with_inhibOTUs$speciesG) ) {
    # n_spl <- c("Anbo","Acinetobacter")
    n_spl <- unlist(strsplit(n, split="_"))
    temp_dat <- mf_treat_with_inhibOTUs %>%
        filter(species==n_spl[1], G==n_spl[2]) %>%
        group_by(SampleID, species, time,toadID, eBD_log, PABD, prepost) %>%
        summarize(sumReads=sum(reads))
    temp_dat_con <- mf_con_with_inhibOTUs %>%
        filter(species==n_spl[1], G==n_spl[2]) %>%
        group_by(SampleID, species, time,toadID, eBD_log, PABD, prepost) %>%
        summarize(sumReads=sum(reads))
    
    if ( (nrow(temp_dat)>0) ) {
        if ( (var(temp_dat$sumReads)>0) & (length(unique(temp_dat$prepost))>1) ) {
            # Get "pre" reads
            a <- temp_dat %>%
                filter(prepost=="Pre") %>%
                group_by(toadID) %>%
                summarize(sumReads_bytoad=mean(sumReads)) %>%
                pull(sumReads_bytoad)
            # Get "post" reads
            b <- temp_dat %>%
                filter(prepost=="Pos") %>%
                group_by(toadID) %>%
                summarize(sumReads_bytoad=mean(sumReads)) %>%
                pull(sumReads_bytoad)
            # Need to compare to control; include control here
            a_c <- temp_dat_con %>%
                filter(prepost=="Pre") %>%
                group_by(toadID) %>%
                summarize(sumReads_bytoad=mean(sumReads)) %>%
                pull(sumReads_bytoad)
            b_c <- temp_dat_con %>%
                filter(prepost=="Pos") %>%
                group_by(toadID) %>%
                summarize(sumReads_bytoad=mean(sumReads)) %>%
                pull(sumReads_bytoad)

            perc_change <- c(mean(b)-mean(a))/mean(a)
            perc_change_con <- c(mean(b_c)-mean(a_c))/mean(a_c)
            
            temp_stat1 <- wilcox.test(b-a)
            temp_stat2 <- wilcox.test(b_c-a_c)

            all_tests_inhib[n,] <- c(n_spl[1]
                                     ,n_spl[2]
                                     ,(ifelse(perc_change==Inf, mean(b), perc_change))
                                     ,ifelse(perc_change==Inf, "Appeared", "Percent_change")
                                     ,temp_stat1$p.value
                                     ,temp_stat1$statistic
                                     ,(ifelse(is.na(perc_change_con),0, ifelse((perc_change_con==Inf),mean(b),perc_change_con)))
                                     ,ifelse(is.na(perc_change_con), "Absent", ifelse((perc_change_con==Inf),"Appeared","Percent_change"))
                                     ,temp_stat2$p.value
                                     ,temp_stat2$statistic
            )
        }
    }
}

# all_tests_inhib <- all_tests_inhib %>%
#     mutate(p_exposure=as.numeric(p_exposure)) %>%
#     filter(!is.na(Species),p_exposure<0.05 ) %>%
#     arrange(Species, p_exposure)
# What if we also removed controls that were significant?

# treat sig; con not
req1 <- (all_tests_inhib$p_exposure<0.05 & all_tests_inhib$p_exposure_control>0.05)
# treat and con are different directions
req2 <- (all_tests_inhib$p_exposure<0.05 & (sign(as.numeric(all_tests_inhib$change))!=sign(as.numeric(all_tests_inhib$change_con))) & (all_tests_inhib$changetype == all_tests_inhib$change_con_type))
# control is NA
req3 <- (all_tests_inhib$p_exposure<0.05 &all_tests_inhib$change_con_type=="Absent")

all_tests_inhib <- all_tests_inhib %>%
    mutate(p_exposure=as.numeric(p_exposure)) %>%
    filter(!is.na(Species),(req1|req2|req3) ) %>%
    arrange(Species, p_exposure)
write_csv(as.data.frame(all_tests_inhib), path="FIGURES/Exposure_corr_inhibAbundance_byspecies.csv")

all_test_exp_sp <- all_tests_inhib

all_sp <- all_tests_inhib %>%
    pull(Species) %>%
    unique()
for ( sp in all_sp ) {
    all_G <- all_tests_inhib %>%
        filter(Species==sp)%>%
        pull(G)
    for ( g in all_G) {
        assign(paste0("gg_",sp,"_bysp",g), mf_treat_with_inhibOTUs %>%
                   filter(species==sp,G==g) %>%
                   group_by(SampleID, species, toadID, time, eBD_log, PABD, prepost, G) %>%
                   summarize(totalReads=sum(reads)) %>%
                   group_by(G) %>%
                   mutate(maxheight=max(totalReads)) %>%
                   ungroup() %>%
                   mutate(eBD_log_zorona = ifelse(eBD_log==0, as.numeric(NA), eBD_log)) %>%
                   ggplot(aes(x=time, y=totalReads, fill=G)) +
                   geom_bar(stat="identity", show.legend = FALSE) +
                   geom_vline(aes(xintercept=5.5)) +
                   geom_point(aes(y=totalReads+maxheight*0.1, col=eBD_log_zorona), show.legend = FALSE) +
                   xlab("Time") +
                   ylab("Proportion of Reads")+
                   facet_grid(toadID~G) +
                   scale_fill_manual(values=color_taxa) +
                   scale_color_continuous(low="white",high="darkred", na.value = NA))
    }
}

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Anbo")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Anbo")%>%pull(toadID)))
pdf("FIGURES/by_EXPOSURE/Anbo_exposure_byspecies.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Anbo_byspAcinetobacter
             ,nrow=1)
dev.off()

# 
# num_G <- length(unique(all_tests_inhib%>%filter(Species=="Anma")%>%pull(G)))
# num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Anma")%>%pull(toadID)))
# pdf("FIGURES/by_EXPOSURE/Anma_exposure_byspecies.pdf", width=2.5*num_G, height=1.5*num_indiv)
# grid.arrange(gg_Anma_byspLactococcus+theme(strip.text.y = element_blank())
#              , gg_Anma_byspRalstonia+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
#              , gg_Anma_byspVariovorax+theme(axis.title.y=element_blank())
#              ,nrow=1)
# dev.off()

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Lica")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Lica")%>%pull(toadID)))
pdf("FIGURES/by_EXPOSURE/Lica_exposure_byspecies.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Lica_byspAcinetobacter+theme(strip.text.y = element_blank())
             , gg_Lica_byspDelftia+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Lica_byspBrevundimonas+theme(axis.title.y=element_blank())
             ,nrow=1)
dev.off()
# 
# num_G <- length(unique(all_tests_inhib%>%filter(Species=="Lipi")%>%pull(G)))
# num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Lipi")%>%pull(toadID)))
# pdf("FIGURES/by_EXPOSURE/Lipi_exposure_byspecies.pdf", width=2.5*num_G, height=1.5*num_indiv)
# grid.arrange(gg_Lipi_byspMitsuaria+theme(strip.text.y = element_blank())
#              , gg_Lipi_byspPseudomonas+theme(axis.title.y=element_blank())
#              ,nrow=1)
# dev.off()

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Osse")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Osse")%>%pull(toadID)))
pdf("FIGURES/by_EXPOSURE/Osse_exposure_byspecies.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Osse_byspAcinetobacter
             ,nrow=1)
dev.off()


#### DO SELECT MICROBES CORRELATE WITH BD LOAD, BY INDIVIDUAL??
mf_treat_with_inhibOTUs <- mf_treat_with_inhibOTUs %>%
    unite(toadIDG, toadID, G, remove=FALSE)

all_tests_inhib <- as.data.frame(matrix(ncol=11, nrow=length(unique(mf_treat_with_inhibOTUs$toadIDG)), dimnames = list(unique(mf_treat_with_inhibOTUs$toadIDG),c("Species","indiv","G","change","changetype","p_exposure","p_stat","change_con","change_con_type","p_exposure_control","p_stat_control"))))
for ( n in unique(mf_treat_with_inhibOTUs$toadIDG) ) {
    # n_spl <- c("Anbo",4,"Delftia")
    n_spl <- unlist(strsplit(n, split="_"))
    temp_dat <- mf_treat_with_inhibOTUs %>%
        separate(toadID, into=c("species","indiv"), remove=FALSE) %>%
        filter(species==n_spl[1], indiv==n_spl[2], G==n_spl[3]) %>%
        group_by(SampleID, species, time,toadID, eBD_log, PABD, prepost) %>%
        summarize(sumReads=sum(reads))
    ### I am comparing to SPECIES averages for individuals, since individuals are not consistent across treatments.
    temp_dat_con <- mf_con_with_inhibOTUs %>%
        separate(toadID, into=c("species","indiv"), remove=FALSE) %>%
        filter(species==n_spl[1], G==n_spl[3]) %>%
        group_by(SampleID, species, time,toadID, eBD_log, PABD, prepost) %>%
        summarize(sumReads=sum(reads))
    
    if ( (nrow(temp_dat)>0) ) {
        if ( (var(temp_dat$sumReads)>0) & (length(unique(temp_dat$prepost))>1) ) {
            a <- temp_dat %>%
                filter(prepost=="Pre") %>%
                pull(sumReads)
            b <- temp_dat %>%
                filter(prepost=="Pos") %>%
                pull(sumReads)
            a_c <- temp_dat_con %>%
                filter(prepost=="Pre") %>%
                pull(sumReads)
            b_c <- temp_dat_con %>%
                filter(prepost=="Pos") %>%
                pull(sumReads)
            
            perc_change<- (mean(b)-mean(a))/mean(a)
            perc_change_con<- (mean(b_c)-mean(a_c))/mean(a_c)
            
            temp_stat1 <- wilcox.test(a,b)
            temp_stat2 <- wilcox.test(a_c,b_c)
            
            all_tests_inhib[n,] <- c(n_spl[1]
                                     ,n_spl[2]
                                     ,n_spl[3]
                                     ,ifelse(perc_change==Inf, mean(b), perc_change)
                                     ,ifelse(perc_change==Inf, "Appear", "Percent_change")
                                     ,temp_stat1$p.value
                                     ,temp_stat1$statistic
                                     ,ifelse(is.na(perc_change_con),0,ifelse(perc_change_con==Inf, mean(b), perc_change_con))
                                     ,ifelse(is.na(perc_change_con),"Absent",ifelse(perc_change_con==Inf, "Appear", "Percent_change"))
                                     ,temp_stat2$p.value
                                     ,temp_stat2$statistic
            )
        }
    }
}
# all_tests_inhib <- all_tests_inhib %>%
#     mutate(p_exposure=as.numeric(p_exposure)) %>%
#     filter(!is.na(Species),p_exposure<0.05 ) %>%
#     arrange(Species, indiv, p_exposure) %>%
#     unite(toadIDG, Species,indiv,G, remove=FALSE)
# treat sig; con not
req1 <- (all_tests_inhib$p_exposure<0.05 & all_tests_inhib$p_exposure_control>0.05)
# treat and con are different directions
req2 <- (all_tests_inhib$p_exposure<0.05 & (sign(as.numeric(all_tests_inhib$change))!=sign(as.numeric(all_tests_inhib$change_con))) & (all_tests_inhib$changetype == all_tests_inhib$change_con_type))
# control is NA
req3 <- (all_tests_inhib$p_exposure<0.05 &all_tests_inhib$change_con_type=="Absent")

all_tests_inhib <- all_tests_inhib %>%
    mutate(p_exposure=as.numeric(p_exposure)) %>%
    filter(!is.na(Species),(req1|req2|req3) ) %>%
    arrange(Species, indiv, p_exposure) %>%
    unite(toadIDG, Species,indiv,G, remove=FALSE)
    
write_csv(as.data.frame(all_tests_inhib), path="FIGURES/Exposure_corr_inhibAbundance_bytoad.csv")

all_test_exp_indiv <- all_tests_inhib

all_sp <- all_tests_inhib %>%
    pull(Species) %>%
    unique()

for ( sp in all_sp ) {
    all_G <- all_tests_inhib %>%
        filter(Species==sp)%>%
        pull(unique(G))
    for ( g in all_G) {
        assign(paste0("gg_",sp,"_bytoad",g), mf_treat_with_inhibOTUs %>%
                   filter(species==sp,G==g) %>%
                   mutate(sig = ifelse(toadIDG %in% all_tests_inhib$toadIDG,5.5,NA))%>%
                   group_by(SampleID, species, toadID, time, eBD_log, PABD, prepost, G, sig) %>%
                   summarize(totalReads=sum(reads)) %>%
                   group_by(G) %>%
                   mutate(maxheight=max(totalReads)) %>%
                   ungroup() %>%
                   mutate(eBD_log_zorona = ifelse(eBD_log==0, as.numeric(NA), eBD_log)) %>%
                   ggplot(aes(x=time, y=totalReads, fill=G)) +
                   geom_bar(stat="identity", show.legend = FALSE) +
                   geom_vline(aes(xintercept=5.5)) +
                   geom_point(aes(y=totalReads+maxheight*0.1, col=eBD_log_zorona), show.legend = FALSE) +
                   geom_point(aes(x=sig, y=max(totalReads)/2), col="red",cex=2, pch=8, show.legend = FALSE)+
                   xlab("Time") +
                   ylab("Proportion of Reads")+
                   facet_grid(toadID~G) +
                   scale_fill_manual(values=color_taxa) +
                   scale_color_continuous(low="white",high="darkred", na.value = NA))
    }
}


num_G <- length(unique(all_tests_inhib%>%filter(Species=="Anbo")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Anbo")%>%pull(toadID)))
pdf("FIGURES/by_EXPOSURE/Anbo_exposure_bytoadID.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Anbo_bytoadAcinetobacter+theme(strip.text.y = element_blank())
             ,gg_Anbo_bytoadBrevundimonas+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             ,gg_Anbo_bytoadDelftia+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             ,gg_Anbo_bytoadPedobacter+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             ,gg_Anbo_bytoadVariovorax+theme(axis.title.y=element_blank())
             ,nrow=1)
dev.off()

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Anma")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Anma")%>%pull(toadID)))
pdf("FIGURES/by_EXPOSURE/Anma_exposure_bytoadID.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Anma_bytoadBrevundimonas+theme(strip.text.y = element_blank())
             , gg_Anma_bytoadFlavobacterium+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Anma_bytoadMitsuaria+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Anma_bytoadSphingobacterium+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Anma_bytoadTerrimonas+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Anma_bytoadVariovorax+theme(axis.title.y=element_blank())
             ,nrow=1)
dev.off()

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Lica")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Lica")%>%pull(toadID)))
pdf("FIGURES/by_EXPOSURE/Lica_exposure_bytoadID.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Lica_bytoadAcinetobacter+theme(strip.text.y = element_blank())
             , gg_Lica_bytoadDuganella+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Lica_bytoadJanthinobacterium+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Lica_bytoadMicrobacterium+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Lica_bytoadMitsuaria+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Lica_bytoadPseudomonas+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Lica_bytoadStenotrophomonas+theme(axis.title.y=element_blank())
             ,nrow=1)
dev.off()
# 
# num_G <- length(unique(all_tests_inhib%>%filter(Species=="Lipi")%>%pull(G)))
# num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Lipi")%>%pull(toadID)))
# pdf("FIGURES/by_EXPOSURE/Lipi_exposure_bytoadID.pdf", width=2.5*num_G, height=1.5*num_indiv)
# grid.arrange(gg_Lipi_bytoadMitsuaria+theme(strip.text.y = element_blank())
#              , gg_Lipi_bytoadPseudomonas+theme(axis.title.y=element_blank())
#              ,nrow=1)
# dev.off()

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Osse")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Osse")%>%pull(toadID)))
pdf("FIGURES/by_EXPOSURE/Osse_exposure_bytoadID.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Osse_bytoadBrevundimonas+theme(strip.text.y = element_blank())
             , gg_Osse_bytoadElizabethkingia+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Osse_bytoadNovosphingobium+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Osse_bytoadPseudomonas+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Osse_bytoadRhizobium+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Osse_bytoadStenotrophomonas+theme(axis.title.y=element_blank())
             ,nrow=1)
dev.off()


### INHIB: trying to identify OTUs that increase with infection ####

dir.create("FIGURES/by_BDLOAD")
#### DO ANY MICROBES CORRELATE WITH BD LOAD, BY SPECIES??
mf_treat_with_inhibOTUs <- mf_treat_with_inhibOTUs %>%
    unite(speciesG, species, G, remove=FALSE)


# all_tests_inhib <- as.data.frame(matrix(ncol=6, nrow=length(unique(mf_treat_with_inhibOTUs$speciesG)), dimnames = list(unique(mf_treat_with_inhibOTUs$speciesG),c("Species","G","direction","kendall_sumReads","kendall_sumReads_nozeros","wilcox_PABD"))))
all_tests_inhib <- as.data.frame(matrix(ncol=6, nrow=length(unique(mf_treat_with_inhibOTUs$speciesG)), dimnames = list(unique(mf_treat_with_inhibOTUs$speciesG),c("Species","G","tau","kendall_sumReads","kendall_sumReads_nozeros","wilcox_PABD"))))
for ( n in unique(mf_treat_with_inhibOTUs$speciesG) ) {
    # n <- unique(mf_treat_with_inhibOTUs$speciesG)[10]
    n_spl <- unlist(strsplit(n, split="_"))
    reads_before <- mf_treat_with_inhibOTUs %>%
                             filter(species==n_spl[1], G==n_spl[2], prepost=="Pre") %>%
                             group_by(toadID, time) %>%
                             summarize(before_reads=sum(reads)) %>%
                             group_by(toadID) %>%
        summarize(ave_before_reads = median(before_reads))
    temp_dat <- mf_treat_with_inhibOTUs %>%
        left_join(reads_before) %>%
        filter(species==n_spl[1], G==n_spl[2], prepost=="Pos") %>%
        group_by(SampleID, species, time,toadID, eBD_log, PABD, prepost, ave_before_reads) %>%
        summarize(sumReads=(sum(reads))) %>%
        mutate(sumReads_adj=sumReads-ave_before_reads)
    # Create a standardized value
    
    
    
    if ( (nrow(temp_dat)>0) ) {
        if ( (var(temp_dat$eBD_log)>0) & (var(temp_dat$sumReads_adj)>0) ) {
            temp_stat1 <- cor.test(temp_dat$eBD_log, temp_dat$sumReads_adj, method="kendall")
            temp_dat2 <- temp_dat %>%
                filter(eBD_log>0)
            temp_stat2 <- cor.test(temp_dat2$eBD_log, temp_dat2$sumReads_adj, method="kendall")
            temp_dat1a <- temp_dat %>%
                filter(eBD_log==0)
            temp_stat3 <- wilcox.test(temp_dat1a$sumReads_adj, temp_dat2$sumReads_adj)
            
            all_tests_inhib[n,] <- c(n_spl[1]
                                     ,n_spl[2]
                                     ,temp_stat1$estimate
                                     ,temp_stat1$p.value
                                     ,temp_stat2$p.value
                                     ,temp_stat3$p.value)
            
        } 
        
    }
}
# all_tests_inhib <- all_tests_inhib %>%
#     mutate(kendall_sumReads=as.numeric(kendall_sumReads),kendall_sumReads_nozeros=as.numeric(kendall_sumReads_nozeros),wilcox_PABD=as.numeric(wilcox_PABD) ) %>%
#     filter(!is.na(Species),!(kendall_sumReads>0.05 & kendall_sumReads_nozeros >0.05 & wilcox_PABD>0.05) ) %>%
#     unite(SpeciesG, Species, G, remove=FALSE) %>%
#     arrange(Species, kendall_sumReads)
all_tests_inhib <- all_tests_inhib %>%
    mutate(kendall_sumReads=as.numeric(kendall_sumReads),kendall_sumReads_nozeros=as.numeric(kendall_sumReads_nozeros),wilcox_PABD=as.numeric(wilcox_PABD) ) %>%
    filter(!is.na(Species), (kendall_sumReads_nozeros<0.05| wilcox_PABD<0.05) ) %>%
    unite(SpeciesG, Species, G, remove=FALSE) %>%
    arrange(Species, kendall_sumReads)
write_csv(as.data.frame(all_tests_inhib), path="FIGURES/BD_corr_inhibAbundance_byspecies_ALLASVS.csv")

all_test_bd_sp <- all_tests_inhib

all_sp <- all_tests_inhib %>%
    pull(Species) %>%
    unique()
for ( sp in all_sp ) {
    all_G <- all_tests_inhib %>%
        filter(Species==sp)%>%
        pull(G)
    for ( g in all_G) {
        assign(paste0("gg_",sp,"_bysp_bd",g), mf_treat_with_inhibOTUs %>%
                   mutate(sig_sumReads=all_tests_inhib[match(speciesG,all_tests_inhib$SpeciesG),"kendall_sumReads"]
                          , sig_sumReads_nozeros=all_tests_inhib[match(speciesG,all_tests_inhib$SpeciesG),"kendall_sumReads_nozeros"]
                          , sig_wilcox_PABD=all_tests_inhib[match(speciesG,all_tests_inhib$SpeciesG),"wilcox_PABD"]) %>%
                   mutate(sig = ifelse((sig_sumReads<sig_sumReads_nozeros)&(sig_sumReads<sig_wilcox_PABD)&(sig_sumReads<0.05),"sumReads"
                                       , ifelse((sig_sumReads_nozeros<sig_wilcox_PABD)&(sig_sumReads_nozeros<0.05), "sumReads_nozeros"
                                                , ifelse(sig_wilcox_PABD<0.05,"wilcox_PABD","notSig")))) %>%
                   filter(species==sp,G==g, prepost=="Pos") %>%
                   group_by(SampleID, species, toadID, time, eBD_log, PABD, sig, G) %>%
                   summarize(totalReads=sum(reads)) %>%
                   unite(speciesG, species, G, remove=FALSE) %>%
                   ggplot(aes(x=eBD_log, y=totalReads)) +
                   geom_point(aes(fill=ifelse(sig=="sumReads",'red',ifelse(sig=="sumReads_nozeros",'blue',ifelse(sig=="wilcox_PABD",'green','black')))),cex=2, pch=21) +
                   scale_fill_identity() +
                   facet_grid(toadID~G, drop = FALSE) +
                   xlab("log BD load") +
                   ylab("Proportion of Reads")
        )
    }
}



### KENDALL RESULTS
num_G <- length(unique(all_tests_inhib%>%filter(Species=="Anbo")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Anbo")%>%pull(toadID)))
pdf("FIGURES/by_BDLOAD/Anbo_BDload_byspecies.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Anbo_bysp_bdElizabethkingia+theme(strip.text.y = element_blank())
             , gg_Anbo_bysp_bdFlavobacterium+theme(axis.title.y=element_blank())
             ,nrow=1)
dev.off()

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Anma")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Anma")%>%pull(toadID)))
pdf("FIGURES/by_BDLOAD/Anma_BDload_byspecies.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Anma_bysp_bdMicrobacterium+theme(strip.text.y = element_blank())
             , gg_Anma_bysp_bdAcinetobacter+theme(axis.title.y=element_blank())
             ,nrow=1)
dev.off()

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Lica")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Lica")%>%pull(toadID)))
pdf("FIGURES/by_BDLOAD/Lica_BDload_byspecies.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Lica_bysp_bdLysobacter+theme(strip.text.y = element_blank())
             , gg_Lica_bysp_bdVariovorax+theme(axis.title.y=element_blank())
             ,nrow=1)
dev.off()

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Osse")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Osse")%>%pull(toadID)))
pdf("FIGURES/by_BDLOAD/Osse_BDload_byspecies.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Osse_bysp_bdMitsuaria
             ,nrow=1)
dev.off()

## TRY 2 at plotting

all_sp <- all_tests_inhib %>%
    pull(Species) %>%
    unique()
for ( sp in all_sp ) {
    all_G <- all_tests_inhib %>%
        filter(Species==sp)%>%
        pull(G)
    for ( g in all_G) {
        assign(paste0("gg_",sp,"_bysp2_bd",g), mf_treat_with_inhibOTUs %>%
                   filter(species==sp,G==g) %>%
                   group_by(SampleID, species, toadID, time, eBD_log, PABD, prepost, G) %>%
                   summarize(totalReads=sum(reads)) %>%
                   group_by(G) %>%
                   mutate(maxheight=max(totalReads)) %>%
                   ungroup() %>%
                   mutate(eBD_log_zorona = ifelse(eBD_log==0, as.numeric(NA), eBD_log)) %>%
                   ggplot(aes(x=time, y=totalReads, group=toadID, fill=eBD_log)) +
                   geom_bar(stat="identity",color="black", show.legend = FALSE) +
                   geom_vline(aes(xintercept=5.5)) +
                   # geom_point(aes(x=time-0.125, y=totalReads+maxheight*0.01, col=eBD_log_zorona), cex=2, show.legend = FALSE) +
                   xlab("Time") +
                   ylab("Proportion of Reads")+
                   facet_grid(species~G) +
                   # scale_fill_manual(values=color_taxa) +
                   # scale_color_continuous(low="white",high="darkred", na.value = NA) 
                   scale_fill_continuous(low="white",high="darkred", na.value = NA)
        )
    }
}

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Anbo")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Anbo")%>%pull(toadID)))
pdf("FIGURES/by_BDLOAD/Anbo_BDload_byspecies2.pdf", width=2.5*num_G, height=4)
grid.arrange(gg_Anbo_bysp2_bdAcinetobacter+theme(strip.text.y = element_blank())
             , gg_Anbo_bysp2_bdElizabethkingia+theme(axis.title.y=element_blank())
             ,nrow=1)
dev.off()

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Anma")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Anma")%>%pull(toadID)))
pdf("FIGURES/by_BDLOAD/Anma_BDload_byspecies2.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Anma_bysp2_bdMicrobacterium
             ,nrow=1)
dev.off()

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Lica")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Lica")%>%pull(toadID)))
pdf("FIGURES/by_BDLOAD/Lica_BDload_byspecies2.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Lica_bysp2_bdLysobacter+theme(strip.text.y = element_blank())
             , gg_Lica_bysp2_bdVariovorax+theme(axis.title.y=element_blank())
             ,nrow=1)
dev.off()

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Osse")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Osse")%>%pull(toadID)))
pdf("FIGURES/by_BDLOAD/Osse_BDload_byspecies2.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Osse_bysp2_bdElizabethkingia
             ,nrow=1)
dev.off()


####### Summary of all
all_G <- unique(all_tests_inhib %>%
                    pull(G))
all_sp <- all_tests_inhib %>%
    pull(Species) %>%
    unique()
num <- 1
for ( g in all_G) {
    for ( sp in all_sp ) {
        if ( (match(sp, all_sp) == 1) & (match(g, all_G) == 1) ) { # top left
            element_striptextx = element_text()
            element_striptexty = element_blank()
            element_axisy = element_blank()
            element_axisx = element_blank()
        } else if ((match(g, all_G) == length(all_G)) & (match(sp, all_sp) == 1)) { # top right
            element_striptextx = element_text()
            element_striptexty = element_text()
            element_axisy = element_blank()
            element_axisx = element_blank()
        } else if ( (match(g, all_G) == length(all_G)) & (match(sp, all_sp) == length(all_sp)) ) { # bottom right
            element_striptextx = element_blank()
            element_striptexty = element_text()
            element_axisy = element_blank()
            element_axisx = element_text()
        } else if ( (match(g, all_G) == 1) & (match(sp, all_sp) == length(all_sp)) ) { # bottom left
            element_striptextx = element_blank()
            element_striptexty = element_blank()
            element_axisy = element_blank()
            element_axisx = element_blank()
        } else if ( (match(g, all_G) == 1) & !(match(sp, all_sp) %in% c(1, length(all_sp))) ) { # left side
            element_striptextx = element_blank()
            element_striptexty = element_blank()
            element_axisy = element_blank()
            element_axisx = element_blank()
        } else if ( !(match(g, all_G) %in% c(1, length(all_G))) & (match(sp, all_sp) ==length(all_sp)) ) { # bottom side
            element_striptextx = element_blank()
            element_striptexty = element_blank()
            element_axisy = element_blank()
            element_axisx = element_blank()
        } else if ( !(match(g, all_G) %in% c(1, length(all_G))) & (match(sp, all_sp) ==1) ) { # top side
            element_striptextx = element_text()
            element_striptexty = element_blank()
            element_axisy = element_blank()
            element_axisx = element_blank()
        } else if ( (match(g, all_G) == length(all_G)) & !(match(sp, all_sp) %in% c(1, length(all_sp))) ) { # right side
            element_striptextx = element_blank()
            element_striptexty = element_text()
            element_axisy = element_blank()
            element_axisx = element_blank()
        } else {
            element_striptextx = element_blank()
            element_striptexty = element_blank()
            element_axisy = element_blank()
            element_axisx = element_blank()
        }
        assign(paste0("gg_",num),mf_treat_with_inhibOTUs %>%
                   filter(species==sp, G==g) %>%
                   group_by(SampleID, species, toadID, time, eBD_log, PABD, prepost, G) %>%
                   summarize(totalReads=sum(reads)) %>%
                   group_by(G) %>%
                   mutate(maxheight=max(totalReads)) %>%
                   ungroup() %>%
                   mutate(eBD_log_zorona = ifelse(eBD_log==0, as.numeric(NA), eBD_log)) %>%
                   ggplot(aes(x=time, y=totalReads, group=toadID, fill=eBD_log)) +
                   geom_bar(stat="identity",color="black", show.legend = FALSE) +
                   geom_vline(aes(xintercept=5.5)) +
                   # geom_point(aes(x=time-0.125, y=totalReads+maxheight*0.01, col=eBD_log_zorona), cex=2, show.legend = FALSE) +
                   xlab("Time") +
                   ylab("Proportion of Reads")+
                   facet_grid(species~G, scales = "free") +
                   # scale_fill_manual(values=color_taxa) +
                   # scale_color_continuous(low="white",high="darkred", na.value = NA) 
                   scale_fill_continuous(low="white",high="darkred", na.value = NA)+
                   theme(strip.text.y = element_striptexty
                         , axis.title.y=element_axisy
                         , axis.title.x=element_axisx
                         , strip.text.x=element_striptextx)
        )
        num <- num+1
    }
}
# ### PEARSON LAYOUT
# lay <- rbind(c(1,5,9,13,17,21,25),
#              c(2,6,10,14,18,22,26),
#              c(3,7,11,15,19,23,27),
#              c(4,8,12,16,20,24,28))
# 
# pdf("FIGURES/by_BDLOAD/bd_corr_bygenus_stacked.pdf",width=2.5*7, height=1.5*4)
# grid.arrange(gg_1, gg_2, gg_3, gg_4, gg_5, gg_6, gg_7
#              , gg_8, gg_9, gg_10, gg_11, gg_12, gg_13, gg_14
#              , gg_15, gg_16, gg_17, gg_18, gg_19, gg_20, gg_21
#              , gg_22, gg_23, gg_24, gg_25, gg_26, gg_27, gg_28
#              , layout_matrix=lay
# )
# dev.off()

#### KENDALL LAYOUT
lay <- rbind(c(1,5,9,13,17),
             c(2,6,10,14,18),
             c(3,7,11,15,19),
             c(4,8,12,16,20))

pdf("FIGURES/by_BDLOAD/bd_corr_bygenus_stacked.pdf",width=2.5*7, height=1.5*4)
grid.arrange(gg_1, gg_2, gg_3, gg_4, gg_5, gg_6, gg_7
             , gg_8, gg_9, gg_10, gg_11, gg_12, gg_13, gg_14
             , gg_15, gg_16, gg_17, gg_18, gg_19, gg_20
             , layout_matrix=lay, left="Proportion of Reads", bottom="Time"
)
dev.off()

######

#### DO ANY MICROBES CORRELATE WITH BD LOAD, BY INDIVIDUAL??
mf_treat_with_inhibOTUs <- mf_treat_with_inhibOTUs %>%
    unite(toadIDG, toadID, G, remove=FALSE)

# Individual counts
sort(unique(mf_treat_with_inhibOTUs$toadID))

all_tests_inhib <- as.data.frame(matrix(ncol=7, nrow=length(unique(mf_treat_with_inhibOTUs$toadIDG)), dimnames = list(unique(mf_treat_with_inhibOTUs$toadIDG),c("Species","indiv","G","tau","kendall_sumReads","kendall_sumReads_nozeros","wilcox_PABD"))))
for ( n in unique(mf_treat_with_inhibOTUs$toadIDG) ) {
    # n <- unique(mf_treat_with_inhibOTUs$toadIDG)[1]
    # n_spl <- c("Anbo","4","Delftia")
    n_spl <- unlist(strsplit(n, split="_"))
    temp_dat <- mf_treat_with_inhibOTUs %>%
        separate(toadID, into=c("species","indiv"), remove=FALSE) %>%
        filter(species==n_spl[1], indiv==n_spl[2],G==n_spl[3], prepost=="Pos") %>%
        group_by(SampleID, species, time,toadID, eBD_log, PABD, prepost) %>%
        summarize(sumReads=sum(reads))
    
    if ( (nrow(temp_dat)>0) ) {
        if ( (var(temp_dat$eBD_log)>0) & (var(temp_dat$sumReads)>0) ) {
            temp_stat1 <- cor.test(temp_dat$eBD_log, temp_dat$sumReads, method="kendall")
            temp_dat2 <- temp_dat %>%
                filter(eBD_log>0)
            if ( nrow(temp_dat2)>2) {
                temp_stat2 <- cor.test(temp_dat2$eBD_log, temp_dat2$sumReads, method="kendall")
            } else {
                temp_stat2 <- list(p.value=NA)
            }
            temp_dat1a <- temp_dat %>%
                filter(eBD_log==0)
            if ( nrow(temp_dat2) > 0 & nrow(temp_dat1a)>0) {
                temp_stat3 <- wilcox.test(temp_dat1a$sumReads, temp_dat2$sumReads)
            } else {
                temp_stat3 <- list(p.value=NA)
            }
            
            all_tests_inhib[n,] <- c(n_spl[1]
                                     ,n_spl[2]
                                     ,n_spl[3]
                                     ,temp_stat1$estimate
                                     ,temp_stat1$p.value
                                     ,temp_stat2$p.value
                                     ,temp_stat3$p.value)
            
        } 
        
    }
}
# all_tests_inhib <- all_tests_inhib %>%
#     mutate(kendall_sumReads=as.numeric(kendall_sumReads),kendall_sumReads_nozeros=as.numeric(kendall_sumReads_nozeros),wilcox_PABD=as.numeric(wilcox_PABD) ) %>%
#     filter(!is.na(Species),!(kendall_sumReads>0.05 & kendall_sumReads_nozeros >0.05 & wilcox_PABD>0.05) ) %>%
#     unite(toadIDG, Species, indiv, G, remove=FALSE) %>%
#     arrange(Species, indiv, kendall_sumReads)
all_tests_inhib <- all_tests_inhib %>%
    mutate(kendall_sumReads=as.numeric(kendall_sumReads),kendall_sumReads_nozeros=as.numeric(kendall_sumReads_nozeros),wilcox_PABD=as.numeric(wilcox_PABD) ) %>%
    filter(!is.na(Species),(kendall_sumReads_nozeros <0.05 | wilcox_PABD<0.05) ) %>%
    unite(toadIDG, Species, indiv, G, remove=FALSE) %>%
    arrange(Species, indiv, kendall_sumReads)
write_csv(as.data.frame(all_tests_inhib), path="FIGURES/BD_corr_inhibAbundance_bytoad_ALLASVS.csv")

all_test_bd_indiv <- all_tests_inhib


all_sp <- all_tests_inhib %>%
    pull(Species) %>%
    unique()
for ( sp in all_sp ) {
    all_G <- all_tests_inhib %>%
        filter(Species==sp)%>%
        pull(G)
    for ( g in all_G) {
        assign(paste0("gg_",sp,"_bytoadID_bd",g), mf_treat_with_inhibOTUs %>%
                   mutate(sig_sumReads=all_tests_inhib[match(toadIDG,all_tests_inhib$toadIDG),"kendall_sumReads"]
                          , sig_sumReads_nozeros=all_tests_inhib[match(toadIDG,all_tests_inhib$toadIDG),"kendall_sumReads_nozeros"]
                          , sig_wilcox_PABD=all_tests_inhib[match(toadIDG,all_tests_inhib$toadIDG),"wilcox_PABD"]) %>%
                   mutate(sig_sumReads=ifelse(is.na(sig_sumReads), 1,sig_sumReads)
                          , sig_sumReads_nozeros=ifelse(is.na(sig_sumReads_nozeros), 1,sig_sumReads_nozeros)
                          , sig_wilcox_PABD=ifelse(is.na(sig_wilcox_PABD), 1,sig_wilcox_PABD)) %>%
                   mutate(sig = ifelse((sig_sumReads<sig_sumReads_nozeros)&(sig_sumReads<sig_wilcox_PABD)&(sig_sumReads<0.05),"sumReads"
                                       , ifelse((sig_sumReads_nozeros<sig_wilcox_PABD)&(sig_sumReads_nozeros<0.05), "sumReads_nozeros"
                                                , ifelse(sig_wilcox_PABD<0.05,"wilcox_PABD","notSig")))) %>%
                   filter(species==sp,G==g, prepost=="Pos") %>%
                   group_by(SampleID, species, toadID, time, eBD_log, PABD, sig, G) %>%
                   summarize(totalReads=sum(reads)) %>%
                   unite(toadIDG, species, G, remove=FALSE) %>%
                   ggplot(aes(x=eBD_log, y=totalReads)) +
                   geom_point(aes(fill=ifelse(sig=="sumReads",'red',ifelse(sig=="sumReads_nozeros",'blue',ifelse(sig=="wilcox_PABD",'green','black')))),cex=2, pch=21) +
                   scale_fill_identity() +
                   facet_grid(toadID~G, drop = FALSE) +
                   xlab("log BD load") +
                   ylab("Proportion of Reads")
        )
    }
}


### KENDALL
num_G <- length(unique(all_tests_inhib%>%filter(Species=="Anbo")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Anbo")%>%pull(toadID)))
pdf("FIGURES/by_BDLOAD//Anbo_BDload_bytoadID.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Anbo_bytoadID_bdPedobacter+theme(strip.text.y = element_blank())
             , gg_Anbo_bytoadID_bdStenotrophomonas+theme(axis.title.y=element_blank())
             ,nrow=1)
dev.off()

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Anma")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Anma")%>%pull(toadID)))
pdf("FIGURES/by_BDLOAD/Anma_BDload_bytoadID.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Anma_bytoadID_bdMicrobacterium+theme(strip.text.y = element_blank())
             , gg_Anma_bytoadID_bdMitsuaria+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Anma_bytoadID_bdNovosphingobium+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Anma_bytoadID_bdRhizobium+theme(strip.text.y = element_blank(), axis.title.y=element_blank())
             , gg_Anma_bytoadID_bdVariovorax+theme(axis.title.y=element_blank())
             ,nrow=1)
dev.off()

num_G <- length(unique(all_tests_inhib%>%filter(Species=="Lica")%>%pull(G)))
num_indiv <- length(unique(mf_treat_with_inhibOTUs%>%filter(species=="Lica")%>%pull(toadID)))
pdf("FIGURES/by_BDLOAD/Lica_BDload_bytoadID.pdf", width=2.5*num_G, height=1.5*num_indiv)
grid.arrange(gg_Lica_bytoadID_bdVariovorax
             ,nrow=1)
dev.off()



####### Join all exposure and BD corr data for heatmap ######
all_G <-sort(unique(c(all_test_exp_sp$G
               ,all_test_exp_indiv$G
               ,all_test_bd_sp$G
               ,all_test_bd_indiv$G)))
all_sp_reps <- mf_treat_without_init_infect %>%
    group_by(toadID, species) %>%
    summarize() %>%
    pull(species)%>%
    table()

temp_g <- rep(all_G, each=sum(all_sp_reps))
temp_sp <- rep(mf_treat_without_init_infect %>%
    group_by(toadID, species) %>%
    summarize() %>%
    pull(species), length(all_G))
temp_indiv <- rep(mf_treat_without_init_infect %>%
    group_by(toadID, species) %>%
    summarize() %>%
    pull(toadID), length(all_G))
row_names_expand <- paste(temp_g, as.character(temp_sp), temp_indiv, sep = "_" )

heatmap_inhib <- as.data.frame(matrix(ncol=7, nrow=length(row_names_expand), dimnames = list(row_names_expand, c("temp_g","temp_sp","temp_indiv","Exposure_by_species","Exposure_by_indiv","BDcorr_by_species","BDcorr_by_indiv"))))
heatmap_inhib[,c("temp_g","temp_sp","temp_indiv")] <- cbind(temp_g, as.character(temp_sp), temp_indiv)
for ( g in all_G ) {
    # By species; + or -
    exp_sp <- all_test_exp_sp %>%
        filter(G == g) %>%
        mutate(dir=sign(as.numeric(change))) %>%
        mutate(dir = ifelse(is.na(dir), +1,dir)) %>%
        dplyr::select(Species, dir)
    bd_sp <- all_test_bd_sp %>%
        filter(G == g) %>%
        mutate(dir=sign(as.numeric(tau))) %>%
        dplyr::select(Species, dir)
    
    # By indiv; + or -
    exp_indiv <- all_test_exp_indiv %>%
        filter(G == g) %>%
        mutate(dir=sign(as.numeric(change))) %>%
        mutate(dir = ifelse(is.na(dir), +1,dir)) %>%
        unite(toadID, Species, indiv, remove=FALSE) %>%
        dplyr::select(Species,indiv, toadID,dir)
    bd_indiv <- all_test_bd_indiv %>%
        filter(G == g) %>%
        mutate(dir=sign(as.numeric(tau))) %>%
        unite(toadID, Species, indiv, remove=FALSE) %>%
        dplyr::select(Species,indiv, toadID,dir)
    
    ## Exp sp column
    if (nrow(exp_sp)>0) {
        for ( r in 1:nrow(exp_sp)) {
            sp <- exp_sp[r,"Species"]
            heatmap_inhib[which((heatmap_inhib$temp_sp==sp) & (heatmap_inhib$temp_g==g)),"Exposure_by_species"] <- exp_sp[r,"dir"]
            
        }
    }
    
    ## BD sp column
    if (nrow(bd_sp)>0) {
        for ( r in 1:nrow(bd_sp)) {
            sp <- bd_sp[r,"Species"]
            heatmap_inhib[which((heatmap_inhib$temp_sp==sp) & (heatmap_inhib$temp_g==g)),"BDcorr_by_species"] <- bd_sp[r,"dir"]
            
        }
    }
    
    
    ## Exp indiv column
    if (nrow(exp_indiv)>0) {
        for ( r in 1:nrow(exp_indiv)) {
            indiv <- exp_indiv[r,"toadID"]
            heatmap_inhib[which((heatmap_inhib$temp_indiv==indiv) & (heatmap_inhib$temp_g==g)),"Exposure_by_indiv"] <- exp_indiv[r,"dir"]
            
        }
    }
    
    ## BDcorr indiv column
    if (nrow(bd_indiv)>0) {
        for ( r in 1:nrow(bd_indiv)) {
            indiv <- bd_indiv[r,"toadID"]
            heatmap_inhib[which((heatmap_inhib$temp_indv==indiv) & (heatmap_inhib$temp_g==g)),"BDcorr_by_indiv"] <- bd_indiv[r,"dir"]
            
        }
    }
    
    
}

heatmap_inhib

#### turn this into barplot

# TEST <- heatmap_inhib %>%
#     unite(g_sp, temp_g, temp_sp, remove=FALSE) %>%
#     unite(g_sp_indiv, temp_g, temp_sp, temp_indiv, remove=FALSE) %>%
#     # dplyr::select(-c("temp_indiv")) %>%
#     # mutate(single=rownames(heatmap_inhib)) %>%
#     gather(key=Category,value=value)

# Get ggplot default colors
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

sp_col <- gg_color_hue(5)
names(sp_col) <- c("Anbo","Anma","Lica","Lipi","Osse")

all_colors_heatmap <- c(color_taxa,sp_col, "1"="blue","-1"="red", "+"="blue","-"="red")

### LONG FORMAT
heatmap_inhib_long <- heatmap_inhib %>%
    mutate(allRows=rownames(heatmap_inhib)) %>%
    arrange(temp_g, temp_sp, temp_indiv) %>%
    mutate(Amphib_Species=temp_sp, `Bacterial Genus`=temp_g, `Change after exposure (Species)`=Exposure_by_species, `Change after exposure (Individual)`=Exposure_by_indiv, `Correlation w BD load (Species)`=BDcorr_by_species, `Correlation w BD load (Indiv)`=BDcorr_by_indiv) %>%
    dplyr::select(allRows,temp_g ,temp_sp,`Bacterial Genus`, Amphib_Species,`Change after exposure (Species)`,`Change after exposure (Individual)`,`Correlation w BD load (Species)`, `Correlation w BD load (Indiv)`) %>%
    gather(key=Test_category, value=value,4:9) %>%
    mutate(Test_category=factor(Test_category, levels=c("Bacterial Genus","Amphib_Species","Change after exposure (Species)","Change after exposure (Individual)","Correlation w BD load (Species)", "Correlation w BD load (Indiv)")))

for ( sp in c("Anbo","Rhma","Osse","Raca","Rapi")) {
    if ( sp=="Anbo") {
        
        # Number of TOTAL individuals for this species
        n1 <- length(heatmap_inhib %>%
                         filter(temp_sp==sp) %>%
                         pull(temp_indiv))
        # Number of total COLUMNS
        n2 <- length(unique(heatmap_inhib_long$Test_category))-1
        # Number of unique individuals
        n3 <- length(unique(heatmap_inhib %>%
                                filter(temp_sp==sp) %>%
                                pull(temp_indiv)))
        
        ## Creating custom horizontal white lines
        hlines_dat <- data.frame(x = c(rep(c(2,3,4,5),times=(n1-1))# for indiv
                                       ,rep(c(1,2,3,4),times=length(seq(n3+1,n1, by = n3))))+0.5
                                 , y = c(rep(2:n1, each = 4) - 0.5 # for indiv
                                         ,rep(seq(n3+1,n1, by = n3), each=4) - 0.5 # for species
                                 )
                                 , group=c(rep(seq(1:(2*(n1-1+ length(seq(n3+1,n1, by = n3))))), each = 2) # for indivand species
                                 )
        )
        
        assign(paste0("gg_heatmap_",sp),
               heatmap_inhib_long %>%
                   filter(temp_sp==sp) %>%
                   filter(!(Test_category%in% c("Amphib_Species"))) %>%
                   mutate(Test_category = factor(Test_category, levels=c("Bacterial Genus","Amphib_Species","Change after exposure (Species)","Change after exposure (Individual)","Correlation w BD load (Species)", "Correlation w BD load (Indiv)"))) %>%
                   ggplot() +
                   geom_tile(aes(x=Test_category, y=allRows, fill=as.character(value)), show.legend=FALSE) +
                   scale_fill_manual(values=all_colors_heatmap)+
                   ylab("") +
                   theme(
                       axis.text.x = element_text(angle=90, hjust=1, vjust=0)
                       , axis.ticks = element_blank()
                       , axis.text.y=element_blank()
                       , panel.grid = element_blank()
                       , axis.title.x = element_blank()
                       , plot.margin = unit(x=c(1,-0.25,0,1), "cm")
                       # , panel.background = element_blank()
                       # , panel.grid.minor.x = element_line()
                       # , panel.grid.major=element_blank()
                   ) +
                   geom_line(
                       # data = data.frame(x = c(2,3,4,5) + 0.5, y = rep(2:n1, each = 4) - 0.5, group=rep(seq(1:(2*(n1-1))), each = 2)), ## Horizontal lines
                       data=hlines_dat,
                       aes(x = x, y = y, group = group),cex=0.5, col="white") +
                   geom_line(data = data.frame(
                       x = rep(4, each=2)-0.5, y=c(0,n1)+0.5),
                       aes(x = x, y = y, group = x), col="grey", cex=3) +
                   geom_line(data = data.frame(x = rep(1:n2, each=2)-0.5, y = c(0,n1)+0.5), ## Vertical lines
                             aes(x = x, y = y, group = x), col="white", cex=1) +
                   facet_grid(.~temp_sp)
               
        )
    }else {
        
        # Number of TOTAL individuals for this species
        n1 <- length(heatmap_inhib %>%
                         filter(temp_sp==sp) %>%
                         pull(temp_indiv))
        # Number of total COLUMNS
        n2 <- length(unique(heatmap_inhib_long$Test_category))-1
        # Number of unique individuals
        n3 <- length(unique(heatmap_inhib %>%
                                filter(temp_sp==sp) %>%
                                pull(temp_indiv)))
        
        ## Creating custom horizontal white lines
        hlines_dat <- data.frame(x = c(rep(c(1,2,3,4),times=(n1-1))# for indiv
                                       ,rep(c(0,1,2,3),times=length(seq(n3+1,n1, by = n3))))+0.5
                                 , y = c(rep(2:n1, each = 4) - 0.5 # for indiv
                                         ,rep(seq(n3+1,n1, by = n3), each=4) - 0.5 # for species
                                 )
                                 , group=c(rep(seq(1:(2*(n1-1+ length(seq(n3+1,n1, by = n3))))), each = 2) # for indivand species
                                 )
        )
        
        assign(paste0("gg_heatmap_",sp),
               heatmap_inhib_long %>%
                   filter(temp_sp==sp) %>%
                   filter(!(Test_category%in% c("Bacterial Genus","Amphib_Species"))) %>%
                   mutate(Test_category = factor(Test_category, levels=c("Bacterial Genus","Amphib_Species","Change after exposure (Species)","Change after exposure (Individual)","Correlation w BD load (Species)", "Correlation w BD load (Indiv)"))) %>%
                   ggplot() +
                   geom_tile(aes(x=Test_category, y=allRows, fill=as.character(value)), show.legend=FALSE) +
                   scale_fill_manual(values=all_colors_heatmap)+
                   ylab("") +
                   theme(
                       axis.text.x = element_text(angle=90, hjust=1, vjust=0)
                       , axis.ticks = element_blank()
                       , axis.text.y=element_blank()
                       , panel.grid = element_blank()
                       , axis.title.x = element_blank()
                       , plot.margin = unit(x=c(1,0,0,-0.25), "cm")
                       # , panel.background = element_blank()
                       # , panel.grid.minor.x = element_line()
                       # , panel.grid.major=element_blank()
                   ) +
                   geom_line(
                       # data = data.frame(x = c(2,3,4,5) + 0.5, y = rep(2:n1, each = 4) - 0.5, group=rep(seq(1:(2*(n1-1))), each = 2)), ## Horizontal lines
                       data=hlines_dat,
                       aes(x = x, y = y, group = group),cex=0.5, col="white") +
                   geom_line(data = data.frame(
                       x = rep(3, each=2)-0.5, y=c(0,n1)+0.5),
                             aes(x = x, y = y, group = x), col="grey", cex=3) +
                   geom_line(data = data.frame(
                       x = rep(1:(n2-1), each=2)-0.5, y = c(0,n1)+0.5), ## Vertical lines
                       aes(x = x, y = y, group = x), col="white", cex=1) +
                   facet_grid(.~temp_sp)
        )
    }
   
    
}


# Get legend
gg_leg_inhib <- g_legend(heatmap_inhib_long %>%
    filter(temp_sp=="Anbo") %>%
    filter((Test_category%in% c("Bacterial Genus"))) %>%
    rename(`Bacterial Genus` = value) %>%
    ggplot() +
    geom_tile(aes(x=Test_category, y=allRows, fill=`Bacterial Genus`), show.legend=TRUE) +
    scale_fill_manual(values=all_colors_heatmap) +
        guides(fill=guide_legend(reverse=TRUE)) +
        theme(
            plot.margin=unit(c(0,0,0,0),"cm")
        )
    )
gg_leg_change <- g_legend(heatmap_inhib_long %>%
                             # filter(temp_sp=="Anbo") %>%
                             filter((Test_category%in% c("Change after exposure (Species)"))) %>%
                             mutate(`Direction of change` = ifelse(value>0, "+","-")) %>%
                             filter(!is.na(`Direction of change`)) %>%
                             ggplot() +
                             geom_tile(aes(x=Test_category, y=allRows, fill=`Direction of change`), show.legend=TRUE) +
                             scale_fill_manual(values=all_colors_heatmap)+
                              theme(
                              plot.margin=unit(c(0,0,0,0),"cm")
                              )
)

pdf("FIGURES/inhib_corr_exposure_bdload.pdf", height=7, width=12)
grid.arrange(gg_heatmap_Anbo
             ,gg_heatmap_Rhma
             ,gg_heatmap_Osse
             ,gg_heatmap_Raca
             ,gg_heatmap_Rapi
             ,gg_leg_change
             , gg_leg_inhib
             , layout_matrix=rbind(c(rep(1,5),rep(c(2,3,4,5,6,6),each=4))
                                   ,c(rep(1,5),rep(c(2,3,4,5,6,6),each=4))
                                   ,c(rep(1,5),rep(c(2,3,4,5,7,7),each=4))
                                   ,c(rep(1,5),rep(c(2,3,4,5,7,7),each=4))
                                   ,c(rep(1,5),rep(c(2,3,4,5,7,7),each=4))
                                   ,c(rep(1,5),rep(c(2,3,4,5,NA,NA),each=4)))
             , bottom="Test Category"
             
)
dev.off()


#### Correlation between Richness and Inhibitory Richness ####

pdf("FIGURES/corr_rich_inhibrich.pdf", width=5, height=4)
all_p_withcon %>%
    dplyr::select(p_inhibRich, p_logRich, toadID) %>%
    separate(toadID, into=c("species","n"), remove=FALSE) %>%
    rename(Species=species) %>%
    ggplot(aes(x=p_logRich, y=p_inhibRich)) +
    geom_point(aes(col=Species), cex=3) +    
    xlab("OTU Richness (Percentile of species)") + 
    ylab("Inhibitory OTU Richness (Percentile of species)")
dev.off()


pdf("FIGURES/corr_rich__prich_inhibrich.pdf", width=10, height=4)
grid.arrange(all_p_withcon %>%
                 dplyr::select(toadID, p_inhibRich, p_logRich, p_percInhib) %>%
                 separate(toadID, into=c("species","n"), remove=FALSE) %>%
                 rename(Species=species) %>%
                 ggplot(aes(x=p_logRich, y=p_inhibRich)) +
                 geom_point(aes(col=Species), cex=2, show.legend = FALSE) +    
                 geom_smooth(method="lm", col="grey", se=FALSE, lty=2) +
                 xlab("log ASV Richness (Percentile)") + 
                 ylab("Inhibitory ASV Richness (Percentile)") +
                 theme_classic()
             , all_p_withcon %>%
                 dplyr::select(toadID, p_inhibRich, p_logRich, p_percInhib) %>%
                 separate(toadID, into=c("species","n"), remove=FALSE) %>%
                 rename(Species=species) %>%
                 ggplot(aes(x=p_percInhib, y=p_inhibRich)) +
                 geom_smooth(method="lm", col="grey", se=FALSE, lty=2)+
                 geom_point(aes(col=Species), cex=2) +    
                 xlab("Percent of ASVs Inhibitory (Percentile)") + 
                 ylab("Inhibitory ASV Richness (Percentile)") +
                 theme_classic()
             , nrow=1
)
dev.off()

mf.all.without_init_infect <- mf.rare %>%
    filter(SampleID %in% c(mf_con_without_init_infect$SampleID, mf_treat_without_init_infect))

pdf("FIGURES/corr_rich__prich_inhibrich_RAW.pdf", width=10, height=4)
grid.arrange(mf.all.without_init_infect %>%
                 dplyr::select(toadID, inhibRich, logRich, percInhib) %>%
                 separate(toadID, into=c("species","n"), remove=FALSE) %>%
                 rename(Species=species) %>%
                 ggplot(aes(x=logRich, y=inhibRich)) +
                 geom_point(aes(col=Species), cex=1, position=position_jitter(height=0.5, width=0.5), show.legend = FALSE) +    
                 ylab("ASV Richness (log observed ASVs)") + 
                 xlab("Inhibitory ASV Richness")
             , mf.all.without_init_infect %>%
                 dplyr::select(toadID, inhibRich, logRich, percInhib) %>%
                 separate(toadID, into=c("species","n"), remove=FALSE) %>%
                 rename(Species=species) %>%
                 ggplot(aes(x=percInhib, y=inhibRich)) +
                 geom_point(aes(col=Species), cex=1, position=position_jitter(height=0.5, width=0.5)) +    
                 ylab("Percent Inhibitory (as proportion of total reads)") + 
                 xlab("Inhibitory ASV Richness")
             , nrow=1
)
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

disp.normfit <- fitdistr(x=mf_con_without_init_infect$disper_bray_curtis, densfun = "normal")
disp.lognormfit <- fitdistr(x=mf_con_without_init_infect$disper_bray_curtis, densfun = "lognormal")
x.pred <- seq(min(mf_con_without_init_infect$disper_bray_curtis, na.rm = TRUE)-sd(mf_con_without_init_infect$disper_bray_curtis, na.rm = TRUE)
              , max(mf_con_without_init_infect$disper_bray_curtis, na.rm = TRUE)+sd(mf_con_without_init_infect$disper_bray_curtis, na.rm = TRUE)
              , length.out = 100)
y.expect.norm <- dnorm(x=x.pred, mean=disp.normfit$estimate[1], sd=disp.normfit$estimate[2])
y.expect.lognorm <- dlnorm(x=x.pred, meanlog=disp.lognormfit$estimate[1], sdlog = disp.lognormfit$estimate[2])

gg_disp_fit <- mf_con_without_init_infect %>%
    ggplot(aes(x=disper_bray_curtis, y=..density..)) +
    geom_histogram(bins=25) +
    geom_line(data=data.frame(x=x.pred, y=y.expect.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.pred, y=y.expect.lognorm), aes(x=x, y=y), col="blue") 


## Instability

instab.normfit <- fitdistr(x=mf_con_without_init_infect$distance_bray_curtis[!is.na(mf_con_without_init_infect$distance_bray_curtis)], densfun = "normal")
instab.beta <- fitdistr(x=mf_con_without_init_infect$distance_bray_curtis[!is.na(mf_con_without_init_infect$distance_bray_curtis)], densfun = "beta", start = list(shape1=6,shape2=6))
x.pred <- seq(min(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)-sd(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)
              , max(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)+sd(mf_con_without_init_infect$distance_bray_curtis, na.rm = TRUE)
              , length.out = 100)
y.expect.norm <- dnorm(x=x.pred, mean=instab.normfit$estimate[1], sd=instab.normfit$estimate[2])
y.expect.beta <- dbeta(x=x.pred, shape1 = instab.beta$estimate[1], shape2 = instab.beta$estimate[2])

gg_instab_fit <- mf_con_without_init_infect %>%
    ggplot(aes(x=distance_bray_curtis, y=..density..)) +
    geom_histogram(bins=25) +
    geom_line(data=data.frame(x=x.pred, y=y.expect.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.pred, y=y.expect.beta), aes(x=x, y=y), col="darkgreen") 

## InhibRich

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

## percInhib

percInhib.normfit <- fitdistr(x=mf_con_without_init_infect$percInhib[!is.na(mf_con_without_init_infect$percInhib)], densfun = "normal")
percInhib.beta <- fitdistr(x=mf_con_without_init_infect$percInhib[!is.na(mf_con_without_init_infect$percInhib)], densfun = "beta", start = list(shape1=1,shape2=6))
x.pred <- seq(0, min(max(mf_con_without_init_infect$percInhib, na.rm = TRUE)+sd(mf_con_without_init_infect$percInhib, na.rm = TRUE),1)
              , length.out = 100)
y.expect.norm <- dnorm(x=x.pred, mean=percInhib.normfit$estimate[1], sd=percInhib.normfit$estimate[2])
y.expect.beta <- dbeta(x=x.pred, shape1=percInhib.beta$estimate[1], shape2=percInhib.beta$estimate[2])


gg_percInhib_fit <- mf_con_without_init_infect %>%
    ggplot(aes(x=percInhib, y=..density..)) +
    geom_histogram(bins=25) +
    geom_line(data=data.frame(x=x.pred, y=y.expect.norm), aes(x=x, y=y), col="red") +
    geom_line(data=data.frame(x=x.pred, y=y.expect.beta), aes(x=x, y=y), col="darkgreen") 


pdf("FIGURES/comparison_fits.pdf", width=6, height=4)
grid.arrange(gg_logrich_fit, gg_disp_fit, gg_instab_fit, gg_inhibRich_fit, gg_percInhib_fit, nrow=2)
dev.off()


#### NMDS of inhibitory OTUs only #####
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

pdf("FIGURES/NMDS_inhibonly_all.pdf")
mf_combined %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=species,alpha=time, pch=BD_infected))+
    scale_shape_manual(values=c(`y`=16, `n`=21))
dev.off()
adonis(dm_bray_combo_inhib ~ species*BD_infected*time, data=mf_combined)


gg_anbo_inhib <- mf_combined %>%
    filter(species=="Anbo") %>%
    mutate(Treatment=BD_infected, Infected=PABD) %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=Treatment,size=(time), pch=Infected), show.legend = FALSE)+
    scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21))
gg_anma_inhib <- mf_combined %>%
    filter(species=="Anma") %>%
    mutate(Treatment=BD_infected, Infected=PABD) %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=Treatment,size=(time), pch=Infected), show.legend = FALSE)+
    scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21))
gg_lica_inhib <-mf_combined %>%
    filter(species=="Lica") %>%
    mutate(Treatment=BD_infected, Infected=PABD) %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=Treatment,size=(time), pch=Infected), show.legend = FALSE)+
    scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21))
gg_lipi_inhib <-mf_combined %>%
    filter(species=="Lipi") %>%
    mutate(Treatment=BD_infected, Infected=PABD) %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=Treatment,size=(time), pch=Infected), show.legend = FALSE)+
    scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21))
gg_osse_inhib <-mf_combined %>%
    filter(species=="Osse") %>%
    mutate(Treatment=BD_infected, Infected=PABD) %>%
    ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
    geom_point(aes(col=Treatment,size=(time), pch=Infected), show.legend = FALSE)+
    scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21))
gg_nmdsinhib_leg <- g_legend(mf_combined %>%
             filter(species=="Osse") %>%
             mutate(Treatment=BD_infected, Infected=PABD) %>%
             ggplot(aes(x=NMDS1_inhib, y=NMDS2_inhib)) +
             geom_point(aes(col=Treatment,size=(time), pch=Infected), show.legend = TRUE)+
             scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21)))

pdf("FIGURES/nmds_allsp_inhib.pdf", width=12, height=10)
grid.arrange(gg_anbo_inhib, gg_anma_inhib, gg_lica_inhib, gg_lipi_inhib, gg_osse_inhib, gg_nmdsinhib_leg, nrow=2)

dev.off()


