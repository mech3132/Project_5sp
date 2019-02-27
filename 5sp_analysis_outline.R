#' ---
#' title: Analysis outline 
#' author: Melissa Chen
#' output: github_document
#' ---
#' 
#' 
#' EDA:
#' - How should you adjust BD estimates (done)
#' - Check for duplicates, missing data, etc (done)
#' - ** Check controls and contaminants; plot **
#' - 
#' 
#' 
#' How many otus are shared over time on the same individual? \
#' How does this compare to how many otus are shared between individuals at the same time point? \
#' Can we use structural equation modelling to determine whether BD influences microbiome or microbiome influences BD? \
#' 
#' Part I: Microbiome state and effect on the microbiome\
#' \
#' (1a) Does overall diversity of microbiome influence BD infection rate?\
#' - OTU richness vs infection rate\
#' - Shannon diversity vs infection rate\
#' - Chao1 vs infection rate\
#' - PD vs infection rate\
#' --> Use distribution of diversity prior to BD exposure as a predictor for infection rate\
#' --> P/A bd infection of entire "after" exposure\
#' 
#'\
#' (1b) Does overall diversity of microbiome influence BD infection intensity?\
#' - OTU richness vs infection intensity\
#' - Shannon diversity vs infection intensity\
#' - Chao1 vs infection intensity\
#' - PD vs infection intensity\
#' --> use distribution of diversity prior to BD exposure as a predictor for infection intensity\
#' --> max bd infection of entire "after" exposure\
#' 
#' \
#' (2a) Does instability of microbiome influence BD infection rate?\
#' - Bray-curtis distance travelled vs infection rate\
#' - Weighted unifrav distance travelled vs infeciton rate\
#' - Manual calculation of seq turnovervs infection rate\
#' --> Use a GLM (binomial) where each distance travelled "X" time points away from BD exposure\
#' --> P/A bd infection of entire "after" exposure\
#' 
#' \
#' (2b) Does instability of microbiome influence BD infection intensity?\
#' - Bray-curtis distance travelled vs infection intensity\
#' - Weighted unifrav distance travelled vs infeciton intensity\
#' - Manual calculation of seq turnovervs infection intensity\
#' --> Use a GLM (binomial) where each distance travelled "X" time points away from BD exposure\
#' --> max bd infection of entire "after" exposure\
#' 
#' \
#' (3a) Does composition of microbiome influence BD infection intensity?\
#' - Inhibitory OTU richness vs infection rate\
#' - Percent inhibitory otus vs infection rate\
#' - "Other" OTUs vs infection rate\
#' --> For inhibitory OTUs, calculate distribution of richness/percent inhibitory prior to exposure as predictor for infection rate
#' --> for "Other" OTUs, try very elaborate linear model of entire microbiome? Use AIC to find best predictor for BD infection?
#' 
#' 
#' Part II: Affect of BD infection on microbiome state
#' 
#' (1a) Does BD infection state affect microbiome diversity?
#' - OTU richness vs BD infection
#' - Chao1 richness vs BD infection
#' - Shannon richness vs BD infection
#' - PD vs BD infection
#' 
#' (1b) Does BD infection intensity affect microbiome diversity?
#' - OTU richness vs BD infection
#' - Chao1 richness vs BD infection
#' - Shannon richness vs BD infection
#' - PD vs BD infection
#' 