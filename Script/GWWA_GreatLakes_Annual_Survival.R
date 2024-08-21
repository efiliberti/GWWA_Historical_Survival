## This script runs a hierarchical multi-population apparent annual survival CJS 
## model for the Great Lakes population of Golden-winged Warbler. Eleven
## CMR historical datasets ranging from 1980-2022 were collected and compiled. 
## We selected vague prior distributions for all parameters to reflect a lack
## of prior knowledge on range-wide annual survival for this declining species. 
## Year is treated as a random effect to account for unknown and variable search
## effort throughout individual study periods. This model aims to produce
## informed regional-, site-, and sex-specific apparent annual survival 
## estimates to identify particular spatial or demographic populations that 
## might be driving declines of these populations. Sex effect was only assessed
## for sites that had at least 1% female representation in the dataset.

## Load packages
library(IPMbook)
library(jagsUI)
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)
library(sjPlot)
library(MCMCvis)

set.seed(1313)


################################################################################
###                              MODEL FUNCTIONS                             ###
################################################################################

## Establish functions for first marking, latent state, and initial values.
get.first <- function(x) min(which(x!=0))

known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

# Create function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}


################################################################################
###                           DATA FORMATTING                                ###
################################################################################

## Read in data file
comboCH <- read.csv("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Script/All_GL_v3.csv")

## Combine sites that are within 40 km spatially
CCR.BCR <- c("CCR", "BCR")

## Separate capture histories into specific sites
CCR<-subset(comboCH, Site %in% CCR.BCR) ## Costa Rica
MMA<-subset(comboCH, Site == "MMA") ## Southeastern Manitoba
RMI<-subset(comboCH, Site == "RMI") ## Northern Michigan
RWI<-subset(comboCH, Site == "RWI") ## Wisconsin
WMI<-subset(comboCH, Site == "WMI") ## Central Michigan
VON<-subset(comboCH, Site == "VON") ## Eastern Ontario
BGU<-subset(comboCH, Site == "BGU") ## Guatemala
BHO<-subset(comboCH, Site == "BHO") ## Honduras
CNI<-subset(comboCH, Site == "CNI") ## Nicaragua
VMA<-subset(comboCH, Site == "VMA") ## Southwestern Manitoba
VONW<-subset(comboCH, Site == "VONW") ## Western Ontario

## Pull out capture histories for each site
CCR_CH <- CCR[,27:37]
MMA_CH <- MMA[,29:35]
RMI_CH <- RMI[,34:38]
RWI_CH <- RWI[,28:41]
WMI_CH <- WMI[,1:5]
VON_CH <- VON[,22:27]
BGU_CH <- BGU[,35:37]
BHO_CH <- BHO[,32:37]
CNI_CH <- CNI[,23:43]
VMA_CH <- VMA[,29:31]
VONW_CH <- VONW[,29:31]

## Get first occassion for each site
f.CCR <- apply(CCR_CH, 1, get.first)
f.MMA <- apply(MMA_CH, 1, get.first)
f.RMI <- apply(RMI_CH, 1, get.first)
f.RWI <- apply(RWI_CH, 1, get.first)
f.WMI <- apply(WMI_CH, 1, get.first)
f.VON <- apply(VON_CH, 1, get.first)
f.BGU <- apply(BGU_CH, 1, get.first)
f.BHO <- apply(BHO_CH, 1, get.first)
f.CNI <- apply(CNI_CH, 1, get.first)
f.VMA <- apply(VMA_CH, 1, get.first)
f.VONW <- apply(VONW_CH, 1, get.first)

## Pull out sex covariate for each site
CCR.Sex <- CCR[,47]
MMA.Sex <- MMA[,47]
RMI.Sex <- RMI[,47]
RWI.Sex <- RWI[,47]
WMI.Sex <- WMI[,47]
VON.Sex <- VON[,47]
BGU.Sex <- BGU[,47]
BHO.Sex <- BHO[,47]
CNI.Sex <- CNI[,47]
VMA.Sex <- VMA[,47]
VONW.Sex <- VONW[,47]

## Convert all sex covariates to numeric
CCR.Sex <- as.numeric(CCR.Sex)
MMA.Sex <- as.numeric(MMA.Sex)
RMI.Sex <- as.numeric(RMI.Sex)
RWI.Sex <- as.numeric(RWI.Sex)
WMI.Sex <- as.numeric(WMI.Sex)
VON.Sex <- as.numeric(VON.Sex)
BGU.Sex <- as.numeric(BGU.Sex)
BHO.Sex <- as.numeric(BHO.Sex)
CNI.Sex <- as.numeric(CNI.Sex)
VMA.Sex <- as.numeric(VMA.Sex)
VONW.Sex <- as.numeric(VONW.Sex)

## Create vector of number of individuals for each site
nind.CCR<-nrow(CCR_CH)
nind.MMA<-nrow(MMA_CH)
nind.RMI<-nrow(RMI_CH)
nind.RWI<-nrow(RWI_CH)
nind.WMI<-nrow(WMI_CH)
nind.VON<-nrow(VON_CH)
nind.BGU<-nrow(BGU_CH)
nind.BHO<-nrow(BHO_CH)
nind.CNI<-nrow(CNI_CH)
nind.VMA<-nrow(VMA_CH)
nind.VONW<-nrow(VONW_CH)


## Create vector of number of years for each site
nyears.CCR <-length(2006:2016)
nyears.MMA <-length(2008:2014)
nyears.RMI <-length(2013:2017)
nyears.RWI <-length(2007:2020)
nyears.WMI <-length(1980:1984)
nyears.VON <-length(2001:2006)
nyears.BGU <-length(2014:2016)
nyears.BHO <-length(2011:2016)
nyears.CNI <-length(2002:2022)
nyears.VMA <-length(2008:2010)
nyears.VONW <-length(2008:2010)

## Find the sum of both sexes 
gl_sex_count <- rbind(table(CCR.Sex), table(MMA.Sex), table(RMI.Sex),
                      table(RWI.Sex), table(WMI.Sex), table(VON.Sex),
                      table(BGU.Sex), table(BHO.Sex), table(CNI.Sex),
                      table(VMA.Sex), table(VONW.Sex))
colnames(gl_sex_count) <- c("males", "females")

gl_sex_labels <- c("CCR", "MMA", "RMI", "RWI", "WMI", "VON", "BGU", "BHO",
                   "CNI", "VMA", "VONW")
gl_sex_count <- cbind(gl_sex_count, gl_sex_labels)

gl_sex_count <- as.data.frame(gl_sex_count)
gl_sex_count$males <- as.numeric(gl_sex_count$males)
gl_sex_count$females <- as.numeric(gl_sex_count$females)

################################################################################
###                             BUILD THE MODEL                              ###
################################################################################

sink("GLMultiSiteNo2year.jags")
cat("

model {

  ##############################################
  # Define model likelihood for mean survival of a site
  for (s in 1:nsite){
    # s = site effect
    mu.GL[s] ~ dnorm(overall.survival.GL, tau.site.GL)
   
  }
 
  # Define prior for overall mean survival

  overall.survival.GL ~ dnorm(0,1.0E-3)

  #overall.survival.GL ~ dbeta(1,1)
  logit(mean.overall.survival.GL) <- overall.survival.GL
  # Define prior for inter-site variance (site random effect)
 
  tau.site.GL <- 1 / (sd.site.GL * sd.site.GL)
  sd.site.GL ~ dunif(0,5)
 
  # Monitor variance
 
  var.site.GL <- 1/tau.site.GL
 

############################## START CCR #######################################

# Priors and constraints
 # for(i in 1:nind.CCR){
   for (t in 1:(n.occasions.CCR-1)){ # time dependent survival
    logit(phi.CCR[1,t]) <- logitphi.CCR[1,t]
    logit(phi.CCR[2,t]) <- logitphi.CCR[2,t]
    logitphi.CCR[1,t] ~ dnorm(mu.GL[1],tau.eps.CCR)
    logitphi.CCR[2,t] ~ dnorm(mu.GL[1],tau.eps.CCR)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.CCR[t]-CCR.beta.sex*CCR.sex[i])) # + CCR.beta.sex*CCR.sex[i]
 
    p.CCR[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.CCR-1)){
  delta.phi.CCR[t] <- phi.CCR[1,t]-phi.CCR[2,t]
  #logitphi.CCR[t] ~ dnorm(mu.App[2],tau.eps.CCR)
  mean.m.CCR.phi[t] <- mean(phi.CCR[1,t])
  mean.f.CCR.phi[t] <- mean(phi.CCR[2,t])
  }
 
overall.f.CCR.phi <- mean(mean.f.CCR.phi)
overall.m.CCR.phi <- mean(mean.m.CCR.phi)
overall.CCR.phi <- mean(phi.CCR)
overall.delta.phi.CCR <- mean(delta.phi.CCR)
overall.CCR.p <- mean(p.CCR)

  # Define prior for variance of temporal random effect
  tau.eps.CCR <- 1 / (sd.eps.CCR*sd.eps.CCR)
  sd.eps.CCR ~ dunif(0,5)
  # Monitor variance
  var.eps.CCR <- 1 / tau.eps.CCR

# Likelihood
for (i in 1:nind.CCR){
   # Define latent state at first capture
   z.CCR[i,f.CCR[i]] ~ dbern(1)
   for (t in (f.CCR[i]+1):n.occasions.CCR){
      # State process
      z.CCR[i,t] ~  dbern(phi.CCR[CCR.sex[i], t-1] * z.CCR[i,t-1])
      #z.CCR[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.CCR[i, t-1] * z.CCR[i,t-1])
      # Observation process
      y.CCR[i,t] ~ dbern(p.CCR[t-1] * z.CCR[i,t])
      #y.CCR[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.CCR[i,t-1] * z.CCR[i,t])
      } #t
   } #i

################################ END CCR #######################################

############################## START MMA #######################################

# Priors and constraints
 # for(i in 1:nind.MMA){
   for (t in 1:(n.occasions.MMA-1)){ # time dependent survival
    logit(phi.MMA[1,t]) <- logitphi.MMA[1,t]
    logit(phi.MMA[2,t]) <- logitphi.MMA[2,t]
    logitphi.MMA[1,t] ~ dnorm(mu.GL[2],tau.eps.MMA)
    logitphi.MMA[2,t] ~ dnorm(mu.GL[2],tau.eps.MMA)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.MMA[t]-MMA.beta.sex*MMA.sex[i])) # + MMA.beta.sex*MMA.sex[i]
 
    p.MMA[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.MMA-1)){
  delta.phi.MMA[t] <- phi.MMA[1,t]-phi.MMA[2,t]
  #logitphi.MMA[t] ~ dnorm(mu.App[3],tau.eps.MMA)
  mean.m.MMA.phi[t] <- mean(phi.MMA[1,t])
  mean.f.MMA.phi[t] <- mean(phi.MMA[2,t])
  }
 
overall.f.MMA.phi <- mean(mean.f.MMA.phi)
overall.m.MMA.phi <- mean(mean.m.MMA.phi)
overall.MMA.phi <- mean(phi.MMA)
overall.delta.phi.MMA <- mean(delta.phi.MMA)
overall.MMA.p <- mean(p.MMA)
 
  # Define prior for variance of temporal random effect
  tau.eps.MMA <- 1 / (sd.eps.MMA*sd.eps.MMA)
  sd.eps.MMA ~ dunif(0,5)
  # Monitor variance
  var.eps.MMA <- 1 / tau.eps.MMA

# Likelihood
for (i in 1:nind.MMA){
   # Define latent state at first capture
   z.MMA[i,f.MMA[i]] ~ dbern(1)
   for (t in (f.MMA[i]+1):n.occasions.MMA){
      # State process
      z.MMA[i,t] ~  dbern(phi.MMA[MMA.sex[i], t-1] * z.MMA[i,t-1])
      #z.MMA[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.MMA[i, t-1] * z.MMA[i,t-1])
      # Observation process
      y.MMA[i,t] ~ dbern(p.MMA[t-1] * z.MMA[i,t])
      #y.MMA[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.MMA[i,t-1] * z.MMA[i,t])
      } #t
   } #i

################################ END MMA #######################################

############################## START RMI #######################################

# Priors and constraints
 # for(i in 1:nind.RMI){
   for (t in 1:(n.occasions.RMI-1)){ # time dependent survival
    logit(phi.RMI[1,t]) <- logitphi.RMI[1,t]
    logit(phi.RMI[2,t]) <- logitphi.RMI[2,t]
    logitphi.RMI[1,t] ~ dnorm(mu.GL[3],tau.eps.RMI)
    logitphi.RMI[2,t] ~ dnorm(mu.GL[3],tau.eps.RMI)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.RMI[t]-RMI.beta.sex*RMI.sex[i])) # + RMI.beta.sex*RMI.sex[i]
 
    p.RMI[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.RMI-1)){
  delta.phi.RMI[t] <- phi.RMI[1,t]-phi.RMI[2,t]
  #logitphi.RMI[t] ~ dnorm(mu.App[4],tau.eps.RMI)
  mean.m.RMI.phi[t] <- mean(phi.RMI[1,t])
  mean.f.RMI.phi[t] <- mean(phi.RMI[2,t])
  }
 
overall.f.RMI.phi <- mean(mean.f.RMI.phi)
overall.m.RMI.phi <- mean(mean.m.RMI.phi)
overall.RMI.phi <- mean(phi.RMI)
overall.delta.phi.RMI <- mean(delta.phi.RMI)
overall.RMI.p <- mean(p.RMI)
 
  # Define prior for variance of temporal random effect
  tau.eps.RMI <- 1 / (sd.eps.RMI*sd.eps.RMI)
  sd.eps.RMI ~ dunif(0,5)
  # Monitor variance
  var.eps.RMI <- 1 / tau.eps.RMI

# Likelihood
for (i in 1:nind.RMI){
   # Define latent state at first capture
   z.RMI[i,f.RMI[i]] ~ dbern(1)
   for (t in (f.RMI[i]+1):n.occasions.RMI){
      # State process
      z.RMI[i,t] ~  dbern(phi.RMI[RMI.sex[i], t-1] * z.RMI[i,t-1])
      #z.RMI[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.RMI[i, t-1] * z.RMI[i,t-1])
      # Observation process
      y.RMI[i,t] ~ dbern(p.RMI[t-1] * z.RMI[i,t])
      #y.RMI[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.RMI[i,t-1] * z.RMI[i,t])
      } #t
   } #i

################################ END RMI #######################################

############################## START RWI #######################################

# Priors and constraints
 # for(i in 1:nind.RWI){
   for (t in 1:(n.occasions.RWI-1)){ # time dependent survival
    logit(phi.RWI[1,t]) <- logitphi.RWI[1,t]
    logit(phi.RWI[2,t]) <- logitphi.RWI[2,t]
    logitphi.RWI[1,t] ~ dnorm(mu.GL[4],tau.eps.RWI)
    logitphi.RWI[2,t] ~ dnorm(mu.GL[4],tau.eps.RWI)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.RWI[t]-RWI.beta.sex*RWI.sex[i])) # + RWI.beta.sex*RWI.sex[i]
 
    p.RWI[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.RWI-1)){
  delta.phi.RWI[t] <- phi.RWI[1,t]-phi.RWI[2,t]
  #logitphi.RWI[t] ~ dnorm(mu.App[5],tau.eps.RWI)
  mean.m.RWI.phi[t] <- mean(phi.RWI[1,t])
  mean.f.RWI.phi[t] <- mean(phi.RWI[2,t])
  }
 
overall.f.RWI.phi <- mean(mean.f.RWI.phi)
overall.m.RWI.phi <- mean(mean.m.RWI.phi)
overall.RWI.phi <- mean(phi.RWI)
overall.delta.phi.RWI <- mean(delta.phi.RWI)
overall.RWI.p <- mean(p.RWI)
 
  # Define prior for variance of temporal random effect
  tau.eps.RWI <- 1 / (sd.eps.RWI*sd.eps.RWI)
  sd.eps.RWI ~ dunif(0,5)
  # Monitor variance
  var.eps.RWI <- 1 / tau.eps.RWI

# Likelihood
for (i in 1:nind.RWI){
   # Define latent state at first capture
   z.RWI[i,f.RWI[i]] ~ dbern(1)
   for (t in (f.RWI[i]+1):n.occasions.RWI){
      # State process
      z.RWI[i,t] ~  dbern(phi.RWI[RWI.sex[i], t-1] * z.RWI[i,t-1])
      #z.RWI[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.RWI[i, t-1] * z.RWI[i,t-1])
      # Observation process
      y.RWI[i,t] ~ dbern(p.RWI[t-1] * z.RWI[i,t])
      #y.RWI[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.RWI[i,t-1] * z.RWI[i,t])
      } #t
   } #i

################################ END RWI #######################################

############################## START WMI #######################################

# Priors and constraints
 # for(i in 1:nind.WMI){
   for (t in 1:(n.occasions.WMI-1)){ # time dependent survival
    logit(phi.WMI[1,t]) <- logitphi.WMI[1,t]
    logit(phi.WMI[2,t]) <- logitphi.WMI[2,t]
    logitphi.WMI[1,t] ~ dnorm(mu.GL[5],tau.eps.WMI)
    logitphi.WMI[2,t] ~ dnorm(mu.GL[5],tau.eps.WMI)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.WMI[t]-WMI.beta.sex*WMI.sex[i])) # + WMI.beta.sex*WMI.sex[i]
 
    p.WMI[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.WMI-1)){
  delta.phi.WMI[t] <- phi.WMI[1,t]-phi.WMI[2,t]
  #logitphi.WMI[t] ~ dnorm(mu.App[6],tau.eps.WMI)
  mean.m.WMI.phi[t] <- mean(phi.WMI[1,t])
  mean.f.WMI.phi[t] <- mean(phi.WMI[2,t])
  }
 
overall.f.WMI.phi <- mean(mean.f.WMI.phi)
overall.m.WMI.phi <- mean(mean.m.WMI.phi)
overall.WMI.phi <- mean(phi.WMI)
overall.delta.phi.WMI <- mean(delta.phi.WMI)
overall.WMI.p <- mean(p.WMI)

  # Define prior for variance of temporal random effect
  tau.eps.WMI <- 1 / (sd.eps.WMI*sd.eps.WMI)
  sd.eps.WMI ~ dunif(0,5)
  # Monitor variance
  var.eps.WMI <- 1 / tau.eps.WMI

# Likelihood
for (i in 1:nind.WMI){
   # Define latent state at first capture
   z.WMI[i,f.WMI[i]] ~ dbern(1)
   for (t in (f.WMI[i]+1):n.occasions.WMI){
      # State process
      z.WMI[i,t] ~  dbern(phi.WMI[WMI.sex[i], t-1] * z.WMI[i,t-1])
      #z.WMI[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.WMI[i, t-1] * z.WMI[i,t-1])
      # Observation process
      y.WMI[i,t] ~ dbern(p.WMI[t-1] * z.WMI[i,t])
      #y.WMI[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.WMI[i,t-1] * z.WMI[i,t])
      } #t
   } #i

################################ END WMI #######################################

############################## START VON #######################################

# Priors and constraints
 # for(i in 1:nind.VON){
   for (t in 1:(n.occasions.VON-1)){ # time dependent survival
    logit(phi.VON[1,t]) <- logitphi.VON[1,t]
    logit(phi.VON[2,t]) <- logitphi.VON[2,t]
    logitphi.VON[1,t] ~ dnorm(mu.GL[6],tau.eps.VON)
    logitphi.VON[2,t] ~ dnorm(mu.GL[6],tau.eps.VON)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.VON[t]-VON.beta.sex*VON.sex[i])) # + VON.beta.sex*VON.sex[i]
 
    p.VON[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.VON-1)){
  delta.phi.VON[t] <- phi.VON[1,t]-phi.VON[2,t]
  #logitphi.VON[t] ~ dnorm(mu.App[7],tau.eps.VON)
  mean.m.VON.phi[t] <- mean(phi.VON[1,t])
  mean.f.VON.phi[t] <- mean(phi.VON[2,t])
  }
 
overall.f.VON.phi <- mean(mean.f.VON.phi)
overall.m.VON.phi <- mean(mean.m.VON.phi)
overall.VON.phi <- mean(phi.VON)
overall.delta.phi.VON <- mean(delta.phi.VON)
overall.VON.p <- mean(p.VON)
 
  # Define prior for variance of temporal random effect
  tau.eps.VON <- 1 / (sd.eps.VON*sd.eps.VON)
  sd.eps.VON ~ dunif(0,5)
  # Monitor variance
  var.eps.VON <- 1 / tau.eps.VON

# Likelihood
for (i in 1:nind.VON){
   # Define latent state at first capture
   z.VON[i,f.VON[i]] ~ dbern(1)
   for (t in (f.VON[i]+1):n.occasions.VON){
      # State process
      z.VON[i,t] ~  dbern(phi.VON[VON.sex[i], t-1] * z.VON[i,t-1])
      #z.VON[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.VON[i, t-1] * z.VON[i,t-1])
      # Observation process
      y.VON[i,t] ~ dbern(p.VON[t-1] * z.VON[i,t])
      #y.VON[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.VON[i,t-1] * z.VON[i,t])
      } #t
   } #i

################################ END VON #######################################

############################## START BGU #######################################

# Priors and constraints
 # for(i in 1:nind.BGU){
   for (t in 1:(n.occasions.BGU-1)){ # time dependent survival
    logit(phi.BGU[1,t]) <- logitphi.BGU[1,t]
    logit(phi.BGU[2,t]) <- logitphi.BGU[2,t]
    logitphi.BGU[1,t] ~ dnorm(mu.GL[7],tau.eps.BGU)
    logitphi.BGU[2,t] ~ dnorm(mu.GL[7],tau.eps.BGU)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.BGU[t]-BGU.beta.sex*BGU.sex[i])) # + BGU.beta.sex*BGU.sex[i]
 
    p.BGU[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.BGU-1)){
  delta.phi.BGU[t] <- phi.BGU[1,t]-phi.BGU[2,t]
  #logitphi.BGU[t] ~ dnorm(mu.App[9],tau.eps.BGU)
  mean.m.BGU.phi[t] <- mean(phi.BGU[1,t])
  mean.f.BGU.phi[t] <- mean(phi.BGU[2,t])
  }
 
overall.f.BGU.phi <- mean(mean.f.BGU.phi)
overall.m.BGU.phi <- mean(mean.m.BGU.phi)
overall.BGU.phi <- mean(phi.BGU)
overall.delta.phi.BGU <- mean(delta.phi.BGU)
overall.BGU.p <- mean(p.BGU)
 
  # Define prior for variance of temporal random effect
  tau.eps.BGU <- 1 / (sd.eps.BGU*sd.eps.BGU)
  sd.eps.BGU ~ dunif(0,5)
  # Monitor variance
  var.eps.BGU <- 1 / tau.eps.BGU

# Likelihood
for (i in 1:nind.BGU){
   # Define latent state at first capture
   z.BGU[i,f.BGU[i]] ~ dbern(1)
   for (t in (f.BGU[i]+1):n.occasions.BGU){
      # State process
      z.BGU[i,t] ~  dbern(phi.BGU[BGU.sex[i], t-1] * z.BGU[i,t-1])
      #z.BGU[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.BGU[i, t-1] * z.BGU[i,t-1])
      # Observation process
      y.BGU[i,t] ~ dbern(p.BGU[t-1] * z.BGU[i,t])
      #y.BGU[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.BGU[i,t-1] * z.BGU[i,t])
      } #t
   } #i

################################ END BGU #######################################

############################## START BHO #######################################

# Priors and constraints
 # for(i in 1:nind.BHO){
   for (t in 1:(n.occasions.BHO-1)){ # time dependent survival
    logit(phi.BHO[1,t]) <- logitphi.BHO[1,t]
    logit(phi.BHO[2,t]) <- logitphi.BHO[2,t]
    logitphi.BHO[1,t] ~ dnorm(mu.GL[8],tau.eps.BHO)
    logitphi.BHO[2,t] ~ dnorm(mu.GL[8],tau.eps.BHO)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.BHO[t]-BHO.beta.sex*BHO.sex[i])) # + BHO.beta.sex*BHO.sex[i]
 
    p.BHO[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.BHO-1)){
  delta.phi.BHO[t] <- phi.BHO[1,t]-phi.BHO[2,t]
  #logitphi.BHO[t] ~ dnorm(mu.App[10],tau.eps.BHO)
  mean.m.BHO.phi[t] <- mean(phi.BHO[1,t])
  mean.f.BHO.phi[t] <- mean(phi.BHO[2,t])
  }
 
overall.f.BHO.phi <- mean(mean.f.BHO.phi)
overall.m.BHO.phi <- mean(mean.m.BHO.phi)
overall.BHO.phi <- mean(phi.BHO)
overall.delta.phi.BHO <- mean(delta.phi.BHO)
overall.BHO.p <- mean(p.BHO)
 
  # Define prior for variance of temporal random effect
  tau.eps.BHO <- 1 / (sd.eps.BHO*sd.eps.BHO)
  sd.eps.BHO ~ dunif(0,5)
  # Monitor variance
  var.eps.BHO <- 1 / tau.eps.BHO

# Likelihood
for (i in 1:nind.BHO){
   # Define latent state at first capture
   z.BHO[i,f.BHO[i]] ~ dbern(1)
   for (t in (f.BHO[i]+1):n.occasions.BHO){
      # State process
      z.BHO[i,t] ~  dbern(phi.BHO[BHO.sex[i], t-1] * z.BHO[i,t-1])
      #z.BHO[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.BHO[i, t-1] * z.BHO[i,t-1])
      # Observation process
      y.BHO[i,t] ~ dbern(p.BHO[t-1] * z.BHO[i,t])
      #y.BHO[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.BHO[i,t-1] * z.BHO[i,t])
      } #t
   } #i

################################ END BHO #######################################

############################## START CNI #######################################

# Priors and constraints
 # for(i in 1:nind.CNI){
   for (t in 1:(n.occasions.CNI-1)){ # time dependent survival
    logit(phi.CNI[1,t]) <- logitphi.CNI[1,t]
    logit(phi.CNI[2,t]) <- logitphi.CNI[2,t]
    logitphi.CNI[1,t] ~ dnorm(mu.GL[9],tau.eps.CNI)
    logitphi.CNI[2,t] ~ dnorm(mu.GL[9],tau.eps.CNI)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.CNI[t]-CNI.beta.sex*CNI.sex[i])) # + CNI.beta.sex*CNI.sex[i]
 
    p.CNI[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.CNI-1)){
  delta.phi.CNI[t] <- phi.CNI[1,t]-phi.CNI[2,t]
  #logitphi.CNI[t] ~ dnorm(mu.App[12],tau.eps.CNI)
  mean.m.CNI.phi[t] <- mean(phi.CNI[1,t])
  mean.f.CNI.phi[t] <- mean(phi.CNI[2,t])
  }
 
overall.f.CNI.phi <- mean(mean.f.CNI.phi)
overall.m.CNI.phi <- mean(mean.m.CNI.phi)
overall.CNI.phi <- mean(phi.CNI)
overall.delta.phi.CNI <- mean(delta.phi.CNI)
overall.CNI.p <- mean(p.CNI)
 
  # Define prior for variance of temporal random effect
  tau.eps.CNI <- 1 / (sd.eps.CNI*sd.eps.CNI)
  sd.eps.CNI ~ dunif(0,5)
  # Monitor variance
  var.eps.CNI <- 1 / tau.eps.CNI

# Likelihood
for (i in 1:nind.CNI){
   # Define latent state at first capture
   z.CNI[i,f.CNI[i]] ~ dbern(1)
   for (t in (f.CNI[i]+1):n.occasions.CNI){
      # State process
      z.CNI[i,t] ~  dbern(phi.CNI[CNI.sex[i], t-1] * z.CNI[i,t-1])
      #z.CNI[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.CNI[i, t-1] * z.CNI[i,t-1])
      # Observation process
      y.CNI[i,t] ~ dbern(p.CNI[t-1] * z.CNI[i,t])
      #y.CNI[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.CNI[i,t-1] * z.CNI[i,t])
      } #t
   } #i

################################ END CNI #######################################

############################## START VMA #######################################

# Priors and constraints
 # for(i in 1:nind.VMA){
   for (t in 1:(n.occasions.VMA-1)){ # time dependent survival
    logit(phi.VMA[1,t]) <- logitphi.VMA[1,t]
    logit(phi.VMA[2,t]) <- logitphi.VMA[2,t]
    logitphi.VMA[1,t] ~ dnorm(mu.GL[10],tau.eps.VMA)
    logitphi.VMA[2,t] ~ dnorm(mu.GL[10],tau.eps.VMA)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.VMA[t]-VMA.beta.sex*VMA.sex[i])) # + VMA.beta.sex*VMA.sex[i]
 
    p.VMA[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.VMA-1)){
  delta.phi.VMA[t] <- phi.VMA[1,t]-phi.VMA[2,t]
  #logitphi.VMA[t] ~ dnorm(mu.App[13],tau.eps.VMA)
  mean.m.VMA.phi[t] <- mean(phi.VMA[1,t])
  mean.f.VMA.phi[t] <- mean(phi.VMA[2,t])
  }
 
overall.f.VMA.phi <- mean(mean.f.VMA.phi)
overall.m.VMA.phi <- mean(mean.m.VMA.phi)
overall.VMA.phi <- mean(phi.VMA)
overall.delta.phi.VMA <- mean(delta.phi.VMA)
overall.VMA.p <- mean(p.VMA)
 
  # Define prior for variance of temporal random effect
  tau.eps.VMA <- 1 / (sd.eps.VMA*sd.eps.VMA)
  sd.eps.VMA ~ dunif(0,5)
  # Monitor variance
  var.eps.VMA <- 1 / tau.eps.VMA

# Likelihood
for (i in 1:nind.VMA){
   # Define latent state at first capture
   z.VMA[i,f.VMA[i]] ~ dbern(1)
   for (t in (f.VMA[i]+1):n.occasions.VMA){
      # State process
      z.VMA[i,t] ~  dbern(phi.VMA[VMA.sex[i], t-1] * z.VMA[i,t-1])
      #z.VMA[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.VMA[i, t-1] * z.VMA[i,t-1])
      # Observation process
      y.VMA[i,t] ~ dbern(p.VMA[t-1] * z.VMA[i,t])
      #y.VMA[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.VMA[i,t-1] * z.VMA[i,t])
      } #t
   } #i

################################ END VMA #######################################

############################## START VONW ######################################

# Priors and constraints

   for (t in 1:(n.occasions.VONW-1)){ # time dependent survival
    logit(phi.VONW[t]) <- logitphi.VONW[t]
    logitphi.VONW[t] ~ dnorm(mu.GL[11],tau.eps.VONW)
 
    p.VONW[t] ~ dbeta(1,1)
    #epsilon[t] ~ dnorm(0, tau)
      } #t

overall.m.VONW.phi <- mean(phi.VONW)
overall.VONW.phi <- mean(phi.VONW)
overall.VONW.p <- mean(p.VONW)

  # Define prior for variance of temporal random effect
  tau.eps.VONW <- 1 / (sd.eps.VONW*sd.eps.VONW)
  sd.eps.VONW ~ dunif(0,5)
  # Monitor variance
  var.eps.VONW <- 1 / tau.eps.VONW

# Likelihood
for (i in 1:nind.VONW){
   # Define latent state at first capture
   z.VONW[i,f.VONW[i]] ~ dbern(1)
   for (t in (f.VONW[i]+1):n.occasions.VONW){
      # State process
      z.VONW[i,t] ~ dbern(phi.VONW[t-1] * z.VONW[i,t-1])
      # Observation process
      y.VONW[i,t] ~ dbern(p.VONW[t-1] * z.VONW[i,t])
      } #t
   } #i

################################ END VONW ######################################

## Get overall male and female survival estimates
overall.female.survival.gl <- mean(c(overall.f.CCR.phi, overall.f.MMA.phi, overall.f.RMI.phi, overall.f.RWI.phi,
  overall.f.WMI.phi, overall.f.VON.phi, overall.f.BGU.phi, overall.f.BHO.phi, overall.f.CNI.phi,
  overall.f.VMA.phi))

overall.male.survival.gl <- mean(c(overall.m.CCR.phi, overall.m.MMA.phi, overall.m.RMI.phi, overall.m.RWI.phi,
  overall.m.WMI.phi, overall.m.VON.phi, overall.m.BGU.phi, overall.m.BHO.phi, overall.m.CNI.phi,
   overall.m.VMA.phi, overall.m.VONW.phi))
}"
,
fill = T)
sink()


################################################################################
###                               RUN THE MODEL                              ###
################################################################################

## Organize data for all sub-models
jags.data <- list(y.CCR = CCR_CH, n.occasions.CCR = ncol(CCR_CH), nind.CCR = nind.CCR, f.CCR = f.CCR, z.CCR = known.state.cjs(CCR_CH), CCR.sex = CCR.Sex,
                  y.MMA = MMA_CH, n.occasions.MMA = ncol(MMA_CH), nind.MMA = nind.MMA, f.MMA = f.MMA, z.MMA = known.state.cjs(MMA_CH), MMA.sex = MMA.Sex,
                  y.RMI = RMI_CH, n.occasions.RMI = ncol(RMI_CH), nind.RMI = nind.RMI, f.RMI = f.RMI, z.RMI = known.state.cjs(RMI_CH), RMI.sex = RMI.Sex,
                  y.RWI = RWI_CH, n.occasions.RWI = ncol(RWI_CH), nind.RWI = nind.RWI, f.RWI = f.RWI, z.RWI = known.state.cjs(RWI_CH), RWI.sex = RWI.Sex,
                  y.WMI = WMI_CH, n.occasions.WMI = ncol(WMI_CH), nind.WMI = nind.WMI, f.WMI = f.WMI, z.WMI = known.state.cjs(WMI_CH), WMI.sex = WMI.Sex,
                  y.VON = VON_CH, n.occasions.VON = ncol(VON_CH), nind.VON = nind.VON, f.VON = f.VON, z.VON = known.state.cjs(VON_CH), VON.sex = VON.Sex,
                  y.BGU = BGU_CH, n.occasions.BGU = ncol(BGU_CH), nind.BGU = nind.BGU, f.BGU = f.BGU, z.BGU = known.state.cjs(BGU_CH), BGU.sex = BGU.Sex,
                  y.BHO = BHO_CH, n.occasions.BHO = ncol(BHO_CH), nind.BHO = nind.BHO, f.BHO = f.BHO, z.BHO = known.state.cjs(BHO_CH), BHO.sex = BHO.Sex,
                  y.CNI = CNI_CH, n.occasions.CNI = ncol(CNI_CH), nind.CNI = nind.CNI, f.CNI = f.CNI, z.CNI = known.state.cjs(CNI_CH), CNI.sex = CNI.Sex,
                  y.VMA = VMA_CH, n.occasions.VMA = ncol(VMA_CH), nind.VMA = nind.VMA, f.VMA = f.VMA, z.VMA = known.state.cjs(VMA_CH), VMA.sex = VMA.Sex,
                  y.VONW = VONW_CH, n.occasions.VONW = ncol(VONW_CH), nind.VONW = nind.VONW, f.VONW = f.VONW, z.VONW = known.state.cjs(VONW_CH), VONW.sex = VONW.Sex,
                  nsite = 11)


## Grab initial values
inits <- function(){list(
)}


## Organize parameters monitored (included trace parameters for appendix figure)
parameters <- c('overall.female.survival.gl', 'overall.male.survival.gl',
                'overall.f.CCR.phi', 'overall.f.MMA.phi', 'overall.f.RMI.phi',
                'overall.f.RWI.phi', 'overall.f.WMI.phi', 'overall.f.VON.phi',
                'overall.f.BGU.phi', 'overall.f.BHO.phi',
                'overall.f.CNI.phi', 'overall.f.VMA.phi', 'overall.m.CCR.phi',
                'overall.m.MMA.phi', 'overall.m.RMI.phi', 'overall.m.RWI.phi',
                'overall.m.WMI.phi', 'overall.m.VON.phi',
                'overall.m.BGU.phi', 'overall.m.BHO.phi', 'overall.m.CNI.phi',
                'overall.m.VMA.phi',
                'overall.m.VONW.phi',
                'phi.CCR', 'p.CCR', 'delta.phi.CCR', 'CCR.beta.sex', 'overall.CCR.phi', 'overall.delta.phi.CCR', 'overall.CCR.p',
                'phi.MMA', 'p.MMA', 'delta.phi.MMA', 'MMA.beta.sex', 'overall.MMA.phi', 'overall.delta.phi.MMA', 'overall.MMA.p',
                'phi.RMI', 'p.RMI', 'delta.phi.RMI', 'RMI.beta.sex', 'overall.RMI.phi', 'overall.delta.phi.RMI', 'overall.RMI.p',
                'phi.RWI', 'p.RWI', 'delta.phi.RWI', 'RWI.beta.sex', 'overall.RWI.phi', 'overall.delta.phi.RWI', 'overall.RWI.p',
                'phi.WMI', 'p.WMI', 'delta.phi.WMI', 'WMI.beta.sex', 'overall.WMI.phi', 'overall.delta.phi.WMI', 'overall.WMI.p',
                'phi.VON', 'p.VON', 'delta.phi.VON', 'VON.beta.sex', 'overall.VON.phi', 'overall.delta.phi.VON', 'overall.VON.p',
                'phi.BGU', 'p.BGU', 'delta.phi.BGU', 'BGU.beta.sex', 'overall.BGU.phi', 'overall.delta.phi.BGU', 'overall.BGU.p',
                'phi.BHO', 'p.BHO', 'delta.phi.BHO', 'BHO.beta.sex', 'overall.BHO.phi', 'overall.delta.phi.BHO', 'overall.BHO.p',
                'phi.CNI', 'p.CNI', 'delta.phi.CNI', 'CNI.beta.sex', 'overall.CNI.phi', 'overall.delta.phi.CNI', 'overall.CNI.p',
                'phi.VMA', 'p.VMA', 'delta.phi.VMA', 'VMA.beta.sex', 'overall.VMA.phi', 'overall.delta.phi.VMA', 'overall.VMA.p',
                'phi.VONW', 'p.VONW', 'delta.phi.VONW', 'VONW.beta.sex', 'overall.VONW.phi', 'overall.VONW.p',
                'mean.overall.survival.GL', 'mu.GL')

trace.parameters <- c('overall.female.survival.gl', 'overall.male.survival.gl',
                      'overall.CCR.phi', 'overall.MMA.phi', 'overall.RMI.phi', 'overall.RWI.phi', 'overall.WMI.phi',
                      'overall.VON.phi', 'overall.BGU.phi', 'overall.BHO.phi',
                      'overall.CNI.phi', 'overall.VMA.phi', 'overall.VONW.phi',
                      'mean.overall.survival.GL') 

## Identify MCMC settings
ni <- 50000
nt <- 50
nb <- 2500
nc <- 3

## Set working directory
setwd('/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Model')

## Run the model
GLMultiSiteNo2year <- jags(jags.data, inits, parameters, "GLMultiSiteNo2year.jags",
                           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)

## Summarize posteriors
print(GLMultiSiteNo2year, digits = 3)
