## Figure A2 in manuscript. This script runs the traceplots for the Appalachian
## apparent annual survival model. Traceplots will show successful convergence 
## for all parameters. 

## Load packages
library(IPMbook)
library(jagsUI)
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)
library(sjPlot)
library(ggpubr)

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
comboCH <- read.csv("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Script/All_App.csv")

## Separate capture histories into specific sites
AWV<-subset(comboCH, Site == "AWV") ## Northeastern West Virginia
CNY<-subset(comboCH, Site == "CNY") ## Southern New York
KNC<-subset(comboCH, Site == "KNC") ## Western North Carolina
PAE<-subset(comboCH, Site == "PAE") ## Eastern Pennsylvania
PAC<-subset(comboCH, Site == "PAC") ## Central Pennsylvania
SNC<-subset(comboCH, Site == "SNC") ## Northwestern North Carolina
CWV<-subset(comboCH, Site == "CWV") ## Southern West Virginia
BTN<-subset(comboCH, Site == "BTN") ## Tennessee

## Take out individual with unknown sex in Tennessee population
BTN <- subset(BTN, Sex==1 | Sex==2) 

## Pull out capture histories for each site
AWV_CH <- AWV[,9:16]
CNY_CH <- CNY[,3:11]
KNC_CH <- KNC[,17:22]
PAE_CH <- PAE[,13:16]
PAC_CH <- PAC[,9:12]
SNC_CH <- SNC[,11:16]
CWV_CH <- CWV[,1:22]
BTN_CH <- BTN[,4:12]

## Get first occassion for each site
f.AWV <- apply(AWV_CH, 1, get.first)
f.CNY <- apply(CNY_CH, 1, get.first)
f.KNC <- apply(KNC_CH, 1, get.first)
f.PAE <- apply(PAE_CH, 1, get.first)
f.PAC <- apply(PAC_CH, 1, get.first)
f.SNC <- apply(SNC_CH, 1, get.first)
f.CWV <- apply(CWV_CH, 1, get.first)
f.BTN <- apply(BTN_CH, 1, get.first)

## Pull out sex covariate for each site
AWV.Sex <- AWV[,25]
CNY.Sex <- CNY[,25]
KNC.Sex <- KNC[,25]
PAE.Sex <- PAE[,25]
PAC.Sex <- PAC[,25]
SNC.Sex <- SNC[,25]
CWV.Sex <- CWV[,25]
BTN.Sex <- BTN[,25]

## Convert all sex covariates to numeric
AWV.Sex <- as.numeric(AWV.Sex)
CNY.Sex <- as.numeric(CNY.Sex)
KNC.Sex <- as.numeric(KNC.Sex)
PAE.Sex <- as.numeric(PAE.Sex)
PAC.Sex <- as.numeric(PAC.Sex)
SNC.Sex <- as.numeric(SNC.Sex)
CWV.Sex <- as.numeric(CWV.Sex)
BTN.Sex <- as.numeric(BTN.Sex)

## Create vector of number of individuals for each site
nind.AWV<-nrow(AWV_CH)
nind.CNY<-nrow(CNY_CH)
nind.KNC<-nrow(KNC_CH)
nind.PAE<-nrow(PAE_CH)
nind.PAC<-nrow(PAC_CH)
nind.SNC<-nrow(SNC_CH)
nind.CWV<-nrow(CWV_CH)
nind.BTN<-nrow(BTN_CH)

## Create vector of number of years for each site
nyears.AWV <-length(2008:2015)
nyears.CNY <-length(2002:2010)
nyears.KNC <-length(2016:2021)
nyears.PAE <-length(2012:2015)
nyears.PAC <-length(2008:2011)
nyears.SNC <-length(2010:2015)
nyears.CWV <-length(2000:2021)
nyears.BTN <-length(2003:2011)

## Find the sum of both sexes
app_sex_count <- rbind(table(AWV.Sex), table(CNY.Sex), table(KNC.Sex), table(PAE.Sex),
                       table(PAC.Sex), table(SNC.Sex), table(CWV.Sex),
                       table(BTN.Sex))
colnames(app_sex_count) <- c("males", "females")

app_sex_labels <- c("AWV", "CNY", "KNC", "PAE", "PAC", "SNC", "CWV", "BTN")
app_sex_count <- cbind(app_sex_count, app_sex_labels)

app_sex_count <- as.data.frame(app_sex_count)
app_sex_count$males <- as.numeric(app_sex_count$males)
app_sex_count$females <- as.numeric(app_sex_count$females)


################################################################################
###                             BUILD THE MODEL                              ###
################################################################################

sink("AppMultiSiteNo2year.jags")
cat("

model {

  ##############################################
  # Define model likelihood for mean survival of a site
  for (s in 1:nsite){
    # s = site effect
    mu.App[s] ~ dnorm(overall.survival.App, tau.site.App)
  }
 
  # Define prior for overall mean survival

  overall.survival.App ~ dnorm(0,1.0E-3)

  #overall.survival.App ~ dbeta(1,1)
  logit(mean.overall.survival.App) <- overall.survival.App
  # Define prior for inter-site variance (site random effect)
 
  tau.site.App <- 1 / (sd.site.App * sd.site.App)
  sd.site.App ~ dunif(0,5)
 
  # Monitor variance
 
  var.site.App <- 1/tau.site.App
 
############################## START AWV #######################################

# Priors and constraints
 # for(i in 1:nind.AWV){
   for (t in 1:(n.occasions.AWV-1)){ # time dependent survival
    logit(phi.AWV[1,t]) <- logitphi.AWV[1,t]
    logit(phi.AWV[2,t]) <- logitphi.AWV[2,t]
    logitphi.AWV[1,t] ~ dnorm(mu.App[1],tau.eps.AWV)
    logitphi.AWV[2,t] ~ dnorm(mu.App[1],tau.eps.AWV)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.AWV[t]-AWV.beta.sex*AWV.sex[i])) + AWV.beta.sex*AWV.sex[i]
 
    p.AWV[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.AWV-1)){
  delta.phi.AWV[t] <- phi.AWV[1,t]-phi.AWV[2,t]
  #logitphi.AWV[t] ~ dnorm(mu.App[1],tau.eps.AWV)
  mean.m.AWV.phi[t] <- mean(phi.AWV[1,t])
  mean.f.AWV.phi[t] <- mean(phi.AWV[2,t])
  }
 
overall.f.AWV.phi <- mean(mean.f.AWV.phi)
overall.m.AWV.phi <- mean(mean.m.AWV.phi)
overall.AWV.phi <- mean(phi.AWV)
overall.delta.phi.AWV <- mean(delta.phi.AWV)

  #AWV.beta.sex ~ dnorm(0, 0.001)I(-10, 10) # Prior for slope parameter
 
  # Define prior for variance of temporal random effect
  tau.eps.AWV <- 1 / (sd.eps.AWV*sd.eps.AWV)
  sd.eps.AWV ~ dunif(0,5)
  # Monitor variance
  var.eps.AWV <- 1 / tau.eps.AWV

# Likelihood
for (i in 1:nind.AWV){
   # Define latent state at first capture
   z.AWV[i,f.AWV[i]] ~ dbern(1)
   for (t in (f.AWV[i]+1):n.occasions.AWV){
      # State process
      z.AWV[i,t] ~  dbern(phi.AWV[AWV.sex[i], t-1] * z.AWV[i,t-1])
      #z.AWV[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.AWV[i, t-1] * z.AWV[i,t-1])
      # Observation process
      y.AWV[i,t] ~ dbern(p.AWV[t-1] * z.AWV[i,t])
      #y.AWV[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.AWV[i,t-1] * z.AWV[i,t])
      } #t
   } #i

################################ END AWV #######################################
 
############################## START CNY #######################################

# Priors and constraints
 # for(i in 1:nind.CNY){
   for (t in 1:(n.occasions.CNY-1)){ # time dependent survival
    logit(phi.CNY[1,t]) <- logitphi.CNY[1,t]
    logit(phi.CNY[2,t]) <- logitphi.CNY[2,t]
    logitphi.CNY[1,t] ~ dnorm(mu.App[2],tau.eps.CNY)
    logitphi.CNY[2,t] ~ dnorm(mu.App[2],tau.eps.CNY)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.CNY[t]-CNY.beta.sex*CNY.sex[i])) # + CNY.beta.sex*CNY.sex[i]
 
    p.CNY[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.CNY-1)){
  delta.phi.CNY[t] <- phi.CNY[1,t]-phi.CNY[2,t]
  #logitphi.CNY[t] ~ dnorm(mu.App[2],tau.eps.CNY)
  mean.m.CNY.phi[t] <- mean(phi.CNY[1,t])
  mean.f.CNY.phi[t] <- mean(phi.CNY[2,t])
  }
 
overall.f.CNY.phi <- mean(mean.f.CNY.phi)
overall.m.CNY.phi <- mean(mean.m.CNY.phi)
overall.CNY.phi <- mean(phi.CNY)
overall.delta.phi.CNY <- mean(delta.phi.CNY)
 
  # Define prior for variance of temporal random effect
  tau.eps.CNY <- 1 / (sd.eps.CNY*sd.eps.CNY)
  sd.eps.CNY ~ dunif(0,5)
  # Monitor variance
  var.eps.CNY <- 1 / tau.eps.CNY

# Likelihood
for (i in 1:nind.CNY){
   # Define latent state at first capture
   z.CNY[i,f.CNY[i]] ~ dbern(1)
   for (t in (f.CNY[i]+1):n.occasions.CNY){
      # State process
      z.CNY[i,t] ~  dbern(phi.CNY[CNY.sex[i], t-1] * z.CNY[i,t-1])
      #z.CNY[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.CNY[i, t-1] * z.CNY[i,t-1])
      # Observation process
      y.CNY[i,t] ~ dbern(p.CNY[t-1] * z.CNY[i,t])
      #y.CNY[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.CNY[i,t-1] * z.CNY[i,t])
      } #t
   } #i

################################ END CNY #######################################

############################## START KNC #######################################

# Priors and constraints
 # for(i in 1:nind.KNC){
   for (t in 1:(n.occasions.KNC-1)){ # time dependent survival
    logit(phi.KNC[1,t]) <- logitphi.KNC[1,t]
    logit(phi.KNC[2,t]) <- logitphi.KNC[2,t]
    logitphi.KNC[1,t] ~ dnorm(mu.App[3],tau.eps.KNC)
    logitphi.KNC[2,t] ~ dnorm(mu.App[3],tau.eps.KNC)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.KNC[t]-KNC.beta.sex*KNC.sex[i])) # + KNC.beta.sex*KNC.sex[i]
 
    p.KNC[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.KNC-1)){
  delta.phi.KNC[t] <- phi.KNC[1,t]-phi.KNC[2,t]
  #logitphi.KNC[t] ~ dnorm(mu.App[3],tau.eps.KNC)
  mean.m.KNC.phi[t] <- mean(phi.KNC[1,t])
  mean.f.KNC.phi[t] <- mean(phi.KNC[2,t])
  }
 
overall.f.KNC.phi <- mean(mean.f.KNC.phi)
overall.m.KNC.phi <- mean(mean.m.KNC.phi)
overall.KNC.phi <- mean(phi.KNC)
overall.delta.phi.KNC <- mean(delta.phi.KNC)

  # Define prior for variance of temporal random effect
  tau.eps.KNC <- 1 / (sd.eps.KNC*sd.eps.KNC)
  sd.eps.KNC ~ dunif(0,5)
  # Monitor variance
  var.eps.KNC <- 1 / tau.eps.KNC

# Likelihood
for (i in 1:nind.KNC){
   # Define latent state at first capture
   z.KNC[i,f.KNC[i]] ~ dbern(1)
   for (t in (f.KNC[i]+1):n.occasions.KNC){
      # State process
      z.KNC[i,t] ~  dbern(phi.KNC[KNC.sex[i], t-1] * z.KNC[i,t-1])
      #z.KNC[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.KNC[i, t-1] * z.KNC[i,t-1])
      # Observation process
      y.KNC[i,t] ~ dbern(p.KNC[t-1] * z.KNC[i,t])
      #y.KNC[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.KNC[i,t-1] * z.KNC[i,t])
      } #t
   } #i

################################ END CNY #######################################

############################## START PAE #######################################

# Priors and constraints

   for (t in 1:(n.occasions.PAE-1)){ # time dependent survival
    logit(phi.PAE[t]) <- logitphi.PAE[t]
    logitphi.PAE[t] ~ dnorm(mu.App[4],tau.eps.PAE)
 
    p.PAE[t] ~ dbeta(1,1)
    #epsilon[t] ~ dnorm(0, tau)
      } #t

overall.m.PAE.phi <- mean(phi.PAE)
overall.PAE.phi <- mean(phi.PAE)

  # Define prior for variance of temporal random effect
  tau.eps.PAE <- 1 / (sd.eps.PAE*sd.eps.PAE)
  sd.eps.PAE ~ dunif(0,5)
  # Monitor variance
  var.eps.PAE <- 1 / tau.eps.PAE

# Likelihood
for (i in 1:nind.PAE){
   # Define latent state at first capture
   z.PAE[i,f.PAE[i]] ~ dbern(1)
   for (t in (f.PAE[i]+1):n.occasions.PAE){
      # State process
      z.PAE[i,t] ~ dbern(phi.PAE[t-1] * z.PAE[i,t-1])
      # Observation process
      y.PAE[i,t] ~ dbern(p.PAE[t-1] * z.PAE[i,t])
      } #t
   } #i

################################ END PAE #######################################

############################## START PAC #######################################

# Priors and constraints
 # for(i in 1:nind.PAC){
   for (t in 1:(n.occasions.PAC-1)){ # time dependent survival
    logit(phi.PAC[1,t]) <- logitphi.PAC[1,t]
    logit(phi.PAC[2,t]) <- logitphi.PAC[2,t]
    logitphi.PAC[1,t] ~ dnorm(mu.App[5],tau.eps.PAC)
    logitphi.PAC[2,t] ~ dnorm(mu.App[5],tau.eps.PAC)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.PAC[t]-PAC.beta.sex*PAC.sex[i])) # + PAC.beta.sex*PAC.sex[i]
 
    p.PAC[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.PAC-1)){
  delta.phi.PAC[t] <- phi.PAC[1,t]-phi.PAC[2,t]
  #logitphi.PAC[t] ~ dnorm(mu.App[5],tau.eps.PAC)
  mean.m.PAC.phi[t] <- mean(phi.PAC[1,t])
  mean.f.PAC.phi[t] <- mean(phi.PAC[2,t])
  }
 
overall.f.PAC.phi <- mean(mean.f.PAC.phi)
overall.m.PAC.phi <- mean(mean.m.PAC.phi)
overall.PAC.phi <- mean(phi.PAC)
overall.delta.phi.PAC <- mean(delta.phi.PAC)
 
  # Define prior for variance of temporal random effect
  tau.eps.PAC <- 1 / (sd.eps.PAC*sd.eps.PAC)
  sd.eps.PAC ~ dunif(0,5)
  # Monitor variance
  var.eps.PAC <- 1 / tau.eps.PAC

# Likelihood
for (i in 1:nind.PAC){
   # Define latent state at first capture
   z.PAC[i,f.PAC[i]] ~ dbern(1)
   for (t in (f.PAC[i]+1):n.occasions.PAC){
      # State process
      z.PAC[i,t] ~  dbern(phi.PAC[PAC.sex[i], t-1] * z.PAC[i,t-1])
      #z.PAC[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.PAC[i, t-1] * z.PAC[i,t-1])
      # Observation process
      y.PAC[i,t] ~ dbern(p.PAC[t-1] * z.PAC[i,t])
      #y.PAC[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.PAC[i,t-1] * z.PAC[i,t])
      } #t
   } #i

################################ END PAC #######################################

############################## START SNC #######################################

# Priors and constraints
 # for(i in 1:nind.SNC){
   for (t in 1:(n.occasions.SNC-1)){ # time dependent survival
    logit(phi.SNC[1,t]) <- logitphi.SNC[1,t]
    logit(phi.SNC[2,t]) <- logitphi.SNC[2,t]
    logitphi.SNC[1,t] ~ dnorm(mu.App[6],tau.eps.SNC)
    logitphi.SNC[2,t] ~ dnorm(mu.App[6],tau.eps.SNC)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.SNC[t]-SNC.beta.sex*SNC.sex[i])) # + SNC.beta.sex*SNC.sex[i]
 
    p.SNC[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.SNC-1)){
  delta.phi.SNC[t] <- phi.SNC[1,t]-phi.SNC[2,t]
  #logitphi.SNC[t] ~ dnorm(mu.App[7],tau.eps.SNC)
  mean.m.SNC.phi[t] <- mean(phi.SNC[1,t])
  mean.f.SNC.phi[t] <- mean(phi.SNC[2,t])
  }
 
overall.f.SNC.phi <- mean(mean.f.SNC.phi)
overall.m.SNC.phi <- mean(mean.m.SNC.phi)
overall.SNC.phi <- mean(phi.SNC)
overall.delta.phi.SNC <- mean(delta.phi.SNC)

  # Define prior for variance of temporal random effect
  tau.eps.SNC <- 1 / (sd.eps.SNC*sd.eps.SNC)
  sd.eps.SNC ~ dunif(0,5)
  # Monitor variance
  var.eps.SNC <- 1 / tau.eps.SNC

# Likelihood
for (i in 1:nind.SNC){
   # Define latent state at first capture
   z.SNC[i,f.SNC[i]] ~ dbern(1)
   for (t in (f.SNC[i]+1):n.occasions.SNC){
      # State process
      z.SNC[i,t] ~  dbern(phi.SNC[SNC.sex[i], t-1] * z.SNC[i,t-1])
      #z.SNC[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.SNC[i, t-1] * z.SNC[i,t-1])
      # Observation process
      y.SNC[i,t] ~ dbern(p.SNC[t-1] * z.SNC[i,t])
      #y.SNC[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.SNC[i,t-1] * z.SNC[i,t])
      } #t
   } #i

################################ END SNC #######################################

############################## START CWV #######################################

# Priors and constraints
 # for(i in 1:nind.CWV){
   for (t in 1:(n.occasions.CWV-1)){ # time dependent survival
    logit(phi.CWV[1,t]) <- logitphi.CWV[1,t]
    logit(phi.CWV[2,t]) <- logitphi.CWV[2,t]
    logitphi.CWV[1,t] ~ dnorm(mu.App[7],tau.eps.CWV)
    logitphi.CWV[2,t] ~ dnorm(mu.App[7],tau.eps.CWV)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.CWV[t]-CWV.beta.sex*CWV.sex[i])) # + CWV.beta.sex*CWV.sex[i]
 
    p.CWV[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.CWV-1)){
  delta.phi.CWV[t] <- phi.CWV[1,t]-phi.CWV[2,t]
  #logitphi.CWV[t] ~ dnorm(mu.App[8],tau.eps.CWV)
  mean.m.CWV.phi[t] <- mean(phi.CWV[1,t])
  mean.f.CWV.phi[t] <- mean(phi.CWV[2,t])
  }
 
overall.f.CWV.phi <- mean(mean.f.CWV.phi)
overall.m.CWV.phi <- mean(mean.m.CWV.phi)
overall.CWV.phi <- mean(phi.CWV)
overall.delta.phi.CWV <- mean(delta.phi.CWV)

  # Define prior for variance of temporal random effect
  tau.eps.CWV <- 1 / (sd.eps.CWV*sd.eps.CWV)
  sd.eps.CWV ~ dunif(0,5)
  # Monitor variance
  var.eps.CWV <- 1 / tau.eps.CWV

# Likelihood
for (i in 1:nind.CWV){
   # Define latent state at first capture
   z.CWV[i,f.CWV[i]] ~ dbern(1)
   for (t in (f.CWV[i]+1):n.occasions.CWV){
      # State process
      z.CWV[i,t] ~  dbern(phi.CWV[CWV.sex[i], t-1] * z.CWV[i,t-1])
      #z.CWV[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.CWV[i, t-1] * z.CWV[i,t-1])
      # Observation process
      y.CWV[i,t] ~ dbern(p.CWV[t-1] * z.CWV[i,t])
      #y.CWV[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.CWV[i,t-1] * z.CWV[i,t])
      } #t
   } #i

################################ END CWV #######################################

############################## START BTN #######################################

# Priors and constraints
 # for(i in 1:nind.BTN){
   for (t in 1:(n.occasions.BTN-1)){ # time dependent survival
    logit(phi.BTN[1,t]) <- logitphi.BTN[1,t]
    logit(phi.BTN[2,t]) <- logitphi.BTN[2,t]
    logitphi.BTN[1,t] ~ dnorm(mu.App[8],tau.eps.BTN)
    logitphi.BTN[2,t] ~ dnorm(mu.App[8],tau.eps.BTN)
    #phi.est[i,t] <- 1 / (1+exp(-logitphi.BTN[t]-BTN.beta.sex*BTN.sex[i])) # + BTN.beta.sex*BTN.sex[i]
 
    p.BTN[t] ~ dbeta(1,1)

      } #t
 # }#i
 
  for (t in 1:(n.occasions.BTN-1)){
  delta.phi.BTN[t] <- phi.BTN[1,t]-phi.BTN[2,t]
  mean.m.BTN.phi[t] <- mean(phi.BTN[1,t])
  mean.f.BTN.phi[t] <- mean(phi.BTN[2,t])
  }
 
overall.f.BTN.phi <- mean(mean.f.BTN.phi)
overall.m.BTN.phi <- mean(mean.m.BTN.phi)
overall.BTN.phi <- mean(phi.BTN)
overall.delta.phi.BTN <- mean(delta.phi.BTN)

  # Define prior for variance of temporal random effect
  tau.eps.BTN <- 1 / (sd.eps.BTN*sd.eps.BTN)
  sd.eps.BTN ~ dunif(0,5)
  # Monitor variance
  var.eps.BTN <- 1 / tau.eps.BTN

# Likelihood
for (i in 1:nind.BTN){
   # Define latent state at first capture
   z.BTN[i,f.BTN[i]] ~ dbern(1)
   for (t in (f.BTN[i]+1):n.occasions.BTN){
      # State process
      z.BTN[i,t] ~  dbern(phi.BTN[BTN.sex[i], t-1] * z.BTN[i,t-1])
      #z.BTN[i,t] ~  dbern(mu1[i,t])
       #mu1[i,t] <- (phi.BTN[i, t-1] * z.BTN[i,t-1])
      # Observation process
      y.BTN[i,t] ~ dbern(p.BTN[t-1] * z.BTN[i,t])
      #y.BTN[i,t] ~ dbern(mu2[i,t])
       #mu2[i,t] <- (p.BTN[i,t-1] * z.BTN[i,t])
      } #t
   } #i

################################ END CWV #######################################

## Get overall male and female survival estimates
overall.female.survival <- mean(c(overall.f.AWV.phi, overall.f.CNY.phi, overall.f.KNC.phi, overall.f.PAC.phi,
  overall.f.SNC.phi,overall.f.CWV.phi,overall.f.BTN.phi))

overall.male.survival <- mean(c(overall.m.AWV.phi, overall.m.CNY.phi, overall.m.KNC.phi, overall.m.PAE.phi, overall.m.PAC.phi,
  overall.m.SNC.phi,overall.m.CWV.phi,overall.m.BTN.phi))
}"
,
fill = T)
sink()


################################################################################
###                               RUN THE MODEL                              ###
################################################################################

## Organize data for all sub-models
jags.data <- list(y.AWV = AWV_CH, n.occasions.AWV = ncol(AWV_CH), nind.AWV = nind.AWV, f.AWV = f.AWV, z.AWV = known.state.cjs(AWV_CH), AWV.sex = AWV.Sex,
                  y.CNY = CNY_CH, n.occasions.CNY = ncol(CNY_CH), nind.CNY = nind.CNY, f.CNY = f.CNY, z.CNY = known.state.cjs(CNY_CH), CNY.sex = CNY.Sex,
                  y.KNC = KNC_CH, n.occasions.KNC = ncol(KNC_CH), nind.KNC = nind.KNC, f.KNC = f.KNC, z.KNC = known.state.cjs(KNC_CH), KNC.sex = KNC.Sex,
                  y.PAE = PAE_CH, n.occasions.PAE = ncol(PAE_CH), nind.PAE = nind.PAE, f.PAE = f.PAE, z.PAE = known.state.cjs(PAE_CH),
                  y.PAC = PAC_CH, n.occasions.PAC = ncol(PAC_CH), nind.PAC = nind.PAC, f.PAC = f.PAC, z.PAC = known.state.cjs(PAC_CH), PAC.sex = PAC.Sex,
                  y.SNC = SNC_CH, n.occasions.SNC = ncol(SNC_CH), nind.SNC = nind.SNC, f.SNC = f.SNC, z.SNC = known.state.cjs(SNC_CH), SNC.sex = SNC.Sex,
                  y.CWV = CWV_CH, n.occasions.CWV = ncol(CWV_CH), nind.CWV = nind.CWV, f.CWV = f.CWV, z.CWV = known.state.cjs(CWV_CH), CWV.sex = CWV.Sex,
                  y.BTN = BTN_CH, n.occasions.BTN = ncol(BTN_CH), nind.BTN = nind.BTN, f.BTN = f.BTN, z.BTN = known.state.cjs(BTN_CH), BTN.sex = BTN.Sex,
                  nsite = 8)

## Grab initial values
inits <- function(){list(
)}

## Organize parameters monitored (included trace parameters for appendix figure)
trace.parameters <- c('overall.female.survival', 'overall.male.survival',
                      'overall.AWV.phi', 'overall.CNY.phi','overall.KNC.phi','overall.PAC.phi','overall.PAE.phi',
                      'overall.SNC.phi','overall.CWV.phi','overall.BTN.phi',
                      'mean.overall.survival.App') 

## Identify MCMC settings
ni <- 50000
nt <- 50
nb <- 2500
nc <- 3

## Set working directory
setwd('/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Model')

## Run the model
AppsMultiSiteNo2year <- jags(jags.data, inits, trace.parameters, "AppMultiSiteNo2year.jags",
                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)

## Summarize posteriors
print(AppsMultiSiteNo2year, digits = 3)

################################################################################
###                  CREATE TRACEPLOT VISUALIZATION                          ###
################################################################################

setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Tables and Figures/Figure A2/")

traceplot <- MCMCtrace(object = AppsMultiSiteNo2year$samples,
                       ISB = FALSE,
                       exact = TRUE,
                       pdf = TRUE)

dev.off()
