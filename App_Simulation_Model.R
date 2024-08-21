## This script runs a series of simulations to see under which scenarios we could 
## expect to increase precision in our parameters. Different scenarios included 
## running survival models that involved different occasion lengths 
## (3, 6, and 9 years), different levels of female representation in the datasets 
## (30% and 50%), different sample sizes of initial capture at each occasion 
## (10, 30, and 50), and different recapture probabilities (30%, 60%, and 90%).

## Load packages
library(IPMbook)
library(jagsUI)
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)
library(sjPlot)
library(MCMCvis)
library(truncnorm)

## Load in Appalachian JAGS model file
load("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Model/AppMultiSiteNo2year.jags")

################################################################################
###                              MODEL FUNCTIONS                             ###
################################################################################

## Define function to simulate a capture-history matrix
simul.cjs <- function(PHI, P, marked){
  n.occasions <- dim(PHI) [2] + 1
  CH <- matrix (0, ncol = n.occasions, nrow = sum(marked))
  # define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  # fill the CH matrix
  for (i in 1:sum(marked)){
    CH[i, mark.occ[i]] <- 1
    if (mark.occ[i] == n.occasions) next
    for (t in (mark.occ[i]+1): n.occasions){
      sur <- rbinom (1,1,PHI[i, t-1])
      if(sur==0) break
      rp <- rbinom(1,1,P[i,t-1])
      if(rp==1) CH[i,t] <- 1
    } #t
  } #i
  return(CH)
}

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

## Indicate simulation preferences
nsim <- 5 #s
n.occ <- c(3, 6, 9) #o  
sex.rat <- c(0.3, 0.5) #f
samp.siz <- c(10, 30, 50) #n
recap <- c(0.3, 0.6, 0.9) #p



for (s in 1:nsim) {
  # Calculate the iteration count within each set of models
  iteration_count <- ((s - 1) %% 30) + 1
  for (o in 1:length(n.occ)){
    for (f in 1:length(sex.rat)){
      for (n in 1:length(samp.siz)){
        for (p in 1:length(recap)){
          ###############################################################################
          ###                 ALL APPLICABLE METADATA FOR SIMULATION                  ###
          ###############################################################################
          
          # Number of capture occasions
          n.occasions <- n.occ[o]
          
          # Annual number of newly marked individuals
          m.marked <- rep(samp.siz[n]-round(samp.siz[n]*sex.rat[f]), n.occasions-1)
          f.marked <- rep(samp.siz[n]-round(samp.siz[n]*(1-sex.rat[f])), n.occasions-1)
          
          marked <- rbind(m.marked, f.marked)
          
          ###############################################################################
          ###                  SITE-SPECIFIC METADATA FOR SIMULATION                  ###
          ###############################################################################
          
          
          ################################   AWV    #####################################
          
          # Sex-specific survival
          AWV.phi.m <- rtruncnorm(n = n.occasions-1, mean = 0.520, sd = 0.041, a = 0, b = 0.999)
          AWV.phi.f <- rtruncnorm(n = n.occasions-1, mean = 0.500, sd = 0.047, a = 0, b = 0.999)
          AWV.p.m <- rep(recap[p], n.occasions-1)
          AWV.p.f <- rep(recap[p], n.occasions-1)
          
          ## Define matrices with survival and recapture probabilities
          AWV.PHI.M <- matrix(AWV.phi.m, ncol = n.occasions-1, nrow = sum(m.marked))
          AWV.PHI.F <- matrix(AWV.phi.f, ncol = n.occasions-1, nrow = sum(f.marked))
          AWV.P.M <- matrix(AWV.p.m, ncol = n.occasions-1, nrow = sum(m.marked))
          AWV.P.F <- matrix(AWV.p.f, ncol = n.occasions-1, nrow = sum(f.marked))
          
          ## Simulate capture histories
          AWV.CH.M <- simul.cjs (AWV.PHI.M, AWV.P.M, m.marked)
          AWV.CH.F <- simul.cjs (AWV.PHI.F, AWV.P.F, f.marked)
          
          ## Merge CH
          AWV_CH <- rbind(AWV.CH.M, AWV.CH.F)
          
          f.AWV <- apply(AWV_CH, 1, get.first)
          AWV.Sex <- c(rep(1, dim(AWV.CH.M)[1]), rep(2, dim(AWV.CH.F)[1]))
          AWV.Sex <- as.numeric(AWV.Sex)
          nyears.AWV <- n.occasions
          nind.AWV<-nrow(AWV_CH)
          
          ################################   AWV    #####################################
          
          ################################   CNY    #####################################
          
          # Sex-specific survival
          CNY.phi.m <- rtruncnorm(n = n.occasions-1, mean = 0.562, sd = 0.049, a = 0, b = 0.999)
          CNY.phi.f <- rtruncnorm(n = n.occasions-1, mean = 0.541, sd = 0.053, a = 0, b = 0.999)
          CNY.p.f <- CNY.p.m  <- AWV.p.f
          
          ## Define matrices with survival and recapture probabilities
          CNY.PHI.M <- matrix(CNY.phi.m, ncol = n.occasions-1, nrow = sum(m.marked))
          CNY.PHI.F <- matrix(CNY.phi.f, ncol = n.occasions-1, nrow = sum(f.marked))
          CNY.P.M <- matrix(CNY.p.m, ncol = n.occasions-1, nrow = sum(m.marked))
          CNY.P.F <- matrix(CNY.p.f, ncol = n.occasions-1, nrow = sum(f.marked))
          
          ## Simulate capture histories
          CNY.CH.M <- simul.cjs (CNY.PHI.M, CNY.P.M, m.marked)
          CNY.CH.F <- simul.cjs (CNY.PHI.F, CNY.P.F, f.marked)
          
          ## Merge CH
          CNY_CH <- rbind(CNY.CH.M, CNY.CH.F)
          
          f.CNY <- apply(CNY_CH, 1, get.first)
          CNY.Sex <- c(rep(1, dim(CNY.CH.M)[1]), rep(2, dim(CNY.CH.F)[1]))
          CNY.Sex <- as.numeric(CNY.Sex)
          nyears.CNY <- n.occasions
          nind.CNY<-nrow(CNY_CH)
          
          ################################   CNY    #####################################
          
          ################################   KNC    #####################################
          
          # Sex-specific survival
          KNC.phi.m <- rtruncnorm(n = n.occasions-1, mean = 0.392, sd = 0.098, a = 0, b = 0.999)
          KNC.phi.f <- rtruncnorm(n = n.occasions-1, mean = 0.438, sd = 0.131, a = 0, b = 0.999)
          KNC.p.m  <- KNC.p.f <- AWV.p.f
          
          ## Define matrices with survival and recapture probabilities
          KNC.PHI.M <- matrix(KNC.phi.m, ncol = n.occasions-1, nrow = sum(m.marked))
          KNC.PHI.F <- matrix(KNC.phi.f, ncol = n.occasions-1, nrow = sum(f.marked))
          KNC.P.M <- matrix(KNC.p.m, ncol = n.occasions-1, nrow = sum(m.marked))
          KNC.P.F <- matrix(KNC.p.f, ncol = n.occasions-1, nrow = sum(f.marked))

          ## Simulate capture histories
          KNC.CH.M <- simul.cjs (KNC.PHI.M, KNC.P.M, m.marked)
          KNC.CH.F <- simul.cjs (KNC.PHI.F, KNC.P.F, f.marked)
          
          ## Merge CH
          KNC_CH <- rbind(KNC.CH.M, KNC.CH.F)
          
          f.KNC <- apply(KNC_CH, 1, get.first)
          KNC.Sex <- c(rep(1, dim(KNC.CH.M)[1]), rep(2, dim(KNC.CH.F)[1]))
          KNC.Sex <- as.numeric(KNC.Sex)
          nyears.KNC <- n.occasions
          nind.KNC<-nrow(KNC_CH)
          
          ################################   KNC    #####################################
          
          ################################   PAE    #####################################
          
          # Sex-specific survival
          PAE.phi.m <- rtruncnorm(n = n.occasions-1, mean = 0.480, sd = 0.075, a = 0, b = 0.999)
          PAE.phi.f <- rtruncnorm(n = n.occasions-1, mean = 0.494, sd = 0.040, a = 0, b = 0.999)
          PAE.p.m  <- PAE.p.f  <- AWV.p.f
          
          ## Define matrices with survival and recapture probabilities
          PAE.PHI.M <- matrix(PAE.phi.m, ncol = n.occasions-1, nrow = sum(m.marked))
          PAE.PHI.F <- matrix(PAE.phi.f, ncol = n.occasions-1, nrow = sum(f.marked))
          PAE.P.M <- matrix(PAE.p.m, ncol = n.occasions-1, nrow = sum(m.marked))
          PAE.P.F <- matrix(PAE.p.f, ncol = n.occasions-1, nrow = sum(f.marked))

          ## Simulate capture histories
          PAE.CH.M <- simul.cjs (PAE.PHI.M, PAE.P.M, m.marked)
          PAE.CH.F <- simul.cjs (PAE.PHI.F, PAE.P.F, f.marked)
          
          ## Merge CH
          PAE_CH <- rbind(PAE.CH.M, PAE.CH.F)

          f.PAE <- apply(PAE_CH, 1, get.first)
          PAE.Sex <- c(rep(1, dim(PAE.CH.M)[1]), rep(2, dim(PAE.CH.F)[1]))
          PAE.Sex <- as.numeric(PAE.Sex)
          nyears.PAE <- n.occasions
          nind.PAE<-nrow(PAE_CH)
          
          ################################   PAE    #####################################
          
          ################################   PAC    #####################################
          
          # Sex-specific survival
          PAC.phi.m <- rtruncnorm(n = n.occasions-1, mean = 0.757, sd = 0.084, a = 0, b = 0.999)
          PAC.phi.f <- rtruncnorm(n = n.occasions-1, mean = 0.608, sd = 0.144, a = 0, b = 0.999)
          PAC.p.m <- PAC.p.f <- AWV.p.f

          ## Define matrices with survival and recapture probabilities
          PAC.PHI.M <- matrix(PAC.phi.m, ncol = n.occasions-1, nrow = sum(m.marked))
          PAC.PHI.F <- matrix(PAC.phi.f, ncol = n.occasions-1, nrow = sum(f.marked))
          PAC.P.M <- matrix(PAC.p.m, ncol = n.occasions-1, nrow = sum(m.marked))
          PAC.P.F <- matrix(PAC.p.f, ncol = n.occasions-1, nrow = sum(f.marked))

          ## Simulate capture histories
          PAC.CH.M <- simul.cjs (PAC.PHI.M, PAC.P.M, m.marked)
          PAC.CH.F <- simul.cjs (PAC.PHI.F, PAC.P.F, f.marked)
          
          ## Merge CH
          PAC_CH <- rbind(PAC.CH.M, PAC.CH.F)

          f.PAC <- apply(PAC_CH, 1, get.first)
          PAC.Sex <- c(rep(1, dim(PAC.CH.M)[1]), rep(2, dim(PAC.CH.F)[1]))
          PAC.Sex <- as.numeric(PAC.Sex)
          nyears.PAC <- n.occasions
          nind.PAC<-nrow(PAC_CH)
          
          ################################   PAC    #####################################
          
          ################################   SNC    #####################################
          
          # Sex-specific survival
          SNC.phi.m <- rtruncnorm(n = n.occasions-1, mean = 0.518, sd = 0.095, a = 0, b = 0.999)
          SNC.phi.f <- rtruncnorm(n = n.occasions-1, mean = 0.478, sd = 0.102, a = 0, b = 0.999)
          SNC.p.m <- SNC.p.f <- AWV.p.f
          
          ## Define matrices with survival and recapture probabilities
          SNC.PHI.M <- matrix(SNC.phi.m, ncol = n.occasions-1, nrow = sum(m.marked))
          SNC.PHI.F <- matrix(SNC.phi.f, ncol = n.occasions-1, nrow = sum(f.marked))
          SNC.P.M <- matrix(SNC.p.m, ncol = n.occasions-1, nrow = sum(m.marked))
          SNC.P.F <- matrix(SNC.p.f, ncol = n.occasions-1, nrow = sum(f.marked))

          ## Simulate capture histories
          SNC.CH.M <- simul.cjs (SNC.PHI.M, SNC.P.M, m.marked)
          SNC.CH.F <- simul.cjs (SNC.PHI.F, SNC.P.F, f.marked)
          
          ## Merge CH
          SNC_CH <- rbind(SNC.CH.M, SNC.CH.F)
          
          f.SNC <- apply(SNC_CH, 1, get.first)
          SNC.Sex <- c(rep(1, dim(SNC.CH.M)[1]), rep(2, dim(SNC.CH.F)[1]))
          SNC.Sex <- as.numeric(SNC.Sex)
          nyears.SNC <- n.occasions
          nind.SNC<-nrow(SNC_CH)
          
          ################################   SNC    #####################################
          
          
          ################################   CWV    #####################################
          
          # Sex-specific survival
          CWV.phi.m <- rtruncnorm(n = n.occasions-1, mean = 0.422, sd = 0.054, a = 0, b = 0.999)
          CWV.phi.f <- rtruncnorm(n = n.occasions-1, mean = 0.473, sd = 0.082, a = 0, b = 0.999)
          CWV.p.m <- CWV.p.f <- AWV.p.f
          
          ## Define matrices with survival and recapture probabilities
          CWV.PHI.M <- matrix(CWV.phi.m, ncol = n.occasions-1, nrow = sum(m.marked))
          CWV.PHI.F <- matrix(CWV.phi.f, ncol = n.occasions-1, nrow = sum(f.marked))
          CWV.P.M <- matrix(CWV.p.m, ncol = n.occasions-1, nrow = sum(m.marked))
          CWV.P.F <- matrix(CWV.p.f, ncol = n.occasions-1, nrow = sum(f.marked))

          ## Simulate capture histories
          CWV.CH.M <- simul.cjs (CWV.PHI.M, CWV.P.M, m.marked)
          CWV.CH.F <- simul.cjs (CWV.PHI.F, CWV.P.F, f.marked)
          
          ## Merge CH
          CWV_CH <- rbind(CWV.CH.M, CWV.CH.F)
          
          f.CWV <- apply(CWV_CH, 1, get.first)
          CWV.Sex <- c(rep(1, dim(CWV.CH.M)[1]), rep(2, dim(CWV.CH.F)[1]))
          CWV.Sex <- as.numeric(CWV.Sex)
          nyears.CWV <- n.occasions
          nind.CWV<-nrow(CWV_CH)
          
          ################################   CWV    #####################################
          
          ################################   BTN    #####################################
          
          # Sex-specific survival
          BTN.phi.m <- rtruncnorm(n = n.occasions-1, mean = 0.487, sd = 0.063, a = 0, b = 0.999)
          BTN.phi.f <- rtruncnorm(n = n.occasions-1, mean = 0.418, sd = 0.084, a = 0, b = 0.999)
          BTN.p.m <- BTN.p.f <- AWV.p.f
          
          ## Define matrices with survival and recapture probabilities
          BTN.PHI.M <- matrix(BTN.phi.m, ncol = n.occasions-1, nrow = sum(m.marked))
          BTN.PHI.F <- matrix(BTN.phi.f, ncol = n.occasions-1, nrow = sum(f.marked))
          BTN.P.M <- matrix(BTN.p.m, ncol = n.occasions-1, nrow = sum(m.marked))
          BTN.P.F <- matrix(BTN.p.f, ncol = n.occasions-1, nrow = sum(f.marked))
          
          ## Simulate capture histories
          BTN.CH.M <- simul.cjs (BTN.PHI.M, BTN.P.M, m.marked)
          BTN.CH.F <- simul.cjs (BTN.PHI.F, BTN.P.F, f.marked)
          
          ## Merge CH
          BTN_CH <- rbind(BTN.CH.M, BTN.CH.F)
          
          f.BTN <- apply(BTN_CH, 1, get.first)
          BTN.Sex <- c(rep(1, dim(BTN.CH.M)[1]), rep(2, dim(BTN.CH.F)[1]))
          BTN.Sex <- as.numeric(BTN.Sex)
          nyears.BTN <- n.occasions
          nind.BTN<-nrow(BTN_CH)
          
          ################################   BTN    #####################################
          
          ####################organize data for all sub-models
          jags.data <- list(y.AWV = AWV_CH, n.occasions.AWV = ncol(AWV_CH), nind.AWV = nind.AWV, f.AWV = f.AWV, z.AWV = known.state.cjs(AWV_CH), AWV.sex = AWV.Sex,
                            y.CNY = CNY_CH, n.occasions.CNY = ncol(CNY_CH), nind.CNY = nind.CNY, f.CNY = f.CNY, z.CNY = known.state.cjs(CNY_CH), CNY.sex = CNY.Sex,
                            y.KNC = KNC_CH, n.occasions.KNC = ncol(KNC_CH), nind.KNC = nind.KNC, f.KNC = f.KNC, z.KNC = known.state.cjs(KNC_CH), KNC.sex = KNC.Sex,
                            y.PAE = PAE_CH, n.occasions.PAE = ncol(PAE_CH), nind.PAE = nind.PAE, f.PAE = f.PAE, z.PAE = known.state.cjs(PAE_CH),
                            y.PAC = PAC_CH, n.occasions.PAC = ncol(PAC_CH), nind.PAC = nind.PAC, f.PAC = f.PAC, z.PAC = known.state.cjs(PAC_CH), PAC.sex = PAC.Sex,
                            y.SNC = SNC_CH, n.occasions.SNC = ncol(SNC_CH), nind.SNC = nind.SNC, f.SNC = f.SNC, z.SNC = known.state.cjs(SNC_CH), SNC.sex = SNC.Sex,
                            y.CWV = CWV_CH, n.occasions.CWV = ncol(CWV_CH), nind.CWV = nind.CWV, f.CWV = f.CWV, z.CWV = known.state.cjs(CWV_CH), CWV.sex = CWV.Sex,
                            y.BTN = BTN_CH, n.occasions.BTN = ncol(BTN_CH), nind.BTN = nind.BTN, f.BTN = f.BTN, z.BTN = known.state.cjs(BTN_CH), BTN.sex = BTN.Sex,
                            nsite = 8)
          
          ## initial values
          
          inits <- function(){list(
          )}
          
          
          ## parameters monitored ########
          parameters <- c('overall.female.survival.app', 'overall.male.survival.app',
                          'overall.AWV.phi', 'overall.f.AWV.phi', 'overall.m.AWV.phi', 'overall.delta.phi.AWV', 'overall.AWV.p',
                          'overall.CNY.phi', 'overall.f.CNY.phi', 'overall.m.CNY.phi', 'overall.delta.phi.CNY', 'overall.CNY.p',
                          'overall.KNC.phi', 'overall.f.KNC.phi', 'overall.m.KNC.phi', 'overall.delta.phi.KNC', 'overall.KNC.p',
                          'overall.PAE.phi', 'overall.f.PAE.phi', 'overall.m.PAE.phi', 'overall.delta.phi.PAE', 'overall.PAE.p',
                          'overall.PAC.phi', 'overall.f.PAC.phi', 'overall.m.PAC.phi', 'overall.delta.phi.PAC', 'overall.PAC.p',
                          'overall.SNC.phi', 'overall.f.SNC.phi', 'overall.m.SNC.phi', 'overall.delta.phi.SNC', 'overall.SNC.p',
                          'overall.CWV.phi', 'overall.f.CWV.phi', 'overall.m.CWV.phi', 'overall.delta.phi.CWV', 'overall.CWV.p',
                          'overall.BTN.phi', 'overall.f.BTN.phi', 'overall.m.BTN.phi', 'overall.delta.phi.BTN', 'overall.BTN.p',
                          'mean.overall.survival.App', 'mu.App')
          
          
          ## MCMC settings
          ni <- 10000
          nt <- 25
          nb <- 2500
          nc <- 3
          
          ##################### Run and summarize the model ##########################
          
          setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Model")
          ## call jags from R
          simulAppMultiSite <- jags(jags.data, inits, parameters, "AppMultiSiteNo2year.jags",
                                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T)
          
          
          save(simulAppMultiSite, file=paste0('/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Simulations/App Simulation Rnorm/sim_nocc', n.occ[o], '_sr', sex.rat[f], '_n', samp.siz[n], '_rc', recap[p], '_', iteration_count, '.RData'))  
          # print(paste(s, o, f, n, p))  
          
        }
      }
    }
  }
  
}
