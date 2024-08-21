## Figure A7 in manuscript. This script shows the mean overall apparent annual
## survival (and 95% credible intervals) estimates for one site (West Virginia)
## as an example of what site-specific survival looks like under different
## simulation scenarios. Taken from mean posterior estimates of 30 simulations 
## for 54 different respective scenarios. 

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
library(tidyverse)
library(stringr)

## Indicate desired simulations 
nsim <- 30 #s
n.occ <- c(3, 6, 9) #o  
sex.rat <- c(0.3, 0.5) #f
samp.siz <- c(10, 30, 50) #n
recap <- c(0.3, 0.6, 0.9) #p

# Read model output in from file where simluations are stored
setwd('/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Simulations/App Simulation Rnorm') 

## Make a list of full file names
filenames <- list.files('.', full.names=T)

## Take the simulation number and ./ off of each file name
simnames <- unique(substr(filenames, start=3, stop=27)) 

## Indicate number of expected posteriors
nmcmc <- 1800

## Identify number of combinations of posteriors (54)
n_unique_datasets <- length(n.occ)*length(sex.rat)*length(samp.siz)*length(recap) 


################################################################################
###                               CREATE FUNCTIONS                           ###
################################################################################

#pull_samples <- function(q){
all_model_samples <- list() # should end up being 54 long
for (q in 1:length(simnames)){
  files_sub <- filenames[which(str_detect(filenames, pattern=simnames[q])==T)]
  mods <- list() # list to store models
  mat_param1 <- mat_param2 <- mat_param3 <- mat_param4 <- mat_param5 <- mat_param6 <- mat_param7 <- matrix(NA, nrow=nsim, ncol=nmcmc) # matrices to store - 1 row per model 1 col per posterior sample
  for (z in 1:nsim){ 
    load(files_sub[z]) # read in each model
    mods[[z]] <- simulAppMultiSite
    # fill parameter matrices with posterior samples for each model
    mat_param1[z,] <- mods[[z]]$sims.list$overall.delta.phi.AWV
    mat_param2[z,] <- mods[[z]]$sims.list$overall.delta.phi.CNY
    mat_param3[z,] <- mods[[z]]$sims.list$overall.delta.phi.KNC
    mat_param4[z,] <- mods[[z]]$sims.list$overall.delta.phi.PAC
    mat_param5[z,] <- mods[[z]]$sims.list$overall.delta.phi.SNC
    mat_param6[z,] <- mods[[z]]$sims.list$overall.delta.phi.CWV
    mat_param7[z,] <- mods[[z]]$sims.list$overall.delta.phi.BTN
  }
  samples <- list(AWV=mat_param1, CNY=mat_param2, KNC=mat_param3, PAC=mat_param4,
                  SNC=mat_param5, CWV=mat_param6, BTN=mat_param7) 
  all_model_samples[[q]] <- samples
}

names(all_model_samples) <- simnames # add names to elements of the list 

## 95% quanile of posterior means - comes from the posterior distribution - mean of means

# example of summarizing means (using all samples)
param_delta <- data.frame(simulation=names(all_model_samples), AWV=NA, CNY=NA, KNC=NA, 
                          PAC=NA, SNC=NA, CWV=NA, BTN=NA)
for (i in 1:length(all_model_samples)){
  param_delta$AWV[i] <- mean(all_model_samples[[i]]$AWV)
  param_delta$CNY[i] <- mean(all_model_samples[[i]]$CNY)
  param_delta$KNC[i] <- mean(all_model_samples[[i]]$KNC)
  param_delta$PAC[i] <- mean(all_model_samples[[i]]$PAC)
  param_delta$SNC[i] <- mean(all_model_samples[[i]]$SNC)
  param_delta$CWV[i] <- mean(all_model_samples[[i]]$CWV)
  param_delta$BTN[i] <- mean(all_model_samples[[i]]$BTN)
}
param_delta

# example of summarizing lower CI

param_lower <- data.frame(simulation=names(all_model_samples), AWV=NA, CNY=NA, KNC=NA, 
                          PAC=NA, SNC=NA, CWV=NA, BTN=NA)
for (i in 1:length(all_model_samples)){
  param_lower$AWV[i] <- quantile(all_model_samples[[i]]$AWV, probs = c(2.5)/100)
  param_lower$CNY[i] <- quantile(all_model_samples[[i]]$CNY, probs = c(2.5)/100)
  param_lower$KNC[i] <- quantile(all_model_samples[[i]]$KNC, probs = c(2.5)/100)
  param_lower$PAC[i] <- quantile(all_model_samples[[i]]$PAC, probs = c(2.5)/100)
  param_lower$SNC[i] <- quantile(all_model_samples[[i]]$SNC, probs = c(2.5)/100)
  param_lower$CWV[i] <- quantile(all_model_samples[[i]]$CWV, probs = c(2.5)/100)
  param_lower$BTN[i] <- quantile(all_model_samples[[i]]$BTN, probs = c(2.5)/100)
}
param_lower


# example of summarizing lower CI
param_upper <- data.frame(simulation=names(all_model_samples), AWV=NA, CNY=NA, KNC=NA, 
                          PAC=NA, SNC=NA, CWV=NA, BTN=NA)
for (i in 1:length(all_model_samples)){
  param_upper$AWV[i] <- quantile(all_model_samples[[i]]$AWV, probs = c(97.5)/100)
  param_upper$CNY[i] <- quantile(all_model_samples[[i]]$CNY, probs = c(97.5)/100)
  param_upper$KNC[i] <- quantile(all_model_samples[[i]]$KNC, probs = c(97.5)/100)
  param_upper$PAC[i] <- quantile(all_model_samples[[i]]$PAC, probs = c(97.5)/100)
  param_upper$SNC[i] <- quantile(all_model_samples[[i]]$SNC, probs = c(97.5)/100)
  param_upper$CWV[i] <- quantile(all_model_samples[[i]]$CWV, probs = c(97.5)/100)
  param_upper$BTN[i] <- quantile(all_model_samples[[i]]$BTN, probs = c(97.5)/100)
}
param_upper

################################################################################
###                       ASSEMBLE SIMULATION ESTIMATES                      ###
################################################################################

## West Virginia (NE)
AWV_simulations <- subset(param_delta, select=simulation)
AWV_means <- subset(param_delta, select=AWV)
AWV_lower <- subset(param_lower, select=AWV)
AWV_upper <- subset(param_upper, select=AWV)

AWV_overall <- cbind(AWV_simulations, AWV_means, AWV_lower, AWV_upper)
colnames(AWV_overall) <- c("Simulation", "Mean", "Lower", "Upper")
AWV_overall$CI_range <- AWV_overall$Upper - AWV_overall$Lower
AWV_overall <- AWV_overall[order(AWV_overall$CI_range), ]


selected_row <- subset(df_app_sex, App.Code == "AWV")
AWV_true <- data.frame(Simulation = "True Model Output",
                       Mean = selected_row$App.Delta,
                       Lower = selected_row$App.Lower,
                       Upper = selected_row$App.Upper,
                       CI_range = selected_row$App.Upper - selected_row$App.Lower)
AWV_overall <- rbind(AWV_overall, AWV_true)

mycolours <- c(values = ifelse(AWV_overall$Simulation=="True Model Output", "red", "black"))

AWV_overall <- AWV_overall %>% mutate(test = ifelse(Simulation == "True Model Output", 1, 0)) %>%
  mutate(test = as.factor(test))


AWV_overall$Occasions <- substr(AWV_overall$Simulation, start = 9, stop = 9)
AWV_overall$Sex_Ratio <- substr(AWV_overall$Simulation, start = 13, stop = 15)
AWV_overall$Sample_Size <- substr(AWV_overall$Simulation, start = 18, stop = 19)
AWV_overall$Recapture_Prob <- substr(AWV_overall$Simulation, start = 23, stop = 25)


################################################################################
###                               CREATE VISUAL                              ###
################################################################################

## You'll need to combine these two figures to replicate the figure in manuscript.
## This figure has the recap probability included on the side. For the figure you'll
## have to manually indicate which column represents what metadata. 
sim_wv <- ggplot(data = AWV_overall, aes(x = Mean, y = reorder(Simulation, -CI_range), color = test))+ 
  geom_point() +
  xlim(-0.6, 1.3) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper, color = test), width = 0.2) +
  scale_color_manual(values = c("black", "red"))+
  theme(panel.grid.major  = element_line(color = "white"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"),
        legend.position = "none")+
  ggtitle("West Virginia (NE)")+
  xlab("Mean Difference Between Sexes")+
  geom_text(aes(x = max(Mean) + 0.7, y = reorder(Simulation, -CI_range), label = Occasions),
            hjust = 0, size = 3)+
  geom_text(aes(x = max(Mean) + 0.85, y = reorder(Simulation, -CI_range), label = Sex_Ratio),
            hjust = 0, size = 3)+
  geom_text(aes(x = max(Mean) + 1.00, y = reorder(Simulation, -CI_range), label = Sample_Size),
            hjust = 0, size = 3)+
  geom_text(aes(x = max(Mean) + 1.15, y = reorder(Simulation, -CI_range), label = Recapture_Prob),
            hjust = 0, size = 3)

## Save part 1 output
setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Tables and Figures/Figure A7")
ggsave("sim_wv.png", dpi = 1200, width = 5, height = 9.5, units = "in")

## This figure has the CI range included on the side
sim_wv_2 <- ggplot(data = AWV_overall, aes(x = Mean, y = reorder(Simulation, -CI_range), color = test))+ 
  geom_point() +
  xlim(-0.6, 1.3) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper, color = test), width = 0.2) +
  scale_color_manual(values = c("black", "red"))+
  theme(panel.grid.major  = element_line(color = "white"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"),
        legend.position = "none")+
  ggtitle("West Virginia (NE)")+
  xlab("Mean Difference Between Sexes")+
  geom_text(aes(x = max(Mean) + 0.7, y = reorder(Simulation, -CI_range), label = Occasions),
            hjust = 0, size = 3)+
  geom_text(aes(x = max(Mean) + 0.85, y = reorder(Simulation, -CI_range), label = Sex_Ratio),
            hjust = 0, size = 3)+
  geom_text(aes(x = max(Mean) + 1.00, y = reorder(Simulation, -CI_range), label = Sample_Size),
            hjust = 0, size = 3)+
  geom_text(aes(x = max(Mean) + 1.15, y = reorder(Simulation, -CI_range), label = sprintf("%.3f", CI_range)),
            hjust = 0, size = 3)

## Save part 2 output
setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Tables and Figures/Figure A7")
ggsave("sim_wv_2.png", dpi = 1200, width = 5, height = 9.5, units = "in")

