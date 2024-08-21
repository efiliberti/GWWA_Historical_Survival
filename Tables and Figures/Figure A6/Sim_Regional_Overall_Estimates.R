## Figure A6 in manuscript. This script shows the mean overall apparent annual
## survival (and 95% credible intervals) estimates for the Appalachian region.
## Taken from mean posterior estimates of 30 simulations for 54 different
## respective scenarios. 

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
  mat_param1 <- matrix(NA, nrow=nsim, ncol=nmcmc) # matrices to store - 1 row per model 1 col per posterior sample
  for (z in 1:nsim){ 
    load(files_sub[z]) # read in each model
    mods[[z]] <- simulAppMultiSite
    # fill parameter matrices with posterior samples for each model
    mat_param1[z,] <- mods[[z]]$sims.list$mean.overall.survival.App
  }
  samples <- list(param1=mat_param1) 
  all_model_samples[[q]] <- samples
}

names(all_model_samples) <- simnames # add names to elements of the list 



## 95% quanile of posterior means - comes from the posterior distribution - mean of means

# example of summarizing means (using all samples)
param_means <- data.frame(simulation=names(all_model_samples), param1=NA)
for (i in 1:length(all_model_samples)){
  param_means$param1[i] <- mean(all_model_samples[[i]]$param1)
}
param_means

# example of summarizing lower CI
param_lower <- data.frame(simulation=names(all_model_samples), param1=NA)
for (i in 1:length(all_model_samples)){
  param_lower$param1[i] <- quantile(all_model_samples[[i]]$param1, probs = c(2.5)/100)
}
param_lower

# example of summarizing lower CI
param_upper <- data.frame(simulation=names(all_model_samples), param1=NA)
for (i in 1:length(all_model_samples)){
  param_upper$param1[i] <- quantile(all_model_samples[[i]]$param1, probs = c(97.5)/100)
}
param_upper


################################################################################
###                       ASSEMBLE SIMULATION ESTIMATES                      ###
################################################################################

App_simulations <- subset(param_means, select=simulation)
App_means <- subset(param_means, select=param1)
App_lower <- subset(param_lower, select=param1)
App_upper <- subset(param_upper, select=param1)

App_overall <- cbind(App_simulations, App_means, App_lower, App_upper)
colnames(App_overall) <- c("Simulation", "Mean", "Lower", "Upper")
App_overall$CI_range <- round(App_overall$Upper - App_overall$Lower, digits = 3)
App_overall$TrueEst <- 0.495

App_true <- data.frame(Simulation = "True Model Output", Mean = 0.495,
                       Lower = 0.366, Upper = 0.619, CI_range = 0.253,
                       TrueEst = 0.495)

App_overall <- rbind(App_overall, App_true)

App_overall$Difference <- sprintf("%.3f", round(App_overall$Mean - App_overall$TrueEst))
App_overall <- App_overall[order(App_overall$CI_range), ]

App_overall$Occasions <- substr(App_overall$Simulation, start = 9, stop = 9)
App_overall$Sex_Ratio <- substr(App_overall$Simulation, start = 13, stop = 15)
App_overall$Sample_Size <- substr(App_overall$Simulation, start = 18, stop = 19)
App_overall$Recapture_Prob <- substr(App_overall$Simulation, start = 23, stop = 25)

mycolours <- c(values = ifelse(App_overall$Simulation=="True Model Output", "red", "black"))

App_overall <- App_overall %>% mutate(test = ifelse(Simulation == "True Model Output", 1, 0)) %>%
  mutate(test = as.factor(test))


################################################################################
###                               CREATE VISUAL                              ###
################################################################################

## You'll need to combine these two figures to replicate the figure in manuscript.
## This figure has the CI-range included on the side. For the figure you'll
## have to manually indicate which column represents what metadata. 
sim_app <- ggplot(data = App_overall, aes(x = Mean, y = reorder(Simulation, -CI_range), color = test))+ 
  geom_point() +
  xlim(0, 1.55) +
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
  ggtitle("Appalachian Simulation Output")+
  xlab("Mean Overall Regional Survival")+
  geom_text(aes(x = max(Mean) + 0.35, y = reorder(Simulation, -CI_range), label = Occasions),
            hjust = 0, size = 3)+
  geom_text(aes(x = max(Mean) + 0.45, y = reorder(Simulation, -CI_range), label = Sex_Ratio),
            hjust = 0, size = 3)+
  geom_text(aes(x = max(Mean) + 0.60, y = reorder(Simulation, -CI_range), label = Sample_Size),
            hjust = 0, size = 3)+
  geom_text(aes(x = max(Mean) + 0.73, y = reorder(Simulation, -CI_range), label = sprintf("%.3f", CI_range)),
            hjust = 0, size = 3)

## Save part 1 output
setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Tables and Figures/Figure A6")
ggsave("simapp.png", dpi = 1200, width = 5, height = 9.5, units = "in")

## This figure has the recapture probability included on the side
sim_app <- ggplot(data = App_overall, aes(x = Mean, y = reorder(Simulation, -CI_range), color = test))+ 
  geom_point() +
  xlim(0, 1.55) +
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
  ggtitle("Appalachian Simulation Output")+
  xlab("Mean Overall Regional Survival")+
  geom_text(aes(x = max(Mean) + 0.35, y = reorder(Simulation, -CI_range), label = Occasions),
            hjust = 0, size = 3)+
  geom_text(aes(x = max(Mean) + 0.45, y = reorder(Simulation, -CI_range), label = Sex_Ratio),
            hjust = 0, size = 3)+
  geom_text(aes(x = max(Mean) + 0.60, y = reorder(Simulation, -CI_range), label = Sample_Size),
            hjust = 0, size = 3)+
  geom_text(aes(x = max(Mean) + 0.73, y = reorder(Simulation, -CI_range), label = Recapture_Prob),
            hjust = 0, size = 3)

## Save part 2 output
setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Tables and Figures/Figure A6")
ggsave("sim_app.png", dpi = 1200, width = 5, height = 9.5, units = "in")
