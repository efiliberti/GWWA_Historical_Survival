## Figure 5 in manuscript. This visualization shows the mean difference in site
## specific annual survival (with 95% credible intervals) for both the Great
## Lakes and the Appalachians. 

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
library(ggpubr)

## Indicate desired simulations
nsim <- 30 #s
n.occ <- c(3, 6, 9) #o  
sex.rat <- c(0.3, 0.5) #f
samp.siz <- c(10, 30, 50) #n
recap <- c(0.3, 0.6, 0.9) #p

## Read model output in from file where simulations are stored
setwd('/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Simulations/App Simulation Rnorm')

## Make a list of full file names
filenames <- list.files('.', full.names=T)

## Take the simulation number and ./ off of each file name
simnames <- unique(substr(filenames, start=3, stop=27))

## Indicate number of expected posteriors
nmcmc <- 1800

## Identify number of combinations of parameters (54)
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
  
###################################  AWV  ###################################### 
AWV_simulations <- subset(param_delta, select=simulation)
AWV_means <- subset(param_delta, select=AWV)
AWV_lower <- subset(param_lower, select=AWV)
AWV_upper <- subset(param_upper, select=AWV)
  
AWV_overall <- cbind(AWV_simulations, AWV_means, AWV_lower, AWV_upper)
colnames(AWV_overall) <- c("Simulation", "Mean", "Lower", "Upper")
AWV_overall$CI_range <- AWV_overall$Upper - AWV_overall$Lower
AWV_overall <- AWV_overall[order(AWV_overall$CI_range), ]

  
## Subset dataframe 9.3.5.9
AWV_9359 <- AWV_overall[grep("nocc9", AWV_overall$Simulation),]
AWV_9359 <- AWV_9359[grep("sr0.3", AWV_9359$Simulation),]
AWV_9359 <- AWV_9359[grep("n50", AWV_9359$Simulation),]
AWV_9359 <- AWV_9359[grep("rc0.9", AWV_9359$Simulation),]
  
## Subset dataframe 9.5.5.9
AWV_9559 <- AWV_overall[grep("nocc9", AWV_overall$Simulation),]
AWV_9559 <- AWV_9559[grep("sr0.5", AWV_9559$Simulation),]
AWV_9559 <- AWV_9559[grep("n50", AWV_9559$Simulation),]
AWV_9559 <- AWV_9559[grep("rc0.9", AWV_9559$Simulation),]
  
## Subset dataframe 3.4.5.9
AWV_3359 <- AWV_overall[grep("nocc3", AWV_overall$Simulation),]
AWV_3359 <- AWV_3359[grep("sr0.3", AWV_3359$Simulation),]
AWV_3359 <- AWV_3359[grep("n50", AWV_3359$Simulation),]
AWV_3359 <- AWV_3359[grep("rc0.9", AWV_3359$Simulation),]
  
selected_row <- subset(df_app_sex, App.Code == "AWV")
AWV_true <- data.frame(Simulation = "True Model Output",
                         Mean = selected_row$App.Delta,
                         Lower = selected_row$App.Lower,
                         Upper = selected_row$App.Upper,
                         CI_range = selected_row$App.Upper - selected_row$App.Lower)
AWV_overall <- rbind(AWV_9359, AWV_3359, AWV_9559, AWV_true)
  
mycolours <- c(values = ifelse(AWV_overall$Simulation=="True Model Output", "red", "black"))
  
AWV_overall <- AWV_overall %>% mutate(test = ifelse(Simulation == "True Model Output", 1, 0)) %>%
    mutate(test = as.factor(test))
  
AWV_overall[AWV_overall$Simulation == "sim_nocc9_sr0.3_n50_rc0.9", "Simulation"] <- "o=9, sr=30%, n=50, p=90%"
AWV_overall[AWV_overall$Simulation == "sim_nocc9_sr0.5_n50_rc0.9", "Simulation"] <- "o=9, sr=50%, n=50, p=90%"
AWV_overall[AWV_overall$Simulation == "sim_nocc3_sr0.3_n50_rc0.9", "Simulation"] <- "o=3, sr=50%, n=50, p=90%"
  
sim_a <- ggplot(data = AWV_overall, aes(x = Mean, y = reorder(Simulation, -CI_range), color = test))+
    geom_vline(aes(xintercept=0), col = "gray", linetype = 'dashed', size = 0.75)+
    geom_point()+
    xlim(-0.5,0.5)+
    geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.05)+
    ylab("Scenario") +
    scale_color_manual(values = c("black", "red"))+
    theme(panel.grid.major  = element_line(color = "white"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(size = 10, hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          legend.title=element_text(colour="white"),
          legend.position = "none")+
    ggtitle("West Virginia (NE)")+
    theme(plot.title = element_text(size = 11, face = "bold"))+
  xlab("Mean Difference Between Sexes")
  
###################################  CNY  ###################################### 
CNY_simulations <- subset(param_delta, select=simulation)
CNY_means <- subset(param_delta, select=CNY)
CNY_lower <- subset(param_lower, select=CNY)
CNY_upper <- subset(param_upper, select=CNY)
  
CNY_overall <- cbind(CNY_simulations, CNY_means, CNY_lower, CNY_upper)
colnames(CNY_overall) <- c("Simulation", "Mean", "Lower", "Upper")
CNY_overall$CI_range <- CNY_overall$Upper - CNY_overall$Lower
CNY_overall <- CNY_overall[order(CNY_overall$CI_range), ]
  
  
## Subset dataframe 9.3.5.9
CNY_9359 <- CNY_overall[grep("nocc9", CNY_overall$Simulation),]
CNY_9359 <- CNY_9359[grep("sr0.3", CNY_9359$Simulation),]
CNY_9359 <- CNY_9359[grep("n50", CNY_9359$Simulation),]
CNY_9359 <- CNY_9359[grep("rc0.9", CNY_9359$Simulation),]
  
## Subset dataframe 9.5.5.9
CNY_9559 <- CNY_overall[grep("nocc9", CNY_overall$Simulation),]
CNY_9559 <- CNY_9559[grep("sr0.5", CNY_9559$Simulation),]
CNY_9559 <- CNY_9559[grep("n50", CNY_9559$Simulation),]
CNY_9559 <- CNY_9559[grep("rc0.9", CNY_9559$Simulation),]
  
## Subset dataframe 3.4.5.9
CNY_3359 <- CNY_overall[grep("nocc3", CNY_overall$Simulation),]
CNY_3359 <- CNY_3359[grep("sr0.3", CNY_3359$Simulation),]
CNY_3359 <- CNY_3359[grep("n50", CNY_3359$Simulation),]
CNY_3359 <- CNY_3359[grep("rc0.9", CNY_3359$Simulation),]
  
selected_row <- subset(df_app_sex, App.Code == "CNY")
CNY_true <- data.frame(Simulation = "True Model Output",
                         Mean = selected_row$App.Delta,
                         Lower = selected_row$App.Lower,
                         Upper = selected_row$App.Upper,
                         CI_range = selected_row$App.Upper - selected_row$App.Lower)
CNY_overall <- rbind(CNY_9359, CNY_3359, CNY_9559, CNY_true)
  
mycolours <- c(values = ifelse(CNY_overall$Simulation=="True Model Output", "red", "black"))
  
CNY_overall <- CNY_overall %>% mutate(test = ifelse(Simulation == "True Model Output", 1, 0)) %>%
    mutate(test = as.factor(test))
  
CNY_overall[CNY_overall$Simulation == "sim_nocc9_sr0.3_n50_rc0.9", "Simulation"] <- "o=9, sr=30%, n=50, p=90%"
CNY_overall[CNY_overall$Simulation == "sim_nocc9_sr0.5_n50_rc0.9", "Simulation"] <- "o=9, sr=50%, n=50, p=90%"
CNY_overall[CNY_overall$Simulation == "sim_nocc3_sr0.3_n50_rc0.9", "Simulation"] <- "o=3, sr=50%, n=50, p=90%"
  
sim_b <- ggplot(data = CNY_overall, aes(x = Mean, y = reorder(Simulation, -CI_range), color = test))+
    geom_vline(aes(xintercept=0), col = "gray", linetype = 'dashed', size = 0.75)+
    geom_point()+
    xlim(-0.5,0.5)+
    geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.05)+
    ylab("Scenario") +
    scale_color_manual(values = c("black", "red"))+
    theme(panel.grid.major  = element_line(color = "white"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(size = 10, hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          legend.title=element_text(colour="white"),
          legend.position = "none")+
    ggtitle("New York (S)")+
    theme(plot.title = element_text(size = 11, face = "bold"))+
    xlab("Mean Difference Between Sexes")
  
###################################  KNC  ###################################### 
KNC_simulations <- subset(param_delta, select=simulation)
KNC_means <- subset(param_delta, select=KNC)
KNC_lower <- subset(param_lower, select=KNC)
KNC_upper <- subset(param_upper, select=KNC)
  
KNC_overall <- cbind(KNC_simulations, KNC_means, KNC_lower, KNC_upper)
colnames(KNC_overall) <- c("Simulation", "Mean", "Lower", "Upper")
KNC_overall$CI_range <- KNC_overall$Upper - KNC_overall$Lower
KNC_overall <- KNC_overall[order(KNC_overall$CI_range), ]
  
  
## Subset dataframe 9.3.5.9
KNC_9359 <- KNC_overall[grep("nocc9", KNC_overall$Simulation),]
KNC_9359 <- KNC_9359[grep("sr0.3", KNC_9359$Simulation),]
KNC_9359 <- KNC_9359[grep("n50", KNC_9359$Simulation),]
KNC_9359 <- KNC_9359[grep("rc0.9", KNC_9359$Simulation),]
  
## Subset dataframe 9.5.5.9
KNC_9559 <- KNC_overall[grep("nocc9", KNC_overall$Simulation),]
KNC_9559 <- KNC_9559[grep("sr0.5", KNC_9559$Simulation),]
KNC_9559 <- KNC_9559[grep("n50", KNC_9559$Simulation),]
KNC_9559 <- KNC_9559[grep("rc0.9", KNC_9559$Simulation),]
  
## Subset dataframe 3.4.5.9
KNC_3359 <- KNC_overall[grep("nocc3", KNC_overall$Simulation),]
KNC_3359 <- KNC_3359[grep("sr0.3", KNC_3359$Simulation),]
KNC_3359 <- KNC_3359[grep("n50", KNC_3359$Simulation),]
KNC_3359 <- KNC_3359[grep("rc0.9", KNC_3359$Simulation),]
  
selected_row <- subset(df_app_sex, App.Code == "KNC")
KNC_true <- data.frame(Simulation = "True Model Output",
                         Mean = selected_row$App.Delta,
                         Lower = selected_row$App.Lower,
                         Upper = selected_row$App.Upper,
                         CI_range = selected_row$App.Upper - selected_row$App.Lower)
KNC_overall <- rbind(KNC_9359, KNC_3359, KNC_9559, KNC_true)
  
mycolours <- c(values = ifelse(KNC_overall$Simulation=="True Model Output", "red", "black"))
  
KNC_overall <- KNC_overall %>% mutate(test = ifelse(Simulation == "True Model Output", 1, 0)) %>%
    mutate(test = as.factor(test))
  
KNC_overall[KNC_overall$Simulation == "sim_nocc9_sr0.3_n50_rc0.9", "Simulation"] <- "o=9, sr=30%, n=50, p=90%"
KNC_overall[KNC_overall$Simulation == "sim_nocc9_sr0.5_n50_rc0.9", "Simulation"] <- "o=9, sr=50%, n=50, p=90%"
KNC_overall[KNC_overall$Simulation == "sim_nocc3_sr0.3_n50_rc0.9", "Simulation"] <- "o=3, sr=50%, n=50, p=90%"
  
sim_c <- ggplot(data = KNC_overall, aes(x = Mean, y = reorder(Simulation, -CI_range), color = test))+
    geom_vline(aes(xintercept=0), col = "gray", linetype = 'dashed', size = 0.75)+
    geom_point()+
    xlim(-0.5,0.5)+
    geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.05)+
    ylab("Scenario") +
    scale_color_manual(values = c("black", "red"))+
    theme(panel.grid.major  = element_line(color = "white"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(size = 10, hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          legend.title=element_text(colour="white"),
          legend.position = "none")+
    ggtitle("North Carolina (W)")+
    theme(plot.title = element_text(size = 11, face = "bold"))+
    xlab("Mean Difference Between Sexes")
  
###################################  PAC  ###################################### 
PAC_simulations <- subset(param_delta, select=simulation)
PAC_means <- subset(param_delta, select=PAC)
PAC_lower <- subset(param_lower, select=PAC)
PAC_upper <- subset(param_upper, select=PAC)
  
PAC_overall <- cbind(PAC_simulations, PAC_means, PAC_lower, PAC_upper)
colnames(PAC_overall) <- c("Simulation", "Mean", "Lower", "Upper")
PAC_overall$CI_range <- PAC_overall$Upper - PAC_overall$Lower
PAC_overall <- PAC_overall[order(PAC_overall$CI_range), ]
  
  
## Subset dataframe 9.3.5.9
PAC_9359 <- PAC_overall[grep("nocc9", PAC_overall$Simulation),]
PAC_9359 <- PAC_9359[grep("sr0.3", PAC_9359$Simulation),]
PAC_9359 <- PAC_9359[grep("n50", PAC_9359$Simulation),]
PAC_9359 <- PAC_9359[grep("rc0.9", PAC_9359$Simulation),]
  
## Subset dataframe 9.5.5.9
PAC_9559 <- PAC_overall[grep("nocc9", PAC_overall$Simulation),]
PAC_9559 <- PAC_9559[grep("sr0.5", PAC_9559$Simulation),]
PAC_9559 <- PAC_9559[grep("n50", PAC_9559$Simulation),]
PAC_9559 <- PAC_9559[grep("rc0.9", PAC_9559$Simulation),]
  
## Subset dataframe 3.4.5.9
PAC_3359 <- PAC_overall[grep("nocc3", PAC_overall$Simulation),]
PAC_3359 <- PAC_3359[grep("sr0.3", PAC_3359$Simulation),]
PAC_3359 <- PAC_3359[grep("n50", PAC_3359$Simulation),]
PAC_3359 <- PAC_3359[grep("rc0.9", PAC_3359$Simulation),]
  
selected_row <- subset(df_app_sex, App.Code == "PAC")
PAC_true <- data.frame(Simulation = "True Model Output",
                         Mean = selected_row$App.Delta,
                         Lower = selected_row$App.Lower,
                         Upper = selected_row$App.Upper,
                         CI_range = selected_row$App.Upper - selected_row$App.Lower)
PAC_overall <- rbind(PAC_9359, PAC_3359, PAC_9559, PAC_true)
  
mycolours <- c(values = ifelse(PAC_overall$Simulation=="True Model Output", "red", "black"))
  
PAC_overall <- PAC_overall %>% mutate(test = ifelse(Simulation == "True Model Output", 1, 0)) %>%
    mutate(test = as.factor(test))
  
PAC_overall[PAC_overall$Simulation == "sim_nocc9_sr0.3_n50_rc0.9", "Simulation"] <- "o=9, sr=30%, n=50, p=90%"
PAC_overall[PAC_overall$Simulation == "sim_nocc9_sr0.5_n50_rc0.9", "Simulation"] <- "o=9, sr=50%, n=50, p=90%"
PAC_overall[PAC_overall$Simulation == "sim_nocc3_sr0.3_n50_rc0.9", "Simulation"] <- "o=3, sr=50%, n=50, p=90%"
  
sim_d <- ggplot(data = PAC_overall, aes(x = Mean, y = reorder(Simulation, -CI_range), color = test))+
    geom_vline(aes(xintercept=0), col = "gray", linetype = 'dashed', size = 0.75)+
    geom_point()+
    xlim(-0.5,0.5)+
    geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.05)+
    ylab("Scenario") +
    scale_color_manual(values = c("black", "red"))+
    theme(panel.grid.major  = element_line(color = "white"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(size = 10, hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          legend.title=element_text(colour="white"),
          legend.position = "none")+
    ggtitle("Pennsylvania (C)")+
    theme(plot.title = element_text(size = 11, face = "bold"))+
    xlab("Mean Difference Between Sexes")
  
  
  
###################################  SNC  ###################################### 
SNC_simulations <- subset(param_delta, select=simulation)
SNC_means <- subset(param_delta, select=SNC)
SNC_lower <- subset(param_lower, select=SNC)
SNC_upper <- subset(param_upper, select=SNC)
  
SNC_overall <- cbind(SNC_simulations, SNC_means, SNC_lower, SNC_upper)
colnames(SNC_overall) <- c("Simulation", "Mean", "Lower", "Upper")
SNC_overall$CI_range <- SNC_overall$Upper - SNC_overall$Lower
SNC_overall <- SNC_overall[order(SNC_overall$CI_range), ]
  
  
## Subset dataframe 9.3.5.9
SNC_9359 <- SNC_overall[grep("nocc9", SNC_overall$Simulation),]
SNC_9359 <- SNC_9359[grep("sr0.3", SNC_9359$Simulation),]
SNC_9359 <- SNC_9359[grep("n50", SNC_9359$Simulation),]
SNC_9359 <- SNC_9359[grep("rc0.9", SNC_9359$Simulation),]
  
## Subset dataframe 9.5.5.9
SNC_9559 <- SNC_overall[grep("nocc9", SNC_overall$Simulation),]
SNC_9559 <- SNC_9559[grep("sr0.5", SNC_9559$Simulation),]
SNC_9559 <- SNC_9559[grep("n50", SNC_9559$Simulation),]
SNC_9559 <- SNC_9559[grep("rc0.9", SNC_9559$Simulation),]
  
## Subset dataframe 3.4.5.9
SNC_3359 <- SNC_overall[grep("nocc3", SNC_overall$Simulation),]
SNC_3359 <- SNC_3359[grep("sr0.3", SNC_3359$Simulation),]
SNC_3359 <- SNC_3359[grep("n50", SNC_3359$Simulation),]
SNC_3359 <- SNC_3359[grep("rc0.9", SNC_3359$Simulation),]
  
selected_row <- subset(df_app_sex, App.Code == "SNC")
SNC_true <- data.frame(Simulation = "True Model Output",
                         Mean = selected_row$App.Delta,
                         Lower = selected_row$App.Lower,
                         Upper = selected_row$App.Upper,
                         CI_range = selected_row$App.Upper - selected_row$App.Lower)
SNC_overall <- rbind(SNC_9359, SNC_3359, SNC_9559, SNC_true)
  
mycolours <- c(values = ifelse(SNC_overall$Simulation=="True Model Output", "red", "black"))
  
SNC_overall <- SNC_overall %>% mutate(test = ifelse(Simulation == "True Model Output", 1, 0)) %>%
    mutate(test = as.factor(test))
  
SNC_overall[SNC_overall$Simulation == "sim_nocc9_sr0.3_n50_rc0.9", "Simulation"] <- "o=9, sr=30%, n=50, p=90%"
SNC_overall[SNC_overall$Simulation == "sim_nocc9_sr0.5_n50_rc0.9", "Simulation"] <- "o=9, sr=50%, n=50, p=90%"
SNC_overall[SNC_overall$Simulation == "sim_nocc3_sr0.3_n50_rc0.9", "Simulation"] <- "o=3, sr=50%, n=50, p=90%"
  
sim_e <- ggplot(data = SNC_overall, aes(x = Mean, y = reorder(Simulation, -CI_range), color = test))+
    geom_vline(aes(xintercept=0), col = "gray", linetype = 'dashed', size = 0.75)+
    geom_point()+
    xlim(-0.5,0.5)+
    geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.05)+
    ylab("Scenario") +
    scale_color_manual(values = c("black", "red"))+
    theme(panel.grid.major  = element_line(color = "white"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(size = 10, hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          legend.title=element_text(colour="white"),
          legend.position = "none")+
    ggtitle("North Carolina (NW)")+
    theme(plot.title = element_text(size = 11, face = "bold"))+
    xlab("Mean Difference Between Sexes")
  
  
###################################  CWV  ###################################### 
CWV_simulations <- subset(param_delta, select=simulation)
CWV_means <- subset(param_delta, select=CWV)
CWV_lower <- subset(param_lower, select=CWV)
CWV_upper <- subset(param_upper, select=CWV)
  
CWV_overall <- cbind(CWV_simulations, CWV_means, CWV_lower, CWV_upper)
colnames(CWV_overall) <- c("Simulation", "Mean", "Lower", "Upper")
CWV_overall$CI_range <- CWV_overall$Upper - CWV_overall$Lower
CWV_overall <- CWV_overall[order(CWV_overall$CI_range), ]
  
  
## Subset dataframe 9.3.5.9
CWV_9359 <- CWV_overall[grep("nocc9", CWV_overall$Simulation),]
CWV_9359 <- CWV_9359[grep("sr0.3", CWV_9359$Simulation),]
CWV_9359 <- CWV_9359[grep("n50", CWV_9359$Simulation),]
CWV_9359 <- CWV_9359[grep("rc0.9", CWV_9359$Simulation),]
  
## Subset dataframe 9.5.5.9
CWV_9559 <- CWV_overall[grep("nocc9", CWV_overall$Simulation),]
CWV_9559 <- CWV_9559[grep("sr0.5", CWV_9559$Simulation),]
CWV_9559 <- CWV_9559[grep("n50", CWV_9559$Simulation),]
CWV_9559 <- CWV_9559[grep("rc0.9", CWV_9559$Simulation),]
  
## Subset dataframe 3.4.5.9
CWV_3359 <- CWV_overall[grep("nocc3", CWV_overall$Simulation),]
CWV_3359 <- CWV_3359[grep("sr0.3", CWV_3359$Simulation),]
CWV_3359 <- CWV_3359[grep("n50", CWV_3359$Simulation),]
CWV_3359 <- CWV_3359[grep("rc0.9", CWV_3359$Simulation),]
  
selected_row <- subset(df_app_sex, App.Code == "CWV")
CWV_true <- data.frame(Simulation = "True Model Output",
                         Mean = selected_row$App.Delta,
                         Lower = selected_row$App.Lower,
                         Upper = selected_row$App.Upper,
                         CI_range = selected_row$App.Upper - selected_row$App.Lower)
CWV_overall <- rbind(CWV_9359, CWV_3359, CWV_9559, CWV_true)
  
mycolours <- c(values = ifelse(CWV_overall$Simulation=="True Model Output", "red", "black"))
  
CWV_overall <- CWV_overall %>% mutate(test = ifelse(Simulation == "True Model Output", 1, 0)) %>%
    mutate(test = as.factor(test))
  
CWV_overall[CWV_overall$Simulation == "sim_nocc9_sr0.3_n50_rc0.9", "Simulation"] <- "o=9, sr=30%, n=50, p=90%"
CWV_overall[CWV_overall$Simulation == "sim_nocc9_sr0.5_n50_rc0.9", "Simulation"] <- "o=9, sr=50%, n=50, p=90%"
CWV_overall[CWV_overall$Simulation == "sim_nocc3_sr0.3_n50_rc0.9", "Simulation"] <- "o=3, sr=50%, n=50, p=90%"
  
sim_f <- ggplot(data = CWV_overall, aes(x = Mean, y = reorder(Simulation, -CI_range), color = test))+
    geom_vline(aes(xintercept=0), col = "gray", linetype = 'dashed', size = 0.75)+
    geom_point()+
    xlim(-0.5,0.5)+
    geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.05)+
    ylab("Scenario") +
    scale_color_manual(values = c("black", "red"))+
    theme(panel.grid.major  = element_line(color = "white"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(size = 10, hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          legend.title=element_text(colour="white"),
          legend.position = "none")+
    ggtitle("West Virginia (S)")+
    theme(plot.title = element_text(size = 11, face = "bold"))+
    xlab("Mean Difference Between Sexes")
  
###################################  BTN  ###################################### 
BTN_simulations <- subset(param_delta, select=simulation)
BTN_means <- subset(param_delta, select=BTN)
BTN_lower <- subset(param_lower, select=BTN)
BTN_upper <- subset(param_upper, select=BTN)
  
BTN_overall <- cbind(BTN_simulations, BTN_means, BTN_lower, BTN_upper)
colnames(BTN_overall) <- c("Simulation", "Mean", "Lower", "Upper")
BTN_overall$CI_range <- BTN_overall$Upper - BTN_overall$Lower
BTN_overall <- BTN_overall[order(BTN_overall$CI_range), ]
  
  
## Subset dataframe 9.3.5.9
BTN_9359 <- BTN_overall[grep("nocc9", BTN_overall$Simulation),]
BTN_9359 <- BTN_9359[grep("sr0.3", BTN_9359$Simulation),]
BTN_9359 <- BTN_9359[grep("n50", BTN_9359$Simulation),]
BTN_9359 <- BTN_9359[grep("rc0.9", BTN_9359$Simulation),]
  
## Subset dataframe 9.5.5.9
BTN_9559 <- BTN_overall[grep("nocc9", BTN_overall$Simulation),]  
BTN_9559 <- BTN_9559[grep("sr0.5", BTN_9559$Simulation),]
BTN_9559 <- BTN_9559[grep("n50", BTN_9559$Simulation),]
BTN_9559 <- BTN_9559[grep("rc0.9", BTN_9559$Simulation),]
  
## Subset dataframe 3.4.5.9
BTN_3359 <- BTN_overall[grep("nocc3", BTN_overall$Simulation),]
BTN_3359 <- BTN_3359[grep("sr0.3", BTN_3359$Simulation),]
BTN_3359 <- BTN_3359[grep("n50", BTN_3359$Simulation),]
BTN_3359 <- BTN_3359[grep("rc0.9", BTN_3359$Simulation),]
  
selected_row <- subset(df_app_sex, App.Code == "BTN")
BTN_true <- data.frame(Simulation = "True Model Output",
                         Mean = selected_row$App.Delta,
                         Lower = selected_row$App.Lower,
                         Upper = selected_row$App.Upper,
                         CI_range = selected_row$App.Upper - selected_row$App.Lower)
BTN_overall <- rbind(BTN_9359, BTN_3359, BTN_9559, BTN_true)
  
mycolours <- c(values = ifelse(BTN_overall$Simulation=="True Model Output", "red", "black"))
  
BTN_overall <- BTN_overall %>% mutate(test = ifelse(Simulation == "True Model Output", 1, 0)) %>%
    mutate(test = as.factor(test))
  
BTN_overall[BTN_overall$Simulation == "sim_nocc9_sr0.3_n50_rc0.9", "Simulation"] <- "o=9, sr=30%, n=50, p=90%"
BTN_overall[BTN_overall$Simulation == "sim_nocc9_sr0.5_n50_rc0.9", "Simulation"] <- "o=9, sr=50%, n=50, p=90%"  
BTN_overall[BTN_overall$Simulation == "sim_nocc3_sr0.3_n50_rc0.9", "Simulation"] <- "o=3, sr=50%, n=50, p=90%"
  
sim_g <- ggplot(data = BTN_overall, aes(x = Mean, y = reorder(Simulation, -CI_range), color = test))+
    geom_vline(aes(xintercept=0), col = "gray", linetype = 'dashed', size = 0.75)+
    geom_point()+
    xlim(-0.5,0.5)+
    geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.05)+
    ylab("Scenario") +
    scale_color_manual(values = c("black", "red"))+
    theme(panel.grid.major  = element_line(color = "white"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(size = 10, hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          legend.title=element_text(colour="white"),
          legend.position = "none")+
    ggtitle("Tennessee")+
    theme(plot.title = element_text(size = 11, face = "bold"))+ # Adjust title size and boldness
    xlab("Mean Difference Between Sexes")
  
################################################################################
###                               FINAL VISUAL                               ###
################################################################################  

## Conbine all figures together  
sim_site <- ggarrange(sim_a, sim_b, sim_c, sim_d, sim_e, sim_f, sim_g,
                  ncol = 2, nrow = 4)  

## Save output  
setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Tables and Figures/Figure 5")
ggsave("sim_site.png", dpi = 1200, width = 8, height = 8, units = "in")
  