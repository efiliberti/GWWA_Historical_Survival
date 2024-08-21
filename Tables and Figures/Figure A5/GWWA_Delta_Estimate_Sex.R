## Figure A5 in manuscript. This visualization shows delta estimates comparing 
## male and female apparent annual survival at each respective site

## Load packages
library(ggplot2)
library(lme4)
library(sjPlot)
library(dplyr)
library(performance)
library(cowplot)

################################################################################
###                              APPALACHIANS                                ###
################################################################################

## Get mean delta estimates for Appalachian population
AWV.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.delta.phi.AWV)
CNY.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.delta.phi.CNY)
KNC.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.delta.phi.KNC)
PAC.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.delta.phi.PAC)
SNC.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.delta.phi.SNC)
CWV.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.delta.phi.CWV)
BTN.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.delta.phi.BTN)

App.Delta <- c(AWV.mean, CNY.mean, KNC.mean, PAC.mean, SNC.mean, CWV.mean, BTN.mean)

## Obtain upper and lower 95% credible intervals for each site and put into vector
AWV.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.AWV, probs = 0.975)
CNY.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.CNY, probs = 0.975)
KNC.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.KNC, probs = 0.975)
PAC.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.PAC, probs = 0.975)
SNC.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.SNC, probs = 0.975)
CWV.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.CWV, probs = 0.975)
BTN.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.BTN, probs = 0.975)

App.Upper <- c(AWV.upper, CNY.upper, KNC.upper, PAC.upper, SNC.upper, CWV.upper, BTN.upper)

AWV.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.AWV, probs = 0.025)
CNY.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.CNY, probs = 0.025)
KNC.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.KNC, probs = 0.025)
PAC.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.PAC, probs = 0.025)
SNC.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.SNC, probs = 0.025)
CWV.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.CWV, probs = 0.025)
BTN.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.delta.phi.BTN, probs = 0.025)

App.Lower <- c(AWV.lower, CNY.lower, KNC.lower, PAC.lower, SNC.lower, CWV.lower, BTN.lower)

## Create site names for visualization
App.Site <- c("West Virginia (NE)", "New York (S)", "North Carolina (W)", "Pennsylvania (C)", 
              "North Carolina (NW)", "West Virginia (S)", "Tennessee")

App.Code <- c("AWV", "CNY", "KNC", "PAC", "SNC", "CWV", "BTN")

df_app_sex <- data.frame(App.Site, App.Delta, App.Lower, App.Upper, App.Code)

## Relevel for viz purposes
df_app_sex$App.Site <- factor(df_app_sex$App.Site, levels = c("North Carolina (NW)", "North Carolina (W)", 
                                                  "Tennessee", "West Virginia (NE)", "West Virginia (S)", 
                                                  "Pennsylvania (C)", "New York (S)"))

## Find the percentage of female representation to put in visual
app_sex_count <- cbind(app_sex_count, Total = app_sex_count[, "males"] + app_sex_count[, "females"])
app_sex_count$Percent_Female <- round((app_sex_count$females / app_sex_count$Total) * 100, 1)

app_rep_1 <- as.character(subset(app_sex_count, app_sex_labels == "SNC")$Percent_Female)
app_rep_2 <- as.character(subset(app_sex_count, app_sex_labels == "KNC")$Percent_Female)
app_rep_3 <- as.character(subset(app_sex_count, app_sex_labels == "BTN")$Percent_Female)
app_rep_4 <- as.character(subset(app_sex_count, app_sex_labels == "AWV")$Percent_Female)
app_rep_5 <- as.character(subset(app_sex_count, app_sex_labels == "CWV")$Percent_Female)
app_rep_6 <- as.character(subset(app_sex_count, app_sex_labels == "PAC")$Percent_Female)
app_rep_7 <- as.character(subset(app_sex_count, app_sex_labels == "CNY")$Percent_Female)



sexefapp <- ggplot(df_app_sex, aes(x=App.Site, y=App.Delta))+
  geom_point(size=2, shape = 16)+
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2), labels = scales::number_format(accuracy = 0.1)) +
  geom_errorbar(aes(ymin=App.Lower, ymax=App.Upper, width=0), size = 0.5)+
  coord_flip()+
  ylab("Mean Difference Between Male and Female Survival Estimates")+
  xlab("")+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  annotate("text", label = app_rep_1, x=1, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = app_rep_2, x=2, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = app_rep_3, x=3, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = app_rep_4, x=4, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = app_rep_5, x=5, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = app_rep_6, x=6, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = app_rep_7, x=7, y=1, size = 3.5, color = "grey40")+
  ggtitle("Appalachian Delta Sex Estimates")+
  geom_hline(aes(yintercept=0), col = "red", linetype = 'dashed', size = 0.75)

## Save output
setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Tables and Figures/Figure A5")
ggsave("sexefapp.png", dpi = 1200, width = 5, height = 4, units = "in")

################################################################################
###                               GREAT LAKES                                ###
################################################################################

## Get mean overall.delta estimates for Great Lakes population
CCR.mean <- mean(GLMultiSiteNo2year$sims.list$overall.delta.phi.CCR)
MMA.mean <- mean(GLMultiSiteNo2year$sims.list$overall.delta.phi.MMA)
RMI.mean <- mean(GLMultiSiteNo2year$sims.list$overall.delta.phi.RMI)
RWI.mean <- mean(GLMultiSiteNo2year$sims.list$overall.delta.phi.RWI)
WMI.mean <- mean(GLMultiSiteNo2year$sims.list$overall.delta.phi.WMI)
VON.mean <- mean(GLMultiSiteNo2year$sims.list$overall.delta.phi.VON)
BGU.mean <- mean(GLMultiSiteNo2year$sims.list$overall.delta.phi.BGU)
BHO.mean <- mean(GLMultiSiteNo2year$sims.list$overall.delta.phi.BHO)
CNI.mean <- mean(GLMultiSiteNo2year$sims.list$overall.delta.phi.CNI)
VMA.mean <- mean(GLMultiSiteNo2year$sims.list$overall.delta.phi.VMA)

GL.Delta <- c(CCR.mean, MMA.mean, RMI.mean, RWI.mean, WMI.mean, VON.mean, BGU.mean,
              BHO.mean, CNI.mean, VMA.mean)

## Obtain upper and lower 95% credible intervals for each site and put into vector
CCR.upper <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.CCR, probs = 0.975)
MMA.upper <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.MMA, probs = 0.975)
RMI.upper <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.RMI, probs = 0.975)
RWI.upper <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.RWI, probs = 0.975)
WMI.upper <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.WMI, probs = 0.975)
VON.upper <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.VON, probs = 0.975)
BGU.upper <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.BGU, probs = 0.975)
BHO.upper <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.BHO, probs = 0.975)
CNI.upper <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.CNI, probs = 0.975)
VMA.upper <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.VMA, probs = 0.975)

GL.Upper <- c(CCR.upper, MMA.upper, RMI.upper, RWI.upper, WMI.upper, VON.upper, BGU.upper,
              BHO.upper, CNI.upper, VMA.upper)

CCR.lower <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.CCR, probs = 0.025)
MMA.lower <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.MMA, probs = 0.025)
RMI.lower <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.RMI, probs = 0.025)
RWI.lower <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.RWI, probs = 0.025)
WMI.lower <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.WMI, probs = 0.025)
VON.lower <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.VON, probs = 0.025)
BGU.lower <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.BGU, probs = 0.025)
BHO.lower <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.BHO, probs = 0.025)
CNI.lower <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.CNI, probs = 0.025)
VMA.lower <- quantile(GLMultiSiteNo2year$sims.list$overall.delta.phi.VMA, probs = 0.025)

GL.Lower <- c(CCR.lower, MMA.lower, RMI.lower, RWI.lower, WMI.lower, VON.lower, BGU.lower,
              BHO.lower, CNI.lower, VMA.lower)

## Create site names for visualization
GL.Site <- c("Costa Rica", "Manitoba (SE)", "Michigan (N)", "Wisconsin", "Michigan (C)",
          "Ontario (E)", "Guatemala", "Honduras", "Nicaragua", "Manitoba (SW)")

GL.Code <- c("CCR", "MMA", "RMI", "RWI", "WMI", "VON", "BGU", "BHO", "CNI", "VMA")

df_gl_sex <- data.frame(GL.Site, GL.Delta, GL.Lower, GL.Upper, GL.Code)

## Create site names for visualization
df_gl_sex$GL.Site <- factor(df_gl_sex$GL.Site, levels = c("Nicaragua", "Honduras", "Guatemala", "Costa Rica",
                                            "Ontario (E)", "Michigan (C)", "Michigan (N)",
                                            "Wisconsin", "Manitoba (SW)", "Manitoba (SE)"))

## Find the percentage of female representation to put in visual
gl_sex_count <- cbind(gl_sex_count, Total = gl_sex_count[, "males"] + gl_sex_count[, "females"])
gl_sex_count$Percent_Female <- round((gl_sex_count$females / gl_sex_count$Total) * 100, 1)

gl_rep_1 <- as.character(subset(gl_sex_count, gl_sex_labels == "CNI")$Percent_Female)
gl_rep_2 <- as.character(subset(gl_sex_count, gl_sex_labels == "BHO")$Percent_Female)
gl_rep_3 <- as.character(subset(gl_sex_count, gl_sex_labels == "BGU")$Percent_Female)
gl_rep_4 <- as.character(subset(gl_sex_count, gl_sex_labels == "CCR")$Percent_Female)
gl_rep_5 <- as.character(subset(gl_sex_count, gl_sex_labels == "VON")$Percent_Female)
gl_rep_6 <- as.character(subset(gl_sex_count, gl_sex_labels == "WMI")$Percent_Female)
gl_rep_7 <- as.character(subset(gl_sex_count, gl_sex_labels == "RMI")$Percent_Female)
gl_rep_8 <- as.character(subset(gl_sex_count, gl_sex_labels == "RWI")$Percent_Female)
gl_rep_9 <- as.character(subset(gl_sex_count, gl_sex_labels == "VMA")$Percent_Female)
gl_rep_10 <- as.character(subset(gl_sex_count, gl_sex_labels == "MMA")$Percent_Female)


sexefgl <- ggplot(df_gl_sex, aes(x=GL.Site, y=GL.Delta))+
  geom_point(size=2, shape = 16)+
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2), labels = scales::number_format(accuracy = 0.1)) +
  geom_errorbar(aes(ymin=GL.Lower, ymax=GL.Upper, width=0), size = 0.5)+
  coord_flip()+
  ylab("Mean Difference Between Male and Female Survival Estimates")+
  xlab("")+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  annotate("text", label = gl_rep_1, x=1, y=1.00, size = 3.5, color = "grey40")+
  annotate("text", label = gl_rep_2, x=2, y=1.00, size = 3.5, color = "grey40")+
  annotate("text", label = gl_rep_3, x=3, y=1.00, size = 3.5, color = "grey40")+
  annotate("text", label = gl_rep_4, x=4, y=1.00, size = 3.5, color = "grey40")+
  annotate("text", label = gl_rep_5, x=5, y=1.00, size = 3.5, color = "grey40")+
  annotate("text", label = gl_rep_6, x=6, y=1.00, size = 3.5, color = "grey40")+
  annotate("text", label = gl_rep_7, x=7, y=1.00, size = 3.5, color = "grey40")+
  annotate("text", label = gl_rep_8, x=8, y=1.00, size = 3.5, color = "grey40")+
  annotate("text", label = gl_rep_9, x=9, y=1.00, size = 3.5, color = "grey40")+
  annotate("text", label = gl_rep_10, x=10, y=1.00, size = 3.5, color = "grey40")+
  ggtitle("Great Lakes Site-specific Estimates")+
  geom_hline(aes(yintercept=0), col = "red", linetype = 'dashed', size = 0.75)

## Save output
setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Tables and Figures/Figure A5")
ggsave("sexefgl.png", dpi = 1200, width = 5, height = 4, units = "in")
