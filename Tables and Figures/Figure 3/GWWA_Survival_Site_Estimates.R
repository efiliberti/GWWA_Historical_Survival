## Figure 3 in manuscript. This visualization shows site-specific mean apparent
## annual survival (both males and females). Points represent the mean posterior
## estimates. Lines represent 95% credibility intervals. 

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

## Obtain site-specific mean estimates and create a vector
AWV.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.AWV.phi) ## Northeastern West Virginia
CNY.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.CNY.phi) ## Southern New York
KNC.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.KNC.phi) ## Western North Carolina
PAC.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.PAC.phi) ## Central Pennsylvania
PAE.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.PAE.phi) ## Eastern Pennsylvania
SNC.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.SNC.phi) ## Northwestern North Carolina
CWV.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.CWV.phi) ## Southern West Virginia
BTN.mean <- mean(AppsMultiSiteNo2year$sims.list$overall.BTN.phi) ## Tennessee

App.Mean <- c(AWV.mean, CNY.mean, KNC.mean, PAC.mean, PAE.mean, SNC.mean, CWV.mean, BTN.mean)

## Obtain upper and lower 95% credible intervals for each site and put into vector
AWV.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.AWV.phi, probs = 0.975)
CNY.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.CNY.phi, probs = 0.975)
KNC.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.KNC.phi, probs = 0.975)
PAC.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.PAC.phi, probs = 0.975)
PAE.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.PAE.phi, probs = 0.975)
SNC.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.SNC.phi, probs = 0.975)
CWV.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.CWV.phi, probs = 0.975)
BTN.upper <- quantile(AppsMultiSiteNo2year$sims.list$overall.BTN.phi, probs = 0.975)

App.Upper <- c(AWV.upper, CNY.upper, KNC.upper, PAC.upper, PAE.upper, SNC.upper, CWV.upper, BTN.upper)

AWV.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.AWV.phi, probs = 0.025)
CNY.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.CNY.phi, probs = 0.025)
KNC.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.KNC.phi, probs = 0.025)
PAC.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.PAC.phi, probs = 0.025)
PAE.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.PAE.phi, probs = 0.025)
SNC.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.SNC.phi, probs = 0.025)
CWV.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.CWV.phi, probs = 0.025)
BTN.lower <- quantile(AppsMultiSiteNo2year$sims.list$overall.BTN.phi, probs = 0.025)

App.Lower <- c(AWV.lower, CNY.lower, KNC.lower, PAC.lower, PAE.lower, SNC.lower, CWV.lower, BTN.lower)

## Put mean and 95% intervals in one dataframe
App.Site <- c("West Virginia (NE)", "New York (S)", "North Carolina (W)", "Pennsylvania (C)", "Pennsylvania (E)", 
          "North Carolina (NW)", "West Virginia (S)", "Tennessee")

df_app <- data.frame(App.Site, App.Mean, App.Lower, App.Upper)

## Relevel to organize by geography in visualization 
df_app$App.Site <- factor(df_app$App.Site, levels = c("North Carolina (NW)", "North Carolina (W)", "Tennessee", "West Virginia (NE)", 
                                              "West Virginia (S)", "Pennsylvania (E)", "Pennsylvania (C)", "New York (S)"))

## Obtain sample sizes of all individuals for each site to include in viz
app_text_8 <- as.character(nrow(CNY_CH))
app_text_7 <- as.character(nrow(PAC_CH))
app_text_6 <- as.character(nrow(PAE_CH))
app_text_5 <- as.character(nrow(CWV_CH))
app_text_4 <- as.character(nrow(AWV_CH))
app_text_3 <- as.character(nrow(BTN_CH))
app_text_2 <- as.character(nrow(KNC_CH))
app_text_1 <- as.character(nrow(SNC_CH))

## Create visualization
appsseci <- ggplot(df_app, aes(x=App.Site, y=App.Mean)) +
  geom_point(size=2, shape = 16) +
  geom_errorbar(aes(ymin=App.Lower, ymax=App.Upper, width=0), size = 0.5) +
  ylab("Mean Apparent Annual Survival") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), labels = scales::number_format(accuracy = 0.1)) +
  coord_flip() +
  xlab("") +
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10, hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title = element_text(colour="white")) +
  annotate("text", label = app_text_1, x=1, y=1, size = 3.5, color = "grey40") +
  annotate("text", label = app_text_2, x=2, y=1, size = 3.5, color = "grey40") +
  annotate("text", label = app_text_3, x=3, y=1, size = 3.5, color = "grey40") +
  annotate("text", label = app_text_4, x=4, y=1, size = 3.5, color = "grey40") +
  annotate("text", label = app_text_5, x=5, y=1, size = 3.5, color = "grey40") +
  annotate("text", label = app_text_6, x=6, y=1, size = 3.5, color = "grey40") +
  annotate("text", label = app_text_7, x=7, y=1, size = 3.5, color = "grey40") +
  annotate("text", label = app_text_8, x=8, y=1, size = 3.5, color = "grey40") +
  ggtitle("Appalachians Site-specific Estimates")

## Save output
setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Tables and Figures/Figure 3")
ggsave("appsseci.png", dpi = 1200, width = 5, height = 4, units = "in")

################################################################################
###                               GREAT LAKES                                ###
################################################################################

## Obtain site-specific mean estimates and create a vector
CCR.mean <- mean(GLMultiSite$sims.list$overall.CCR.phi) ## Costa Rica
MMA.mean <- mean(GLMultiSite$sims.list$overall.MMA.phi) ## Southeastern Manitoba
RMI.mean <- mean(GLMultiSite$sims.list$overall.RMI.phi) ## Northern Michigan
RWI.mean <- mean(GLMultiSite$sims.list$overall.RWI.phi) ## Wisconsin
WMI.mean <- mean(GLMultiSite$sims.list$overall.WMI.phi) ## Central Michigan
VON.mean <- mean(GLMultiSite$sims.list$overall.VON.phi) ## Eastern Ontario
BGU.mean <- mean(GLMultiSite$sims.list$overall.BGU.phi) ## Guatemala
BHO.mean <- mean(GLMultiSite$sims.list$overall.BHO.phi) ## Honduras
CNI.mean <- mean(GLMultiSite$sims.list$overall.CNI.phi) ## Nicaragua
VMA.mean <- mean(GLMultiSite$sims.list$overall.VMA.phi) ## Southwestern Manitoba
BTN.mean <- mean(GLMultiSite$sims.list$overall.VONW.phi) ## Western Ontario

GL.Mean <- c(CCR.mean, MMA.mean, RMI.mean, RWI.mean, WMI.mean, VON.mean, BGU.mean, 
          BHO.mean, CNI.mean, VMA.mean, BTN.mean)

## Obtain upper and lower 95% credible intervals for each site and put into vector
CCR.upper <- quantile(GLMultiSite$sims.list$overall.CCR.phi, probs = 0.975)
MMA.upper <- quantile(GLMultiSite$sims.list$overall.MMA.phi, probs = 0.975)
RMI.upper <- quantile(GLMultiSite$sims.list$overall.RMI.phi, probs = 0.975)
RWI.upper <- quantile(GLMultiSite$sims.list$overall.RWI.phi, probs = 0.975)
WMI.upper <- quantile(GLMultiSite$sims.list$overall.WMI.phi, probs = 0.975)
VON.upper <- quantile(GLMultiSite$sims.list$overall.VON.phi, probs = 0.975)
BGU.upper <- quantile(GLMultiSite$sims.list$overall.BGU.phi, probs = 0.975)
BHO.upper <- quantile(GLMultiSite$sims.list$overall.BHO.phi, probs = 0.975)
CNI.upper <- quantile(GLMultiSite$sims.list$overall.CNI.phi, probs = 0.975)
VMA.upper <- quantile(GLMultiSite$sims.list$overall.VMA.phi, probs = 0.975)
VONW.upper <- quantile(GLMultiSite$sims.list$overall.VONW.phi, probs = 0.975)

GL.Upper <- c(CCR.upper, MMA.upper, RMI.upper, RWI.upper, WMI.upper, VON.upper, BGU.upper, 
          BHO.upper, CNI.upper, VMA.upper, BTN.upper)

CCR.lower <- quantile(GLMultiSite$sims.list$overall.CCR.phi, probs = 0.025)
MMA.lower <- quantile(GLMultiSite$sims.list$overall.MMA.phi, probs = 0.025)
RMI.lower <- quantile(GLMultiSite$sims.list$overall.RMI.phi, probs = 0.025)
RWI.lower <- quantile(GLMultiSite$sims.list$overall.RWI.phi, probs = 0.025)
WMI.lower <- quantile(GLMultiSite$sims.list$overall.WMI.phi, probs = 0.025)
VON.lower <- quantile(GLMultiSite$sims.list$overall.VON.phi, probs = 0.025)
BGU.lower <- quantile(GLMultiSite$sims.list$overall.BGU.phi, probs = 0.025)
BHO.lower <- quantile(GLMultiSite$sims.list$overall.BHO.phi, probs = 0.025)
CNI.lower <- quantile(GLMultiSite$sims.list$overall.CNI.phi, probs = 0.025)
VMA.lower <- quantile(GLMultiSite$sims.list$overall.VMA.phi, probs = 0.025)
VONW.lower <- quantile(GLMultiSite$sims.list$overall.VONW.phi, probs = 0.025)

GL.Lower <- c(CCR.lower, MMA.lower, RMI.lower, RWI.lower, WMI.lower, VON.lower, BGU.lower, 
           BHO.lower, CNI.lower, VMA.lower, BTN.lower)

## Put mean and 95% intervals in one dataframe
GL.Site <- c("Costa Rica", "Manitoba (SE)", "Michigan (N)", "Wisconsin", "Michigan (C)",
          "Ontario (E)", "Guatemala", "Honduras", "Nicaragua", "Manitoba (SW)",
          "Ontario (W)")

df_gl <- data.frame(GL.Site, GL.Mean, GL.Lower, GL.Upper)

## Relevel to organize by geography in visualization 
df_gl$GL.Site <- factor(df_gl$GL.Site, levels = c("Nicaragua", "Honduras", "Guatemala", "Costa Rica",
                                            "Ontario (E)", "Michigan (C)", "Michigan (N)",
                                            "Wisconsin", "Ontario (W)", "Manitoba (SW)", "Manitoba (SE)"))

## Obtain sample sizes of all individuals for each site to include in viz
gl_text_11 <- as.character(nrow(MMA_CH))
gl_text_10 <- as.character(nrow(VMA_CH))
gl_text_9 <- as.character(nrow(VONW_CH))
gl_text_8 <- as.character(nrow(RWI_CH))
gl_text_7 <- as.character(nrow(RMI_CH))
gl_text_6 <- as.character(nrow(WMI_CH))
gl_text_5 <- as.character(nrow(VON_CH))
gl_text_4 <- as.character(nrow(CCR_CH))
gl_text_3 <- as.character(nrow(BGU_CH))
gl_text_2 <- as.character(nrow(BHO_CH))
gl_text_1 <- as.character(nrow(CNI_CH))

## Create visualization
glsseci <- ggplot(df_gl, aes(x=GL.Site, y=GL.Mean))+
  geom_point(size=2, shape = 16)+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), labels = scales::number_format(accuracy = 0.1)) +
  geom_errorbar(aes(ymin=GL.Lower, ymax=GL.Upper, width=0), size = 0.5)+
  coord_flip()+
  ylab("Mean Apparent Annual Survival")+
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
  annotate("text", label = gl_text_1, x=1, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = gl_text_2, x=2, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = gl_text_3, x=3, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = gl_text_4, x=4, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = gl_text_5, x=5, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = gl_text_6, x=6, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = gl_text_7, x=7, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = gl_text_8, x=8, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = gl_text_9, x=9, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = gl_text_10, x=10, y=1, size = 3.5, color = "grey40")+
  annotate("text", label = gl_text_11, x=11, y=1, size = 3.5, color = "grey40")+
  ggtitle("Great Lakes Site-specific Estimates")+
  theme(
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent"),
    panel.grid.major = element_line(color = "transparent"),
    panel.grid.minor = element_line(color = "transparent")
    )

## Save output
setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Tables and Figures/Figure 3")
ggsave("glsseci.png", dpi = 1200, width = 5, height = 4, units = "in")
