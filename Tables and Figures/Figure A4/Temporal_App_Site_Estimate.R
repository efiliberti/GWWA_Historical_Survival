## Figure A4 in manuscript. This figure shows annual survival estimates for each
## time period in each site representing the Great Lakes population.

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
###                        INDIVIDUAL SITE TEMPORAL PLOTS                    ###
################################################################################

## Tennessee
phi.BTN <- as.data.frame(AppsMultiSiteNo2year$sims.list$phi.BTN)
colnames(phi.BTN) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b",
                       "5a", "5b", "6a", "6b", "7a", "7b", "8a", "8b")


phi.BTN <- phi.BTN  %>%
  pivot_longer(1:16)

names.BTN <- c("1a" = "2004", "1b" = "2004",
               "2a" = "2005", "2b" = "2005",
               "3a" = "2006", "3b" = "2006",
               "4a" = "2007", "4b" = "2007",
               "5a" = "2008", "5b" = "2008",
               "6a" = "2009", "6b" = "2009",
               "7a" = "2010", "7b" = "2010",
               "8a" = "2011", "8b" = "2011")

phi.BTN <- phi.BTN  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.BTN) ~ names.BTN[name],
    TRUE ~ name
  ))



phi.BTN$site <- "Tennessee"

phi.BTN$name <- as.factor(phi.BTN$name)

label_positions <- seq(1, length(levels(phi.BTN$name)), by = 2)


s1 <- ggplot(phi.BTN, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Tennessee")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.BTN$name),
    labels = ifelse(seq_along(levels(phi.BTN$name)) %in% label_positions, levels(phi.BTN$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.BTN$name)) %in% label_positions, 1, 0.5))
  )

## West Virginia (NE)
phi.AWV <- as.data.frame(AppsMultiSiteNo2year$sims.list$phi.AWV)
colnames(phi.AWV) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b",
                       "5a", "5b", "6a", "6b", "7a", "7b")


phi.AWV <- phi.AWV  %>%
  pivot_longer(1:14)

names.AWV <- c("1a" = "2009", "1b" = "2009",
               "2a" = "2010", "2b" = "2010",
               "3a" = "2011", "3b" = "2011",
               "4a" = "2012", "4b" = "2012",
               "5a" = "2013", "5b" = "2013",
               "6a" = "2014", "6b" = "2014",
               "7a" = "2015", "7b" = "2015",
               "8a" = "2016", "8b" = "2016")

phi.AWV <- phi.AWV  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.AWV) ~ names.AWV[name],
    TRUE ~ name
  ))



phi.AWV$site <- "West Virginia (NE)"

phi.AWV$name <- as.factor(phi.AWV$name)

label_positions <- seq(1, length(levels(phi.AWV$name)), by = 2)

s2 <- ggplot(phi.AWV, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("West Virginia (NE)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.AWV$name),
    labels = ifelse(seq_along(levels(phi.AWV$name)) %in% label_positions, levels(phi.AWV$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.AWV$name)) %in% label_positions, 1, 0.5))
  )


## New York (S)
phi.CNY <- as.data.frame(AppsMultiSiteNo2year$sims.list$phi.CNY)
colnames(phi.CNY) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b",
                       "5a", "5b", "6a", "6b", "7a", "7b", "8a", "8b")


phi.CNY <- phi.CNY  %>%
  pivot_longer(1:16)

names.CNY <- c("1a" = "2003", "1b" = "2003",
               "2a" = "2004", "2b" = "2004",
               "3a" = "2005", "3b" = "2005",
               "4a" = "2006", "4b" = "2006",
               "5a" = "2007", "5b" = "2007",
               "6a" = "2008", "6b" = "2008",
               "7a" = "2009", "7b" = "2009",
               "8a" = "2010", "8b" = "2010")

phi.CNY <- phi.CNY  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.CNY) ~ names.CNY[name],
    TRUE ~ name
  ))



phi.CNY$site <- "New York (S)"

phi.CNY$name <- as.factor(phi.CNY$name)

label_positions <- seq(1, length(levels(phi.CNY$name)), by = 2)

s3 <- ggplot(phi.CNY, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("New York (S)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.CNY$name),
    labels = ifelse(seq_along(levels(phi.CNY$name)) %in% label_positions, levels(phi.CNY$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.CNY$name)) %in% label_positions, 1, 0.5))
  )


## North Carolina (W)
phi.KNC <- as.data.frame(AppsMultiSiteNo2year$sims.list$phi.KNC)
colnames(phi.KNC) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b",
                       "5a", "5b")


phi.KNC <- phi.KNC  %>%
  pivot_longer(1:10)

names.KNC <- c("1a" = "2009", "1b" = "2009",
               "2a" = "2010", "2b" = "2010",
               "3a" = "2011", "3b" = "2011",
               "4a" = "2020", "4b" = "2020",
               "5a" = "2021", "5b" = "2021")

phi.KNC <- phi.KNC  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.KNC) ~ names.KNC[name],
    TRUE ~ name
  ))



phi.KNC$site <- "North Carolina (W)"

phi.KNC$name <- as.factor(phi.KNC$name)

label_positions <- seq(1, length(levels(phi.KNC$name)), by = 2)

s4 <- ggplot(phi.KNC, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("North Carolina (W)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.KNC$name),
    labels = ifelse(seq_along(levels(phi.KNC$name)) %in% label_positions, levels(phi.KNC$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.KNC$name)) %in% label_positions, 1, 0.5))
  )


## Pennsylvania (C)
phi.PAC <- as.data.frame(AppsMultiSiteNo2year$sims.list$phi.PAC)
colnames(phi.PAC) <- c("1a", "1b", "2a", "2b", "3a", "3b")


phi.PAC <- phi.PAC  %>%
  pivot_longer(1:6)

names.PAC <- c("1a" = "2009", "1b" = "2009",
               "2a" = "2010", "2b" = "2010",
               "3a" = "2011", "3b" = "2011")

phi.PAC <- phi.PAC  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.PAC) ~ names.PAC[name],
    TRUE ~ name
  ))



phi.PAC$site <- "Pennsylvania (C)"

phi.PAC$name <- as.factor(phi.PAC$name)

label_positions <- seq(1, length(levels(phi.PAC$name)), by = 2)

s5 <- ggplot(phi.PAC, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Pennsylvania (C)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.PAC$name),
    labels = ifelse(seq_along(levels(phi.PAC$name)) %in% label_positions, levels(phi.PAC$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.PAC$name)) %in% label_positions, 1, 0.5))
  )




## North Carolina (NW)
phi.SNC <- as.data.frame(AppsMultiSiteNo2year$sims.list$phi.SNC)
colnames(phi.SNC) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b",
                       "5a", "5b")


phi.SNC <- phi.SNC  %>%
  pivot_longer(1:10)

names.SNC <- c("1a" = "2011", "1b" = "2011",
               "2a" = "2012", "2b" = "2012",
               "3a" = "2013", "3b" = "2013",
               "4a" = "2014", "4b" = "2014",
               "5a" = "2015", "5b" = "2015")

phi.SNC <- phi.SNC  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.SNC) ~ names.SNC[name],
    TRUE ~ name
  ))



phi.SNC$site <- "North Carolina (NW)"

phi.SNC$name <- as.factor(phi.SNC$name)

label_positions <- seq(1, length(levels(phi.SNC$name)), by = 2)

s6 <- ggplot(phi.SNC, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("North Carolina (NW)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.SNC$name),
    labels = ifelse(seq_along(levels(phi.SNC$name)) %in% label_positions, levels(phi.SNC$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.SNC$name)) %in% label_positions, 1, 0.5))
  )


## West Virginia (S)
phi.CWV <- as.data.frame(AppsMultiSiteNo2year$sims.list$phi.CWV)
colnames(phi.CWV) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b",
                       "5a", "5b", "6a", "6b", "7a", "7b", "8a", "8b",
                       "9a", "9b", "10a", "10b", "11a", "11b", "12a", "12b",
                       "13a", "13b", "14a", "14b", "15a", "15b", "16a", "16b",
                       "17a", "17b", "18a", "18b", "19a", "19b", "20a", "20b",
                       "21a", "21b")


phi.CWV <- phi.CWV  %>%
  pivot_longer(1:42)

names.CWV <- c("1a" = "2001", "1b" = "2001",
               "2a" = "2002", "2b" = "2002",
               "3a" = "2003", "3b" = "2003",
               "4a" = "2004", "4b" = "2004",
               "5a" = "2005", "5b" = "2005",
               "6a" = "2006", "6b" = "2006",
               "7a" = "2007", "7b" = "2007",
               "8a" = "2008", "8b" = "2008",
               "9a" = "2009", "9b" = "2009",
               "10a" = "2010", "10b" = "2010",
               "11a" = "2011", "11b" = "2011",
               "12a" = "2012", "12b" = "2012",
               "13a" = "2013", "13b" = "2013",
               "14a" = "2014", "14b" = "2014",
               "15a" = "2015", "15b" = "2015",
               "16a" = "2016", "16b" = "2016",
               "17a" = "2017", "17b" = "2017",
               "18a" = "2018", "18b" = "2018",
               "19a" = "2019", "19b" = "2019",
               "20a" = "2020", "20b" = "2020",
               "21a" = "2021", "21b" = "2021")

phi.CWV <- phi.CWV  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.CWV) ~ names.CWV[name],
    TRUE ~ name
  ))



phi.CWV$site <- "West Virginia (S)"

phi.CWV$name <- as.factor(phi.CWV$name)

label_positions <- seq(1, length(levels(phi.CWV$name)), by = 2)

s7 <- ggplot(phi.CWV, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("West Virginia (S)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.CWV$name),
    labels = ifelse(seq_along(levels(phi.CWV$name)) %in% label_positions, levels(phi.CWV$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.CWV$name)) %in% label_positions, 1, 0.5))
  )



## Pennsylvania (E)
phi.PAE <- as.data.frame(AppsMultiSiteNo2year$sims.list$phi.PAE)
colnames(phi.PAE) <- c("1a", "2a", "3a")


phi.PAE <- phi.PAE  %>%
  pivot_longer(1:3)

names.PAE <- c("1a" = "2013", "2a" = "2014", "3a" = "2015")

phi.PAE <- phi.PAE  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.PAE) ~ names.PAE[name],
    TRUE ~ name
  ))



phi.PAE$site <- "Pennsylvania (E)"

phi.PAE$name <- as.factor(phi.PAE$name)

label_positions <- seq(1, length(levels(phi.PAE$name)), by = 2)

s8 <- ggplot(phi.PAE, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Pennsylvania (E)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.PAE$name),
    labels = ifelse(seq_along(levels(phi.PAE$name)) %in% label_positions, levels(phi.PAE$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.PAE$name)) %in% label_positions, 1, 0.5))
  )


################################################################################
###                       COMBINE ALL PLOTS FOR FINAL PLOT                   ###
################################################################################

aa <- ggarrange(s1, s2, s3, s4, s5, s6, s8,
                ncol = 3, nrow = 3) 

aa1 <- ggarrange(aa, s7,
                 ncol = 1,
                 heights = c(2.5, 1),
                 labels = c("", "H")) 

## Save output
setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Tables and Figures/Figure A4")
ggsave("aa1.png", dpi = 1200, width = 8.5, height = 15, units = "in")
