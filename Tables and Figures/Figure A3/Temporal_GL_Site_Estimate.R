## Figure A3 in manuscript. This figure shows annual survival estimates for each
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

## Costa Rica
phi.CCR <- as.data.frame(GLMultiSiteNo2year$sims.list$phi.CCR)
colnames(phi.CCR) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b",
                       "5a", "5b", "6a", "6b", "7a", "7b", "8a", "8b",
                       "9a", "9b", "10a", "10b")


phi.CCR <- phi.CCR  %>%
  pivot_longer(1:20)

names.CCR <- c("1a" = "2007", "1b" = "2007",
               "2a" = "2008", "2b" = "2008",
               "3a" = "2009", "3b" = "2009",
               "4a" = "2010", "4b" = "2010",
               "5a" = "2011", "5b" = "2011",
               "6a" = "2012", "6b" = "2012",
               "7a" = "2013", "7b" = "2013",
               "8a" = "2014", "8b" = "2014",
               "9a" = "2015", "9b" = "2015",
               "10a" = "2016", "10b" = "2016")

phi.CCR <- phi.CCR  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.CCR) ~ names.CCR[name],
    TRUE ~ name
  ))



phi.CCR$site <- "Costa Rica"

phi.CCR$name <- as.factor(phi.CCR$name)

label_positions <- seq(1, length(levels(phi.CCR$name)), by = 2)

sa <- ggplot(phi.CCR, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Costa Rica")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.CCR$name),
    labels = ifelse(seq_along(levels(phi.CCR$name)) %in% label_positions, levels(phi.CCR$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.CCR$name)) %in% label_positions, 1, 0.5))
  )

## Manitoba (SE)
phi.MMA <- as.data.frame(GLMultiSiteNo2year$sims.list$phi.MMA)
colnames(phi.MMA) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b",
                       "5a", "5b", "6a", "6b")


phi.MMA <- phi.MMA  %>%
  pivot_longer(1:12)

names.MMA <- c("1a" = "2012", "1b" = "2012",
               "2a" = "2013", "2b" = "2013",
               "3a" = "2014", "3b" = "2014",
               "4a" = "2015", "4b" = "2015",
               "5a" = "2016", "5b" = "2016",
               "6a" = "2017", "6b" = "2017")

phi.MMA <- phi.MMA  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.MMA) ~ names.MMA[name],
    TRUE ~ name
  ))



phi.MMA$site <- "Manitoba (SE)"

phi.MMA$name <- as.factor(phi.MMA$name)

label_positions <- seq(1, length(levels(phi.MMA$name)), by = 2)

sb <- ggplot(phi.MMA, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Manitoba (SE)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.MMA$name),
    labels = ifelse(seq_along(levels(phi.MMA$name)) %in% label_positions, levels(phi.MMA$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.MMA$name)) %in% label_positions, 1, 0.5))
  )


## Michigan (N)
phi.RMI <- as.data.frame(GLMultiSiteNo2year$sims.list$phi.RMI)
colnames(phi.RMI) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b")


phi.RMI <- phi.RMI  %>%
  pivot_longer(1:8)

names.RMI <- c("1a" = "2014", "1b" = "2014",
               "2a" = "2015", "2b" = "2015",
               "3a" = "2016", "3b" = "2016",
               "4a" = "2017", "4b" = "2017")

phi.RMI <- phi.RMI  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.RMI) ~ names.RMI[name],
    TRUE ~ name
  ))



phi.RMI$site <- "Michigan (N)"

phi.RMI$name <- as.factor(phi.RMI$name)

label_positions <- seq(1, length(levels(phi.RMI$name)), by = 2)

sc <- ggplot(phi.RMI, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Michigan (N)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.RMI$name),
    labels = ifelse(seq_along(levels(phi.RMI$name)) %in% label_positions, levels(phi.RMI$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.RMI$name)) %in% label_positions, 1, 0.5))
  )

## Wisconsin
phi.RWI <- as.data.frame(GLMultiSiteNo2year$sims.list$phi.RWI)
colnames(phi.RWI) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b",
                       "5a", "5b", "6a", "6b", "7a", "7b", "8a", "8b",
                       "9a", "9b", "10a", "10b", "11a", "11b", "12a", "12b",
                       "13a", "13b")


phi.RWI <- phi.RWI  %>%
  pivot_longer(1:26)

names.RWI <- c("1a" = "2008", "1b" = "2008",
               "2a" = "2009", "2b" = "2009",
               "3a" = "2010", "3b" = "2010",
               "4a" = "2011", "4b" = "2011",
               "5a" = "2012", "5b" = "2012",
               "6a" = "2013", "6b" = "2013",
               "7a" = "2014", "7b" = "2014",
               "8a" = "2015", "8b" = "2015",
               "9a" = "2016", "9b" = "2016",
               "10a" = "2017", "10b" = "2017",
               "11a" = "2018", "11b" = "2018",
               "12a" = "2019", "12b" = "2019",
               "13a" = "2020", "13b" = "2020")

phi.RWI <- phi.RWI  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.RWI) ~ names.RWI[name],
    TRUE ~ name
  ))



phi.RWI$site <- "Wisconsin"

phi.RWI$name <- as.factor(phi.RWI$name)

label_positions <- seq(1, length(levels(phi.RWI$name)), by = 2)

sd <- ggplot(phi.RWI, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Wisconsin")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.RWI$name),
    labels = ifelse(seq_along(levels(phi.RWI$name)) %in% label_positions, levels(phi.RWI$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.RWI$name)) %in% label_positions, 1, 0.5))
  )

## Michigan (C)
phi.WMI <- as.data.frame(GLMultiSiteNo2year$sims.list$phi.WMI)
colnames(phi.WMI) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b")


phi.WMI <- phi.WMI  %>%
  pivot_longer(1:8)

names.WMI <- c("1a" = "1981", "1b" = "1981",
               "2a" = "1982", "2b" = "1982",
               "3a" = "1983", "3b" = "1983",
               "4a" = "1984", "4b" = "1984")

phi.WMI <- phi.WMI  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.WMI) ~ names.WMI[name],
    TRUE ~ name
  ))



phi.WMI$site <- "Michigan (C)"

phi.WMI$name <- as.factor(phi.WMI$name)

label_positions <- seq(1, length(levels(phi.WMI$name)), by = 2)

se <- ggplot(phi.WMI, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Michigan (C)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.WMI$name),
    labels = ifelse(seq_along(levels(phi.WMI$name)) %in% label_positions, levels(phi.WMI$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.WMI$name)) %in% label_positions, 1, 0.5))
  )


## Ontario (E)
phi.VON <- as.data.frame(GLMultiSiteNo2year$sims.list$phi.VON)
colnames(phi.VON) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b",
                       "5a", "5b")


phi.VON <- phi.VON  %>%
  pivot_longer(1:10)

names.VON <- c("1a" = "2002", "1b" = "2002",
               "2a" = "2003", "2b" = "2003",
               "3a" = "2004", "3b" = "2004",
               "4a" = "2005", "4b" = "2005",
               "5a" = "2006", "5b" = "2006")

phi.VON <- phi.VON  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.VON) ~ names.VON[name],
    TRUE ~ name
  ))



phi.VON$site <- "Ontario (E)"

phi.VON$name <- as.factor(phi.VON$name)

label_positions <- seq(1, length(levels(phi.VON$name)), by = 2)

sf <- ggplot(phi.VON, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Ontario (E)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.VON$name),
    labels = ifelse(seq_along(levels(phi.VON$name)) %in% label_positions, levels(phi.VON$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.VON$name)) %in% label_positions, 1, 0.5))
  )


## Guatemala
phi.BGU <- as.data.frame(GLMultiSiteNo2year$sims.list$phi.BGU)
colnames(phi.BGU) <- c("1a", "1b", "2a", "2b")


phi.BGU <- phi.BGU  %>%
  pivot_longer(1:4)

names.BGU <- c("1a" = "2015", "1b" = "2015",
               "2a" = "2016", "2b" = "2016")

phi.BGU <- phi.BGU  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.BGU) ~ names.BGU[name],
    TRUE ~ name
  ))



phi.BGU$site <- "Guatemala"

phi.BGU$name <- as.factor(phi.BGU$name)

label_positions <- seq(1, length(levels(phi.BGU$name)), by = 2)

sg <- ggplot(phi.BGU, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Guatemala")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.BGU$name),
    labels = ifelse(seq_along(levels(phi.BGU$name)) %in% label_positions, levels(phi.BGU$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.BGU$name)) %in% label_positions, 1, 0.5))
  )


## Honduras
phi.BHO <- as.data.frame(GLMultiSiteNo2year$sims.list$phi.BHO)
colnames(phi.BHO) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b",
                       "5a", "5b")


phi.BHO <- phi.BHO  %>%
  pivot_longer(1:10)

names.BHO <- c("1a" = "2012", "1b" = "2012",
               "2a" = "2013", "2b" = "2013",
               "3a" = "2014", "3b" = "2014",
               "4a" = "2015", "4b" = "2015",
               "5a" = "2016", "5b" = "2016")

phi.BHO <- phi.BHO  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.BHO) ~ names.BHO[name],
    TRUE ~ name
  ))



phi.BHO$site <- "Honduras"

phi.BHO$name <- as.factor(phi.BHO$name)

label_positions <- seq(1, length(levels(phi.BHO$name)), by = 2)

sh <- ggplot(phi.BHO, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Honduras")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.BHO$name),
    labels = ifelse(seq_along(levels(phi.BHO$name)) %in% label_positions, levels(phi.BHO$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.BHO$name)) %in% label_positions, 1, 0.5))
  )

## Nicaragua
phi.CNI <- as.data.frame(GLMultiSiteNo2year$sims.list$phi.CNI)
colnames(phi.CNI) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b",
                       "5a", "5b", "6a", "6b", "7a", "7b", "8a", "8b",
                       "9a", "9b", "10a", "10b", "11a", "11b", "12a", "12b",
                       "13a", "13b", "14a", "14b", "15a", "15b", "16a", "16b",
                       "17a", "17b", "18a", "18b", "19a", "19b", "20a", "20b")


phi.CNI <- phi.CNI  %>%
  pivot_longer(1:40)

names.CNI <- c("1a" = "2003", "1b" = "2003",
               "2a" = "2004", "2b" = "2004",
               "3a" = "2005", "3b" = "2005",
               "4a" = "2006", "4b" = "2006",
               "5a" = "2007", "5b" = "2007",
               "6a" = "2008", "6b" = "2008",
               "7a" = "2009", "7b" = "2009",
               "8a" = "2010", "8b" = "2010",
               "9a" = "2011", "9b" = "2011",
               "10a" = "2012", "10b" = "2011",
               "11a" = "2013", "11b" = "2012",
               "12a" = "2014", "12b" = "2013",
               "13a" = "2015", "13b" = "2014",
               "14a" = "2016", "14b" = "2015",
               "15a" = "2017", "15b" = "2016",
               "16a" = "2018", "16b" = "2017",
               "17a" = "2019", "17b" = "2018",
               "18a" = "2020", "18b" = "2019",
               "19a" = "2021", "19b" = "2020",
               "20a" = "2022", "20b" = "2021")

phi.CNI <- phi.CNI  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.CNI) ~ names.CNI[name],
    TRUE ~ name
  ))



phi.CNI$site <- "Nicaragua"

phi.CNI$name <- as.factor(phi.CNI$name)

label_positions <- seq(1, length(levels(phi.CNI$name)), by = 2)

si <- ggplot(phi.CNI, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Nicaragua")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.CNI$name),
    labels = ifelse(seq_along(levels(phi.CNI$name)) %in% label_positions, levels(phi.CNI$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.CNI$name)) %in% label_positions, 1, 0.5))
  )


## Manitoba (SW)
phi.VMA <- as.data.frame(GLMultiSiteNo2year$sims.list$phi.VMA)
colnames(phi.VMA) <- c("1a", "1b", "2a", "2b")


phi.VMA <- phi.VMA  %>%
  pivot_longer(1:4)

names.VMA <- c("1a" = "2009", "1b" = "2009",
               "2a" = "2010", "2b" = "2010")

phi.VMA <- phi.VMA  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.VMA) ~ names.VMA[name],
    TRUE ~ name
  ))



phi.VMA$site <- "Manitoba (SW)"

phi.VMA$name <- as.factor(phi.VMA$name)

label_positions <- seq(1, length(levels(phi.VMA$name)), by = 2)

sj <- ggplot(phi.VMA, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Manitoba (SW)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.VMA$name),
    labels = ifelse(seq_along(levels(phi.VMA$name)) %in% label_positions, levels(phi.VMA$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.VMA$name)) %in% label_positions, 1, 0.5))
  )


ggplot(phi.VMA, aes(x = name, y = value)) +
  geom_smooth()+
  ylim(0,1)


## Ontario (W)
phi.VONW <- as.data.frame(GLMultiSiteNo2year$sims.list$phi.VONW)
colnames(phi.VONW) <- c("1a", "2a")


phi.VONW <- phi.VONW  %>%
  pivot_longer(1:2)

names.VONW <- c("1a" = "2009", "2a" = "2010")

phi.VONW <- phi.VONW  %>%
  mutate(sex = ifelse(grepl("a", name), "M",
                      ifelse(grepl("b", name), "F", NA))) %>%
  mutate(name = case_when(
    name %in% names(names.VONW) ~ names.VONW[name],
    TRUE ~ name
  ))



phi.VONW$site <- "Ontario (W)"

phi.VONW$name <- as.factor(phi.VONW$name)

label_positions <- seq(1, length(levels(phi.VONW$name)), by = 2)

sk <- ggplot(phi.VONW, aes(x = name, y = value)) +
  geom_boxplot()+
  theme(panel.grid.major  = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title=element_text(colour="white"))+
  xlab("Year")+
  ylab("Apparent Annual Survival")+
  ggtitle("Ontario (W)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(
    breaks = levels(phi.VONW$name),
    labels = ifelse(seq_along(levels(phi.VONW$name)) %in% label_positions, levels(phi.VONW$name), "")
  )+
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Adjust text size, angle, and position
    axis.ticks.x = element_line(size = ifelse(seq_along(levels(phi.VONW$name)) %in% label_positions, 1, 0.5))
  )


################################################################################
###                       COMBINE ALL PLOTS FOR FINAL PLOT                   ###
################################################################################

ag <- ggarrange(sb, sc, se, sf, sg, sh, sj, sk,
                ncol = 3, nrow = 3)

ag1 <- ggarrange(ag, sa, sd, si,
                 ncol = 1,
                 heights = c(2.5, 1, 1, 1),
                 labels = c("", "I", "J", "K"))

## Save output
setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Tables and Figures/Figure A3")
ggsave("ag1.png", dpi = 1200, width = 8.5, height = 15, units = "in")
