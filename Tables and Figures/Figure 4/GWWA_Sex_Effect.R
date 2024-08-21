## Figure 4 in manuscript. This visualization shows sex-specific mean apparent
## annual survival for the Great Lakes and the Appalachian populations. Boxes 
## represent the interquartile range.

## Load packages
library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)

## Create dataframe for Appalachian sex-specific survival
app.f.survival <- as.data.frame(AppsMultiSiteNo2year$sims.list$overall.female.survival)
colnames(app.f.survival) <- "App Female"
app.m.survival <- as.data.frame(AppsMultiSiteNo2year$sims.list$overall.male.survival)
colnames(app.f.survival) <- "App Male"

## Create dataframe for Great Lakes sex-specific survival
gl.f.survival <- as.data.frame(GLMultiSiteNo2year$sims.list$overall.female.survival.gl)
colnames(app.f.survival) <- "GL Female"
gl.m.survival <- as.data.frame(GLMultiSiteNo2year$sims.list$overall.male.survival.gl)
colnames(app.f.survival) <- "GL Male"

## Combine dataframes and change column names
sex_df <- cbind(gl.m.survival, gl.f.survival, app.m.survival, app.f.survival)
colnames(sex_df) <- c("GL Male", "GL Female", "App Male", "App Female")

## Pivot dataframe
sex_df <- sex_df %>%
  pivot_longer(cols = starts_with("GL") | starts_with("App"),
               names_to = c("Region", "Sex"),
               names_pattern = "^(\\w+)\\s(\\w+)$") %>%
  mutate(Region = case_when(
    Region == "GL" ~ "Great Lakes",
    Region == "App" ~ "Appalachian",
    TRUE ~ as.character(Region) 
  )) %>%
  mutate(Region = factor(Region, levels = c("Great Lakes", "Appalachian"))) %>%  # Specify the order of regions
  arrange(Region, Sex)

## Create visual
sexef <- ggplot(sex_df, aes(x = Region, y = value, fill = Sex)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Female" = "white", "Male" = "gray")) +
  theme(panel.grid.major = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.title = element_text(colour = "transparent"),
        legend.background = element_rect(fill = "transparent"), # Make legend background transparent
        legend.key = element_rect(fill = "transparent", colour = NA)) +
  ylab("Annual Apparent Survival") +
  xlab("Region") +
  ylim(0.25, 0.7) +
  theme(
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent"),
    panel.grid.major = element_line(color = "transparent"),
    panel.grid.minor = element_line(color = "transparent")
  )

## Save output
setwd("/Users/emilyfiliberti/Desktop/OneDrive - University of Maine System/Manuscripts/Historical GWWA Survival/Tables and Figures/Figure 4")

ggsave("sexef.png", dpi = 1200, width = 4, height = 4, units = "in")
