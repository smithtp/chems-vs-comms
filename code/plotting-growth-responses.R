#################################################################
# Plots to visualise growth responses
# of chemical impacted communities
#
# T. Smith 2025
#################################################################

# load some packages
library(tidyverse)

library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggbeeswarm)

setwd("~/Documents/chems-vs-comms/code/")

### --- Load the data --- ###

tidy_data <- read.csv("../data/tidy-data.csv") # these are the raw growth curves

################

AUC_data <- read.csv("../data/spline_fits.csv") # these are the AUCs from spline curve fits

# calculate dAUC 
control_data <- AUC_data %>%
  filter(Complexity == 0) %>%
  group_by(Community, Rep) %>% # unique plate for biological replicate
  summarise(Control.AUC = mean(AUC))

# add controls to the data and calculate change
AUC_data <- AUC_data %>%
  filter(Complexity > 0) %>%
  left_join(control_data, by = c("Community", "Rep")) %>%
  mutate(dAUC = AUC/Control.AUC)


# tidy up the community names in both datasets:
comm.numbers <- data.frame(Community = c("R_ID_7.16_10_R10(1)_TSS",
                                         "R_ID_7.16_11_R11(1)_TSS",
                                         "R_ID_7.16_12_R5(1)_TSS",
                                         "R_ID_7.16_1_R10(1)_TSS",
                                         "R_ID_7.16_4_R10(1)_TSS",
                                         "R_ID_7.16_5_R10(1)_TSS",
                                         "R_ID_7.16_6_R5(1)_TSS",
                                         "R_ID_7.16_7_R1(1)_TSS",
                                         "R_ID_7.16_8_R5(1)_TSS",
                                         "R_ID_7.16_9_R1(1)_TSS"),
                           CommunityNumber = c(10, 11, 12, 1, 4, 5, 6, 7, 8, 9))

tidy_data <- left_join(tidy_data, comm.numbers)

AUC_data <- left_join(AUC_data, comm.numbers)

# whats going on in Community 11, 19F5 and Comm 5 1C2?
ggplot(tidy_data %>% filter(CommunityNumber == 11, PlateWell == "4F5"), 
       aes(x = Time, y = OD600, Group = TruePlateWell)) + geom_line()
ggplot(tidy_data %>% filter(CommunityNumber == 5, PlateWell == "1C2"), 
       aes(x = Time, y = OD600, Group = TruePlateWell)) + geom_line()
# something went wrong with these replicates - aren't representative, remove from downstream analysis

tidy_data <- tidy_data %>% filter(!(CommunityNumber == 11 & TruePlateWell == "19F5"),
                                  !(CommunityNumber == 5 & TruePlateWell == "1C2"))
AUC_data <- AUC_data %>% filter(!(CommunityNumber == 11 & TruePlateWell == "1C2"),
                                !(CommunityNumber == 5 & TruePlateWell == "1C2"))

### --- Plotting! --- ###

png("../results/growth-curve-complexity.png", width = 2000, height = 300)
ggplot(tidy_data, aes(x = Time, y = OD600, group = TruePlateWell)) +
  geom_line(aes(col = Complexity)) +
  coord_cartesian(expand = FALSE, xlim = c(0, 72)) + # prevent plot area being wider than the axes
  labs(x = "Time (Hours)") +
  scale_colour_gradient(low = "black", high = "red") +
  facet_wrap(~CommunityNumber, nrow = 1) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.03, 0.7),
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))
dev.off()

png("../results/growth-curve-oxytet.png", width = 2000, height = 300)
ggplot(tidy_data, aes(x = Time, y = OD600, group = TruePlateWell)) +
  geom_line(aes(col = as.factor(Oxytetracycline))) +
  coord_cartesian(expand = FALSE, xlim = c(0, 72)) + # prevent plot area being wider than the axes
  labs(x = "Time (Hours)") +
  scale_colour_manual(values = c("black", "orange")) +
  #scale_colour_gradient(low = "black", high = "red") +
  facet_wrap(~CommunityNumber, nrow = 1) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.03, 0.7),
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_blank())
dev.off()

# this looks like an "ecological suicide" effect, where a non-lethal amount of antibiotic
# has actually improved growth by the end of the experiment in the majority of cases.

png("../results/dAUC-oxytet.png", width = 2000, height = 300)
ggplot(AUC_data, aes(x = Complexity, y = dAUC, fill = as.factor(Oxytetracycline))) +
  geom_beeswarm(shape = 21, cex = 0.5, alpha = 0.6) +
  scale_fill_manual(values = c("black", "orange")) +
  geom_smooth(method = lm, col = "black") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~CommunityNumber, nrow = 1) +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.12, 0.8),
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_blank())
dev.off()

# represent this in a linear mixed model
# Q: does dAUC go down with increasing chemicals, if we account for the oxytet intercept
# and random effect of community?

dauc_model <- lmer(dAUC ~ Complexity + as.factor(Oxytetracycline) + (1|Community), data = AUC_data)
summary(dauc_model)

# general negative effect of the oxytetracycline
# and then a very small positive effect of increasing chemical complexity