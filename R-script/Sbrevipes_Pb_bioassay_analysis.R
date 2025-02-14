# Load libraries
library(ggplot2)
library(dplyr)
library(scales)
library(rstatix)
library(tidyr)
library(multcompView)

# Load data
Pbdata <- read.csv("~/Desktop/Branco Lab/Results/Suillus_Pb_Bioassays/Sbrevipes/Sbrevipes_Pb_bioassay_data.csv", header = TRUE)
Pbdata$Treatment <- as.character(Pbdata$Treatment)

# "Template" to put data back into
PbdataFrame <- Pbdata %>%
  select(Seedling, Fungus, Treatment)

###
###
###

# Plant Length analysis

# Select data for root/shoot length graph
shootdata <- Pbdata %>%
  select(Seedling, Shoot_length_cm)
rootdata <- Pbdata %>%
  select(Seedling, Root_length_cm)

# rm outliers
shoot_outliers <- Pbdata %>%
  group_by(Fungus, Treatment) %>%
  identify_outliers(Shoot_length_cm)

root_outliers <- Pbdata %>%
  group_by(Fungus, Treatment) %>%
  identify_outliers(Root_length_cm)

shootdata <- shootdata[!(shootdata$Seedling %in% shoot_outliers$Seedling),]
rootdata <- rootdata[!(rootdata$Seedling %in% root_outliers$Seedling),]

# double check for outliers
shoot_outliers2 <- shootdata %>%
  identify_outliers(Shoot_length_cm)

root_outliers2 <- rootdata %>%
  identify_outliers(Root_length_cm)

shootdata <- shootdata[!(shootdata$Seedling %in% shoot_outliers2$Seedling),]
rootdata <- rootdata[!(rootdata$Seedling %in% root_outliers2$Seedling),]

# repeat until no outliers
# still outliers in shoot
shoot_outliers3 <- shootdata %>%
  identify_outliers(Shoot_length_cm)

shootdata <- shootdata[!(shootdata$Seedling %in% shoot_outliers3$Seedling),]

# no more outliers present

# merge back together
merge1 <- merge(PbdataFrame, shootdata, by = "Seedling", all = T)
merge1
merge2 <- merge(merge1, rootdata, by = "Seedling", all = T)
merge2


# Determine assumptions met for stats (SHOOT LENGTH)
model <- lm(Shoot_length_cm ~ interaction(Fungus, Treatment), data = merge2)
shapiro_test(residuals(model))

merge2 %>% levene_test(Shoot_length_cm ~ Fungus * Treatment)

# p values for both Shapiro (normality) and Levene (homogeneity of variance) tests should be ns

# perform stats (SHOOT LENGTH)
anova <- aov(Shoot_length_cm ~ Fungus * Treatment, data = merge2)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
multcompLetters4(anova, tukey)

# Determine assumptions met for stats (ROOT LENGTH)
model <- lm(Root_length_cm ~ interaction(Fungus, Treatment), data = merge2)
shapiro_test(residuals(model))

merge2 %>% levene_test(Root_length_cm ~ Fungus * Treatment)

# perform stats (ROOT LENGTH)
anova <- aov(Root_length_cm ~ Fungus * Treatment, data = merge2)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
multcompLetters4(anova, tukey)

# create df for graphing
data_long <- merge2 %>%
  pivot_longer(cols = c(Shoot_length_cm, Root_length_cm),
               names_to = "Type",
               values_to = "Length")
data_long <- na.omit(data_long)

# Summarize data for mean and standard deviation
data_summary <- data_long %>%
  group_by(Fungus, Treatment, Type) %>%
  summarise(
    Mean = mean(Length),
    SD = sd(Length),
    .groups = "drop"
  )
data_summary

# mutate the data so that root goes below x axis
# so that the SD remains positice 
# and make "fillgroup" so we can specify colours later
data_summary <- data_summary %>%
  mutate(
    StackedOffset = ifelse(Type == "Root_length_cm", -Mean, Mean),
    SD_Pos = SD,
    FillGroup = paste(Fungus, Type, sep = "_")
  )
data_summary

# Plot
ggplot(data_summary, aes(
  x = Treatment, y = StackedOffset,
  fill = FillGroup, group = Fungus)) +
  geom_bar(
    stat = "identity", position = position_dodge(),
    width = 0.5, colour = "black") +
  geom_errorbar(
    aes(
      ymin = StackedOffset - SD_Pos,
      ymax = StackedOffset + SD_Pos),
    position = position_dodge(0.5), width = 0.2, colour = "black") +
  scale_y_continuous(
    name = "Root & Shoot Length (cm)",
    breaks = scales::breaks_width(1),
    limits = c(),
    labels = ~scales::label_number(accuracy = 1)(abs(.x)))+
  scale_fill_manual(
    values = c(
      "S.b. 134_Shoot_length_cm" = "blue", 
      "S.b. 134_Root_length_cm" = "red", 
      "none_Shoot_length_cm" = "darkblue", 
      "none_Root_length_cm" = "darkred" 
    )
  ) +
  labs(
    x = "Treatment",
    y = "Length",
    fill = "Sample & Type",
    title = "Shoot and Root Lengths by Treatment and Sample"
  )

# select data for TOTAL LENGTH
Pbdata$totalLength <- Pbdata$Shoot_length_cm + Pbdata$Root_length_cm

totalLength <- Pbdata %>%
  select(Seedling, totalLength)

# rm outliers
totalLength_outliers <- Pbdata %>%
  group_by(Fungus, Treatment) %>%
  identify_outliers(totalLength)

totalLength <- totalLength[!(totalLength$Seedling %in% totalLength_outliers$Seedling),]

# repeat until no outliers
totalLength_outliers2 <- totalLength %>%
  identify_outliers(totalLength)

# no more outliers present

merge1 <- merge(PbdataFrame, totalLength, by = "Seedling", all = T)
merge1 <- na.omit(merge1)


# Determine assumptions met for stats (TOTAL LENGTH)
model <- lm(totalLength ~ interaction(Fungus, Treatment), data = merge1)
shapiro_test(residuals(model))

merge1 %>% levene_test(totalLength ~ Fungus * Treatment)

# perform stats (TOTAL LENGTH)
anova <- aov(totalLength ~ Fungus * Treatment, data = merge1)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
multcompLetters4(anova, tukey)


####
####
####

# Biomass - wet

# Select data for root/shoot biomass graph
shootWBdata <- Pbdata %>%
  select(Seedling, Shoot_biomass_wet_mg)
rootWBdata <- Pbdata %>%
  select(Seedling, Root_biomass_wet_mg)

# rm outliers
shootWB_outliers <- Pbdata %>%
  group_by(Fungus, Treatment) %>%
  identify_outliers(Shoot_biomass_wet_mg)

rootWB_outliers <- Pbdata %>%
  group_by(Fungus, Treatment) %>%
  identify_outliers(Root_biomass_wet_mg)

shootWBdata <- shootWBdata[!(shootWBdata$Seedling %in% shootWB_outliers$Seedling),]
rootWBdata <- rootWBdata[!(rootWBdata$Seedling %in% rootWB_outliers$Seedling),]

# double check for outliers - repeat 
shootWB_outliers2 <- shootWBdata %>%
  identify_outliers(Shoot_biomass_wet_mg)

rootWB_outliers2 <- rootWBdata %>%
  identify_outliers(Root_biomass_wet_mg)

shootWBdata <- shootWBdata[!(shootWBdata$Seedling %in% shootWB_outliers2$Seedling),]
rootWBdata <- rootWBdata[!(rootWBdata$Seedling %in% rootWB_outliers2$Seedling),]

# no more shoot outliers, but one more root outlier
rootWB_outliers3 <- rootWBdata %>%
  identify_outliers(Root_biomass_wet_mg)
rootWBdata <- rootWBdata[!(rootWBdata$Seedling %in% rootWB_outliers3$Seedling),]

rootWB_outliers4 <- rootWBdata %>%
  identify_outliers(Root_biomass_wet_mg)
rootWBdata <- rootWBdata[!(rootWBdata$Seedling %in% rootWB_outliers4$Seedling),]

rootWB_outliers5 <- rootWBdata %>%
  identify_outliers(Root_biomass_wet_mg)
rootWBdata <- rootWBdata[!(rootWBdata$Seedling %in% rootWB_outliers5$Seedling),]

rootWB_outliers6 <- rootWBdata %>%
  identify_outliers(Root_biomass_wet_mg)
rootWBdata <- rootWBdata[!(rootWBdata$Seedling %in% rootWB_outliers6$Seedling),]

# no more outliers

# merge back together
merge1 <- merge(PbdataFrame, shootWBdata, by = "Seedling", all = T)
merge1
merge2 <- merge(merge1, rootWBdata, by = "Seedling", all = T)
merge2

# Determine assumptions met for stats (SHOOT WB)
model <- lm(Shoot_biomass_wet_mg ~ interaction(Fungus, Treatment), data = na.omit(merge2))
shapiro_test(residuals(model))

na.omit(merge2) %>% levene_test(Shoot_biomass_wet_mg ~ Fungus * Treatment)

# perform stats (SHOOT WB)
anova <- aov(Shoot_biomass_wet_mg ~ Fungus * Treatment, data = na.omit(merge2))
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
multcompLetters4(anova, tukey)

## Determine assumptions met for stats (ROOT WB)
model <- lm(Root_biomass_wet_mg ~ interaction(Fungus, Treatment), data = na.omit(merge2))
shapiro_test(residuals(model))

na.omit(merge2) %>% levene_test(Root_biomass_wet_mg ~ Fungus * Treatment)

# perform stats (ROOT WB)
anova <- aov(Root_biomass_wet_mg ~ Fungus * Treatment, data = na.omit(merge2))
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
multcompLetters4(anova, tukey)

# make df for plotting
data_long <- merge2 %>%
  pivot_longer(cols = c(Shoot_biomass_wet_mg, Root_biomass_wet_mg),
               names_to = "Type",
               values_to = "BW")
data_long <- na.omit(data_long)

# Summarize data for mean and standard deviation
data_summary <- data_long %>%
  group_by(Fungus, Treatment, Type) %>%
  summarise(
    Mean = mean(BW),
    SD = sd(BW),
    .groups = "drop"
  )
data_summary

data_summary <- data_summary %>%
  mutate(
    StackedOffset = ifelse(Type == "Root_biomass_wet_mg", -Mean, Mean), 
    SD_Pos = SD,
    FillGroup = paste(Fungus, Type, sep = "_")
  )
data_summary

# Plot
ggplot(data_summary, aes(
  x = Treatment, y = StackedOffset,
  fill = FillGroup, group = Fungus)) +
  geom_bar(
    stat = "identity", position = position_dodge(),
    width = 0.5, colour = "black") +
  geom_errorbar(
    aes(
      ymin = StackedOffset - SD_Pos,
      ymax = StackedOffset + SD_Pos),
    position = position_dodge(0.5), width = 0.2, colour = "black") +
  scale_y_continuous(
    name = "Wet Biomass (mg)",
    breaks = scales::breaks_width(10),
    limits = c(),
    labels = ~scales::label_number(accuracy = 1)(abs(.x)))+
  scale_fill_manual(
    values = c(
      "S.b. 134_Shoot_biomass_wet_mg" = "blue",
      "S.b. 134_Root_biomass_wet_mg" = "red",
      "none_Shoot_biomass_wet_mg" = "darkblue",
      "none_Root_biomass_wet_mg" = "darkred"
    )
  ) +
  labs(
    x = "Treatment",
    y = "Length",
    fill = "Sample & Type",
    title = "Shoot and Root Biomass by Treatment and Sample"
  )

## total wet biomass

# select data for TOTAL BIOMASS WET
Pbdata$totalBiomassWet <- Pbdata$Shoot_biomass_wet_mg + Pbdata$Root_biomass_wet_mg

totalBiomassWet <- Pbdata %>%
  select(Seedling, totalBiomassWet)

totalBiomassWet_outliers <- Pbdata %>%
  group_by(Fungus, Treatment) %>%
  identify_outliers(totalBiomassWet)

# rm outliers
totalBiomassWet <- totalBiomassWet[!(totalBiomassWet$Seedling %in% totalBiomassWet_outliers$Seedling),]

# Double check for outliers
totalBiomassWet_outliers2 <- totalBiomassWet %>%
  identify_outliers(totalBiomassWet)

totalBiomassWet <- totalBiomassWet[!(totalBiomassWet$Seedling %in% totalBiomassWet_outliers2$Seedling),]


merge1 <- merge(PbdataFrame, totalBiomassWet, by = "Seedling", all = T)
merge1 <- na.omit(merge1)

## Determine assumptions met for stats (TOTAL BIOMASS WET)
model <- lm(totalBiomassWet ~ interaction(Fungus, Treatment), data = merge1)
shapiro_test(residuals(model))

merge1 %>% levene_test(totalBiomassWet ~ Fungus * Treatment)

# perform stats (TOTAL BIOMASS WET)
anova <- aov(totalBiomassWet ~ Fungus * Treatment, data = merge1)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
multcompLetters4(anova, tukey)

###
###
###

# biomass - DRY

# Select data for root/shoot dry biomass graph
shootDBdata <- Pbdata %>%
  select(Seedling, Shoot_biomass_dry_mg)
rootDBdata <- Pbdata %>%
  select(Seedling, Root_biomass_dry_mg)

# rm outliers
shootDB_outliers <- Pbdata %>%
  group_by(Fungus, Treatment) %>%
  identify_outliers(Shoot_biomass_dry_mg)

rootDB_outliers <- Pbdata %>%
  group_by(Fungus, Treatment) %>%
  identify_outliers(Root_biomass_dry_mg)

shootDBdata <- shootDBdata[!(shootDBdata$Seedling %in% shootDB_outliers$Seedling),]
rootDBdata <- rootDBdata[!(rootDBdata$Seedling %in% rootDB_outliers$Seedling),]

# double check for outliers - repeat 
shootDB_outliers2 <- shootDBdata %>%
  identify_outliers(Shoot_biomass_dry_mg)

rootDB_outliers2 <- rootDBdata %>%
  identify_outliers(Root_biomass_dry_mg)

shootDBdata <- shootDBdata[!(shootDBdata$Seedling %in% shootDB_outliers2$Seedling),]
rootDBdata <- rootDBdata[!(rootDBdata$Seedling %in% rootDB_outliers2$Seedling),]

# no more shoot outliers
# repeat root outliers check
rootDB_outliers3 <- rootDBdata %>%
  identify_outliers(Root_biomass_dry_mg)
rootDBdata <- rootDBdata[!(rootDBdata$Seedling %in% rootDB_outliers3$Seedling),]

# still one root outlier
rootDB_outliers4 <- rootDBdata %>%
  identify_outliers(Root_biomass_dry_mg)
rootDBdata <- rootDBdata[!(rootDBdata$Seedling %in% rootDB_outliers4$Seedling),]

merge1 <- merge(PbdataFrame, shootDBdata, by = "Seedling", all = T)
merge1
merge2 <- merge(merge1, rootDBdata, by = "Seedling", all = T)
merge2

## Determine assumptions met for stats (SHOOT DB)
model <- lm(Shoot_biomass_dry_mg ~ interaction(Fungus, Treatment), data = na.omit(merge2))
shapiro_test(residuals(model))

na.omit(merge2) %>% levene_test(Shoot_biomass_dry_mg ~ Fungus * Treatment)

# perform stats (SHOOT DB)
anova <- aov(Shoot_biomass_dry_mg ~ Fungus * Treatment, data = na.omit(merge2))
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
multcompLetters4(anova, tukey)

## Determine assumptions met for stats (ROOT DB)
model <- lm(Root_biomass_dry_mg ~ interaction(Fungus, Treatment), data = na.omit(merge2))
shapiro_test(residuals(model))

na.omit(merge2) %>% levene_test(Root_biomass_dry_mg ~ Fungus * Treatment)

# perform stats (ROOT DB)
anova <- aov(Root_biomass_dry_mg ~ Fungus * Treatment, data = na.omit(merge2))
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
multcompLetters4(anova, tukey)


## 
data_long <- merge2 %>%
  pivot_longer(cols = c(Shoot_biomass_dry_mg, Root_biomass_dry_mg),
               names_to = "Type",
               values_to = "BD")
data_long <- na.omit(data_long)

# Summarize data for mean and standard deviation
data_summary <- data_long %>%
  group_by(Fungus, Treatment, Type) %>%
  summarise(
    Mean = mean(BD),
    SD = sd(BD),
    .groups = "drop"
  )
data_summary

data_summary <- data_summary %>%
  mutate(
    StackedOffset = ifelse(Type == "Root_biomass_dry_mg", -Mean, Mean), 
    SD_Pos = SD,
    FillGroup = paste(Fungus, Type, sep = "_")
  )
data_summary

# Plot
ggplot(data_summary, aes(
  x = Treatment, y = StackedOffset,
  fill = FillGroup, group = Fungus)) +
  geom_bar(
    stat = "identity", position = position_dodge(),
    width = 0.5, colour = "black") +
  geom_errorbar(
    aes(
      ymin = StackedOffset - SD_Pos,
      ymax = StackedOffset + SD_Pos),
    position = position_dodge(0.5), width = 0.2, colour = "black") +
  scale_y_continuous(
    name = "Dry Biomass (mg)",
    breaks = scales::breaks_width(10),
    limits = c(),
    labels = ~scales::label_number(accuracy = 1)(abs(.x)))+
  scale_fill_manual(
    values = c(
      "S.b. 134_Shoot_biomass_dry_mg" = "blue",
      "S.b. 134_Root_biomass_dry_mg" = "red",
      "none_Shoot_biomass_dry_mg" = "darkblue",
      "none_Root_biomass_dry_mg" = "darkred" 
    )
  ) +
  labs(
    x = "Treatment",
    y = "Length",
    fill = "Sample & Type",
    title = "Shoot and Root Biomass by Treatment and Sample"
  )

## total dry biomass

# select data for TOTAL BIOMASS DRY
Pbdata$totalBiomassDry <- Pbdata$Shoot_biomass_dry_mg + Pbdata$Root_biomass_dry_mg

totalBiomassDry <- Pbdata %>%
  select(Seedling, totalBiomassDry)

totalBiomassDry_outliers <- Pbdata %>%
  group_by(Fungus, Treatment) %>%
  identify_outliers(totalBiomassDry)

# rm outliers
totalBiomassDry <- totalBiomassDry[!(totalBiomassDry$Seedling %in% totalBiomassDry_outliers$Seedling),]

# Double check for outliers
totalBiomassDry_outliers2 <- totalBiomassDry %>%
  identify_outliers(totalBiomassDry)

totalBiomassDry <- totalBiomassDry[!(totalBiomassDry$Seedling %in% totalBiomassDry_outliers2$Seedling),]

#repeat until no outliers
totalBiomassDry_outliers3 <- totalBiomassDry %>%
  identify_outliers(totalBiomassDry)

totalBiomassDry <- totalBiomassDry[!(totalBiomassDry$Seedling %in% totalBiomassDry_outliers3$Seedling),]


merge1 <- merge(PbdataFrame, totalBiomassDry, by = "Seedling", all = T)
merge1 <- na.omit(merge1)

## Determine assumptions met for stats (TOTAL BIOMASS DRY)
model <- lm(totalBiomassDry ~ interaction(Fungus, Treatment), data = merge1)
shapiro_test(residuals(model))

merge1 %>% levene_test(totalBiomassDry ~ Fungus * Treatment)

# perform stats (TOTAL BIOMASS DRY)
anova <- aov(totalBiomassDry ~ Fungus * Treatment, data = merge1)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
multcompLetters4(anova, tukey)

###
###
###

## ROOT TO SHOOT RATIO

Pbdata$root_shoot_ratio <- Pbdata$Root_biomass_dry_mg/(Pbdata$Shoot_biomass_dry_mg)

ratio <- Pbdata %>%
  select(Seedling, Fungus, Treatment, root_shoot_ratio)

ratio_outliers <- ratio %>%
  identify_outliers(root_shoot_ratio)
ratio <- ratio[!(ratio$Seedling %in% ratio_outliers$Seedling),]

# check again for outliers
ratio_outliers2 <- ratio %>%
  identify_outliers(root_shoot_ratio)
ratio <- ratio[!(ratio$Seedling %in% ratio_outliers2$Seedling),]

# and again
ratio_outliers3 <- ratio %>%
  identify_outliers(root_shoot_ratio)
ratio <- ratio[!(ratio$Seedling %in% ratio_outliers3$Seedling),]

ratio <- na.omit(ratio)

# Determine assumptions met for stats
model <- lm(root_shoot_ratio ~ interaction(Fungus, Treatment), data = ratio)
shapiro_test(residuals(model))

ratio %>% levene_test(root_shoot_ratio ~ Fungus * Treatment)

# perform stats (ratio)
anova <- aov(root_shoot_ratio ~ Fungus * Treatment, data = ratio)
summary(anova)
tukey <- TukeyHSD(anova)
print(tukey)
multcompLetters4(anova, tukey)

# summary to graph 
data_summary <- ratio %>% group_by(Fungus, Treatment) %>%
  summarise(mean=mean(root_shoot_ratio), sd=sd(root_shoot_ratio)) %>%
  arrange(desc(mean))
print(data_summary)

# plot
ggplot(data_summary, aes(
  x = Treatment, y = mean,
  fill = Fungus, group = Fungus)) +
  geom_bar(
    stat = "identity", position = position_dodge(),
    width = 0.5, colour = "black") +
  geom_errorbar(
    aes(
      ymin=mean-sd, ymax=mean+sd),
    position = position_dodge(0.5), width = 0.2, colour = "black") +
  scale_y_continuous(
    name = "Root:Shoot Ratio",
    limits = c(),
  ) +
  scale_fill_manual(
    values = c("#663399", "#9966CC")) +
  labs(
    x = "Treatment",
    title = "Root:Shoot Ratio"
  )


###
###
###

# Colonization rate

colonization <- Pbdata[!Pbdata$Root_colinization_percent=="n/a",]
colonization <- colonization[!colonization$Fungus=="none",]

colonization %>%
  group_by(Fungus, Treatment) %>%
  get_summary_stats(Root_colinization_percent, type = "mean_sd")

colonization %>%
  group_by(Fungus, Treatment) %>%
  identify_outliers(Root_colinization_percent) 
# there are no outliers

# confirm assumptions for stats
model_col <- lm(Root_colinization_percent ~ interaction(Fungus, Treatment), data = colonization)
shapiro_test(residuals(model_col))

colonization %>% levene_test(Root_colinization_percent ~ Fungus * Treatment)

# two sample t-test as only 2 groups to be compared
t.test(Root_colinization_percent ~ Treatment, data = colonization)

# summary for graphing
data_summary <- colonization %>% group_by(Treatment) %>%
  summarise(mean=mean(Root_colinization_percent), sd=sd(Root_colinization_percent)) %>%
  arrange(desc(mean))
print(data_summary)

# plot
ggplot(data_summary, aes(
  x = Treatment, y = mean,
  fill = Treatment, group = Treatment)) +
  geom_bar(
    stat = "identity", position = position_dodge(),
    width = 0.5, colour = "black") +
  geom_errorbar(
    aes(
      ymin=mean-sd, ymax=mean+sd),
    position = position_dodge(0.5), width = 0.2, colour = "black") +
  scale_y_continuous(
    name = "Percent Colonization (%)",
    breaks = scales::breaks_width(5),
    limits = c(),
  ) +
  scale_fill_manual(values = c("#9966CC", "#9966CC")) +
  scale_colour_manual(values = c("#9966CC", "#9966CC")) +  
  labs(
    x = "Treatment",
    title = "Fungal Root Colonization"
  )
