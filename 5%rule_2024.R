###############################################################

# Re-weighing the 5% tagging recommendation: assessing the potential impacts of tags on the behavior and body condition of bats

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R studio (v. 4.2.2)
# # Additional graphic edits were performed in Inkscape

###############################################################

# setwd

# clean the workspace -----------------------------------------------------
rm(list = ls())

# Loading R package -------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(irr, lpSolve, metafor, ggplot2, tidyverse, scatterpie, ggpubr, dplyr, grid, gridExtra, gridtext, maps, svglite, cowplot)

# Loading the Databases ----------------------------------------------------
cohens <- read.csv("cohens.csv")
meta <- read.csv(file = "Meta5Final.csv")
data <- read.csv(file = "new2022_1.csv", stringsAsFactors = T)

# Set custom colors for plotting
my_colors <- c("#B46AAA","#401240")
my_colors2 <- c("#B46AAA","#401240", "#E0E0E0")

## Data preparation for cohens:

###############################################################

# Cohen's Kappa for inter-rater reliability ----------------------
kappa_result <- kappa2(cohens[,3:4])
print(kappa_result)

## Data preparation for meta-analysis and plot:

###############################################################

# Meta-analysis ----------------------------------------------------
meta$Pearson.s_r <- as.numeric(meta$Pearson.s_r)
meta$N <- as.numeric(meta$N)
meta <- meta %>% select(ID, N, Family, Genus, Response_Group, r = Pearson.s_r)

# Derive Fischer's Z and its variance
meta <- metafor::escalc(measure = "COR", ri = r, ni = N, data = meta)
meta <- na.omit(meta)

# meta-analysis for effects of Tag on "health" "behavior"
Hmeta <- rma.mv(yi, vi, random =  ~ 1 | ID, data = meta[meta$Response_Group == "Health",] )
Bmeta <- rma.mv(yi, vi, random =  ~ 1 | ID, data = meta[meta$Response_Group == "Behavior",] ) 

# Evaluation of publication bias via Rosenthal’s method.
failsafe_H <- fsn(yi, vi, type = "Rosenthal", data = meta[meta$Response_Group == "Health",]) 
failsafe_B <- fsn(yi, vi, type = "Rosenthal", data = meta[meta$Response_Group == "Behavior",])


# Extracting estimates for body condition plot
plotHealth <- data.frame(label = c("Pup mass", "Pup condition", "% body condition", "Body mass (f)", "Body mass (e)", "Body mass (d)", "Body mass (c)", 
                                   "Body mass (b)", "Body mass (a)", "BCI", "Body condition"),
                         
                         ES    = c( 0.040, -0.719, -0.401, -0.825, -0.192, -0.001, -0.192, -0.558, -0.015, -0.004, ((exp(Hmeta$b)-1))/((exp(Hmeta$b)+1)) ),
                         L     = c( NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, ((exp(Hmeta$ci.lb)-1)/(exp(Hmeta$ci.lb)+1)) ),
                         U     = c( NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, ((exp(Hmeta$ci.ub)-1)/(exp(Hmeta$ci.ub)+1)) )
                         
)


plotHealth$label <- factor(as.character(plotHealth$label), levels = unique(plotHealth$label))


# Plot for body condition
h <- ggplot(plotHealth, aes(y=label, x=ES, xmin=L, xmax=U)) + 
  geom_point(color = '#0B3C49', shape=16, size=2.5) + 
  geom_errorbarh(height=0, size= 0.5, color = '#0B3C49') +
  geom_vline(xintercept=0, color='black', linetype='dashed') +  
  scale_x_continuous(limits=c(-1,1), name=("")) +
  ylab('') +
  labs(subtitle="Body condition") +
  theme(text = element_text(size = 15),
        axis.text.y = element_text(face = c(rep('plain', 10), 'bold')),
        plot.margin = unit(c(0.5,0.2,0,0.85), 'lines'), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Extracting estimates for behavior plot
plotBehavior <- data.frame(label = c("Roost use", "Night activity","Maneuverability", "Foraging time", "Flight duration", "Flight distance", "Emergence time", "Behavior"),
                         
                         ES    = c( -0.11, -0.284, -0.857, -0.356, 0.046, -0.10, -0.740, ((exp(Bmeta$b)-1))/((exp(Bmeta$b)+1)) ) ,
                         L     = c( NA, NA, NA, NA, NA, NA, NA, ((exp(Bmeta$ci.lb)-1)/(exp(Bmeta$ci.lb)+1)) ),
                         U     = c( NA, NA, NA, NA, NA, NA, NA, ((exp(Bmeta$ci.ub)-1)/(exp(Bmeta$ci.ub)+1)) )
)


plotBehavior$label <- factor(as.character(plotBehavior$label), levels = unique(plotBehavior$label))


# Plot for behavior
b <- ggplot(plotBehavior, aes(y=label, x=ES, xmin=L, xmax=U)) + 
  geom_point(color = "#0B3C49", shape=16, size=2.5) + 
  geom_errorbarh(height=0, size= 0.5, color = "#0B3C49") +
  geom_vline(xintercept=0, color='black', linetype='dashed') + 
  scale_x_continuous(limits=c(-1,1), name=("")) +
  ylab('') +
  labs(subtitle="Behavior") +
  theme(text = element_text(size = 15),
        axis.text.y = element_text(face = c(rep('plain', 7), 'bold')),
        plot.margin = unit(c(0.5,0.2,0,0.85), 'lines'), 
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# View and save
Meta <- grid.arrange(b, h, ncol = 2)
ggsave("Meta.svg", Meta, dpi = 300, width = 250, height = 100, units = "mm")


## qnorm (calculating Z value) for Wilcoxon Meta analyses

###############################################################
# these data were entered into the meta-analysis csv file
# don't run

#WoS_0014
qnorm(0.001/2)
#Other_other_001
qnorm(0.84/2)
#Other_other_002
qnorm(0.9/2) 
#Other_other_004
qnorm(0.14/2) 



## Data preparation for Regression and Bar plot

###############################################################

# Loading the Database 
trend <- data %>% select(ID, Year_publication, Higher_geography, H_5)

# retain distinct papers based on ID
trend <- distinct(trend, ID, .keep_all = TRUE) 

# remove NA values
trend <- na.omit(trend)

# just for model
regress <- trend[trend$Year_publication > 1998 & trend$Year_publication < 2023, ]

# Overall temporal trend in exceeding the 5% guideline
regress <- data.frame(table(regress$Year_publication, regress$H_5)) ; colnames(regress) <- c("yr","Higher than 5","N")
regress$yr <- as.numeric(as.character(regress$yr))

trend <- data.frame(table(trend$Year_publication, trend$H_5)) ; colnames(trend) <- c("yr","Higher than 5","N")
trend$yr <- as.numeric(as.character(trend$yr))

# Has the frequency of exceeding the 5% rule changed over time?
glm_data <- data.frame(yr = unique(regress$yr),
                       exceed = regress[regress$`Higher than 5`=="Yes",]$N, 
                       under = regress[regress$`Higher than 5`=="No",]$N)

# create a data frame for bar plot
trend_yr <- data.frame(yr = unique(trend$yr),
                       exceed = trend[trend$`Higher than 5`=="Yes",]$N, 
                       under = trend[trend$`Higher than 5`=="No",]$N)

# summary stats
sum(table(trend_yr$exceed))

# Fit the model
m1 <- glm(cbind(exceed, under) ~ yr, data = glm_data, family = "binomial")
performance::check_overdispersion(m1)
rsq::rsq(m1)

# Convert yr to numeric
trend_yr$yr <- as.numeric(as.character(trend_yr$yr))

# Reshape data frame to long format for ggplot2
trend_long <- tidyr::gather(trend_yr, key = "status", value = "value", -yr)

# Rename levels of "status" variable
trend_long$status <- factor(trend_long$status, levels = c("exceed", "under"),
                            labels = c(">5%", "<5%"))



## Plot for trend of studies and 5% guidance

###############################################################

# Create stacked bar plot with custom colors and renamed labels, including missing years
p <- ggplot(trend_long, aes(x = yr, y = value, fill = status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors, labels = c("> 5%", "< 5%")) + 
  labs(y = "Number of studies", title = "B") +  
  theme_minimal(base_size = 12) +
  theme(plot.margin = unit(c(1,0,0,3),"lines")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  
        panel.grid = element_blank(),  
        panel.border = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.ticks = element_line(color = "black"), 
        axis.ticks.length = unit(0.2, "cm"),  
        axis.ticks.x = element_line(),  
        axis.ticks.y = element_line(),  
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  guides(fill = FALSE)


p <- p + geom_vline(xintercept = 2014, linetype="dotted")

# calculate the proportion
glm_data$sum <- glm_data$exceed + glm_data$under
glm_data$prop5 <- glm_data$exceed/glm_data$sum
glm_data$prop5 <- as.numeric(glm_data$prop5)

(pM1 <- parameters::model_parameters(m1))

# create plot annotation
round(pM1$Coefficient[2],2)
round(pM1$SE[2],2)
round(pM1$p[2],3)

label <- "<span style='font-size:12pt; color:black'>GLM (n = 24): -0.03 ± 0.02; <i>p</i> = 0.081</span>"
#N = 24 for the 24 years (proportion > 5%)

# Plot
x <- ggplot(data = glm_data, aes(x = yr, y = prop5)) +
  geom_smooth(method = "glm", se = TRUE, color = "black", fill = "#C0C0C0", size = 0.5) +
  geom_point(size = 2) +
  labs(y = "Proportion > 5 percent (%)", title = "C") +
  theme_minimal(base_size = 12) +
  theme(plot.margin = unit(c(1,1,0,1),"lines")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = "black"),  
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

x <- x +
  annotation_custom(richtext_grob(label), xmin = 1998, xmax = 2020, ymin = 0.00, ymax = 0.55) +
  geom_vline(xintercept = 2014, linetype="dotted")


## Plot for devices and 5% guidance

###############################################################

# Selecting relevant columns and removing NAs
device <- data %>%
  select(ID, Device_type_category, H_5) %>%
  na.omit()

# Group by ID and keep unique combinations of Device_type_category
device <- device %>%
  group_by(ID) %>%
  distinct(Device_type_category, .keep_all = TRUE)

# Group by Device_type_category and H_5, then count the occurrences
device_summary <- device %>%
  group_by(Device_type_category, H_5) %>%
  summarise(count = n())

# Convert H_5 to a factor and specify the levels to reorder
device_summary$H_5 <- factor(device_summary$H_5, levels = c("Yes", "No"))

# Create a complete data frame with all combinations of Device_type_category and H_5
complete_df <- expand.grid(
  Device_type_category = unique(device_summary$Device_type_category),
  H_5 = unique(device_summary$H_5)
)

# Left join the complete data frame with the summarized data
device_summary_complete <- left_join(complete_df, device_summary, by = c("Device_type_category", "H_5"))

# Fill in NAs with 0
device_summary_complete$count[is.na(device_summary_complete$count)] <- 0

# Plot
t <- ggplot(device_summary_complete, aes(x = Device_type_category, y = count, fill = H_5)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 3) + 
  scale_fill_manual(values = my_colors, labels = c("Yes", "No")) +
  labs(y = "Frequency", title = "D") +
  theme_minimal() +
  theme(plot.margin = unit(c(1,0,-2,3),"lines"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.ticks = element_line(color = "black"), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  guides(fill = FALSE)  



## Plot for bat Families and 5% guidance

###############################################################

# Selecting relevant columns and removing NAs
family <- data %>%
  select(ID, Family, H_5) %>%
  na.omit()

# Group by ID and keep unique combinations of Family
family <- family%>%
  group_by(ID) %>%
  distinct(Family, .keep_all = TRUE)

# summarize data - used for percentages
sum(table(family$Family)) #total
table(family$Family)/sum(table(family$Family))*100


# Group by Family and H_5, then count the occurrences
family_summary <- family %>%
  group_by(Family, H_5) %>%
  summarise(count = n())

# Convert H_5 to a factor and specify the levels to reorder
family_summary$H_5 <- factor(family_summary$H_5, levels = c("Yes", "No"))


# Create a complete data frame with all combinations of Family and H_5
complete_family_df <- expand.grid(
  Family = unique(family_summary$Family),
  H_5 = unique(family_summary$H_5)
)

# Left join the complete data frame with the summarized data
family_summary_complete <- left_join(complete_family_df, family_summary, by = c("Family", "H_5"))

# Fill in NAs with 0
family_summary_complete$count[is.na(family_summary_complete$count)] <- 0

# Plot
f <- ggplot(family_summary_complete, aes(x = Family, y = count, fill = H_5)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 3) +  # Add count labels to bars
  scale_fill_manual(values = my_colors, labels = c("Yes", "No")) + # Adjust labels
  labs(x = "", y = "", title = "E") +  # Add title here
  theme_minimal() +
  theme(plot.margin = unit(c(1,1,0,1),"lines"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),  # Remove panel border
        axis.line = element_line(color = "black"),  # Add axis lines
        axis.ticks = element_line(color = "black"),  # Add axis ticks
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  guides(fill = FALSE) 

# View and save
trend <- grid.arrange(p, x, t, f, ncol = 2)
ggsave("NEWESTtrend.svg", trend, dpi = 600, width = 250, height = 200, units = "mm") 



## Data preparation for world plot and pie charts

###############################################################

# Select specific columns of data
world <- data %>% select(ID, Year_publication, Higher_geography, H_5)

# Convert 'H_5' column to character type
world$H_5 <- as.character(world$H_5)

# Replace 'NA' with 'Not reported' in column 'H_5'
world$H_5[is.na(world$H_5)] <- "Not reported"

world <- na.omit(world)

# separate data by ;
world <- tidyr::separate_rows(world, ID, Year_publication, Higher_geography, H_5, sep = " ;", convert = T)

# trim white space
trim.leading <- function (x)  sub("^\\s+", "", x)
world$Higher_geography <- trim.leading(world$Higher_geography)

# create database by distinct region
pie <- world %>%
  group_by(ID) %>%
  distinct(Higher_geography, .keep_all = TRUE)

levels(pie$Higher_geography)

pie <- pie[,-c(1,2)]

# summarize data - used for percentages for where studies occurred
sum(table(pie$Higher_geography)) #total
table(pie$Higher_geography)/sum(table(pie$Higher_geography))

summary_df <- pie %>%
  group_by(Higher_geography, H_5) %>%
  summarise(Frequency = n()) %>%
  ungroup()

# create coordinates for the pie chart locations
coordinates_df <- data.frame(
  Higher_geography = c("Afrotropical", "Australasian", "Indomalayan", "Nearctic", "Neotropical", "Palearctic"),
  long = c(21, 133.4, 79, -104, -57, 45),
  lat = c(8, -23.4, 27, 44, -16.8, 56)
)

# Merge the coordinates with the summary data
merged_df <- merge(summary_df, coordinates_df, by = "Higher_geography")

# Reshape the df
reshaped_df <- merged_df %>%
  pivot_wider(names_from = H_5, values_from = c("Frequency"), names_sep = "_") %>%
  replace_na(list(`Not reported` = 0))

# Calculate the total number of studies per region
df_sum <- reshaped_df %>%
  group_by(Higher_geography) %>%
  summarise(total_studies = sum(Yes + No + `Not reported`, na.rm = TRUE)) 

# Add the total number of studies to the original data frame
df <- left_join(reshaped_df, df_sum, by = "Higher_geography")

# Calculate proportional area for the pie chart
df$proportional_area <- sqrt(df$total_studies)

# Pull in world map data
world <- map_data("world")

# Create world map with colors
map1 <- ggplot(world, aes(long, lat)) +
  geom_map(map = world, aes(map_id = region), 
           color = "#ADADAD", fill = "#ADADAD", size = 0.3) +
  coord_quickmap() 

# Append the pie charts
map1 <- map1 + 
  geom_scatterpie(data = df, cols = c("Yes", "No", "Not reported"), color =NA, alpha=.9,
                          aes(x = long, y= lat, group = Higher_geography, r = proportional_area)) 

# plot pie charts on map
map1 <- map1 + scale_fill_manual("",labels=c("> 5%", "< 5%", "Not reported"), values=my_colors2)+
  labs(title="A") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank()) 

# Specify min and max values for the legend based on total_studies
min_total_studies <- min(df$total_studies)
max_total_studies <- max(df$total_studies)

# Calculate proportional area for the legend based on the range of total_studies
legend_proportional_area <- sqrt(c(min_total_studies, max_total_studies))

# Add the scatterpie legend with fixed min and max values and proportional area size
map1 <- map1 + 
  geom_scatterpie_legend(legend_proportional_area, x = -150, y = -45, n = 2, 
                         labeller = function(x) c(min_total_studies, max_total_studies))

# save the file
ggsave("Map.svg", map1, dpi = 600, width = 250, height = 100, units = "mm")



## JUST FOR VHF glm trend

###############################################################

# Loading the Database ----------------------------------------------------
trendx <- data %>%
  select(ID, Year_publication, Higher_geography, H_5, VHF) %>%
  filter(VHF == "Yes")

# retain distinct papers based on ID
trendx <- distinct(trendx, ID, .keep_all = TRUE) 

# remove NA values
trendx <- na.omit(trendx)

# just for model
regressx <- trendx[trendx$Year_publication > 1998 & trendx$Year_publication < 2023, ]

# Overall temporal trend in exceeding the 5% guideline
regressx <- data.frame(table(regressx$Year_publication, regressx$H_5)) ; colnames(regressx) <- c("yr","Higher than 5","N")
regressx$yr <- as.numeric(as.character(regressx$yr))

# Has the frequency of exceeding the 5% rule changed over time?
glm_datax <- data.frame(yr = unique(regressx$yr),
                       exceed = regressx[regressx$`Higher than 5`=="Yes",]$N, 
                       under = regressx[regressx$`Higher than 5`=="No",]$N)

# summary stats
sum(table(glm_datax$exceed))

# Fit the model
m2 <- glm(cbind(exceed, under) ~ yr, data = glm_datax, family = "binomial")
performance::check_overdispersion(m2)
rsq::rsq(m2)

# calculate the proportion
glm_datax$sum <- glm_datax$exceed + glm_datax$under
glm_datax$prop5 <- glm_datax$exceed/glm_datax$sum
glm_datax$prop5 <- as.numeric(glm_datax$prop5)

(pM1 <- parameters::model_parameters(m2))

# round values for manuscript
round(pM1$Coefficient[2],2)
round(pM1$SE[2],2) 
round(pM1$p[2],3) 
