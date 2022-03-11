# set working drive, personal mac
setwd("___") #Edit for your working drive

# load packages
library(tidyverse) # data tidying, subsetting, transformations, etc.
library(emmeans) # estimated marginal means, good for unbalanced data
library(ggpubr)
library(data.table)
library(rstatix) # pretty p-values
library(lme4) # mixed linear models
require(lattice)
library(predictmeans)
library(outliers) # Dixon test
library(plotrix) # standard error


# Color palette for plotting
pal_okabe_ito <- c(
  "#E69F00", #orange
  "#56B4E9", #blue
  "#009E73", #green
  "#F0E442", #yellow
  "#0072B2", #drkblue
  "#D55E00", #drkorange
  "#CC79A7" #pink
)
pal_okabe_ito_blue <- pal_okabe_ito[c(5,6,1,2,3,7,4)] 
pal_okabe_ito_red <- pal_okabe_ito[c(6,5,3,1,2,7,4)] 
pal_okabe_ito_2 <- pal_okabe_ito[c(5,6)]
pal_okabe_ito_3 <- pal_okabe_ito[c(5,6,7)]
pal_okabe_ito_3_light <- pal_okabe_ito[c(1,2,7)]
pal_okabe_ito_4 <- pal_okabe_ito[c(5,6,7,2)]

# reads your CSV files
gels<-read.csv("_____.csv", header=T, na.strings=c("","NA")) #Add raw data file
head(gels)
gels<-select(gels, sample_id, crop, sample.weight, ext.rep, gel.rep, ip6.ug, notes)
head(gels)

lodloq<-read.csv("____.csv", header=T)  # Raw data file for samples that were assayed to determine LOD and LOQ
head(lodloq)
lodloq<-select(lodloq, Lane, abs.quant)
head(lodloq)



######################################################
############### LOD LOQ calculations ##################
######################################################

# LOD and LOQ are calculated using one sample replicated on one gel
# The gel results used for the following LOD and LOQ calculations was Gel #1 analyzed on 1/23/2022
# File name: Gel 1_2022-01-23 10hr 37min 14sec.xlxs
# File location: /Users/juliannekellogg/Documents/2. Ecuador Research/1_Phytic Acid analyses/Gel images/Ecuador samples/22_01_23_JAKstained/Gel 1_Chemidoc Touch Images 2022-01-23_10.37.55

# Calibration curve
# R^2 calculated using Image Lab software (Bio-Rad)
# R^2 = 0.919711

# Calculate standard deviation for use in LOD and LOQ calculations 
stats.lodloq<-lodloq %>%  
  filter(Lane == "7" | Lane == "8" | Lane == "9"|Lane == "10"|Lane == "11"|Lane == "12") %>%  
  mutate(ip6.mg.g = abs.quant*(1/0.01)*(1/(0.05*1000))) %>%   #to get mg/g which is often reported
  summarise( # summarize by group
    mean.ug = mean(abs.quant, na.rm = TRUE),
    sd.ug = sd(abs.quant, na.rm = TRUE),
    mean.mg.g = mean(ip6.mg.g, na.rm = TRUE),
    sd.mg.g = sd(ip6.mg.g, na.rm = TRUE)) %>%
  mutate(rsd.ug = (sd.ug/mean.ug)*100) %>%
  mutate(rsd.mg.g = (sd.mg.g/mean.mg.g)*100) #RSD %
stats.lodloq
# SD ug = 0.164218
#SD mg/g = 0.328436

# Limit of Detection (linear calibration curve)
# the lowest concentration of an analyte in a sample that can be detected, 
# but not necessarily quantified, under the stated conditions of the test
# LOD=3Sa/b
# Sa is the standard deviation of the response
# b is the slope of the calibration curve
3*0.164218/0.919711
# LOD = 0.5356617 in ug
3*0.328436/0.919711
# LOD = 1.071323 in mg/g

# Limit of Quantification (linear calibration curve)
# the lowest concentration of an analyte in a sample that can be determined 
# with acceptable precision and accuracy under the stated conditions of test
# LOQ=10Sa/b 
10*0.164218/0.919711
# LOQ = 1.785539 in ug
10*0.328436/0.919711
# LOQ = 3.571078 in mg/g



########################################################
################### Data tidying ###################
########################################################

# Goal: 
#remove extraction 3 replicate (this extraction is for future analyses)
gels.ext1and2 <- gels %>%
  filter(ext.rep == "1"|ext.rep == "2") 

### IMPUTATION ####
# Imputation help: Palarea-Albaladejo and Martin-Fernandez, 2013
# impute 65% LOD for the NAs that are non-detect ("ND" in this data set)
# impute 65% LOD for values < LOD
# LOD = 0.5356617 ug
0.5356617*0.65
# 65% LOD = 0.3481801

gels.ext1and2.imputed <- gels.ext1and2 %>%
  mutate(
    imputed.ip6.ug = case_when(
      ip6.ug < 0.5356617 ~ 0.3481801,
      ip6.ug == "ND" ~ 0.3481801,
      TRUE ~ as.numeric(ip6.ug)))

gels.ext1and2.imputed %>%
  count(imputed.ip6.ug == 0.3481801)

gels.ext1and2.imputed %>%
  count(imputed.ip6.ug == "NA")

25/255*100 # percentage of data that was <LOD
#9.8%


########################################################
################### Intermediate precision ##################
########################################################

# The precision determined from replicate determinations 
# conducted within a single laboratory not simultaneously
# To be estimated using samples with at least 5 replicates

# Sample 102 has 5 replicates
# This sample will cover all possible changes in methodology 
# (e.g., changes in tips, acrylamide, stain shaking, person imaging, etc.)
# Use RSD % for precision 

interm.precision <- gels.ext1and2.imputed %>%
  filter(sample_id == "102") %>%
  summarise( 
    mean = mean(imputed.ip6.ug, na.rm = TRUE),
    sd = sd(imputed.ip6.ug, na.rm = TRUE)) %>%
  mutate(rsd = (sd/mean)*100) #RSD %
# RSD% = 52.5%
# This is high; RSD % should be between 2-3% for samples with 
# analyte concentrations of 0.1-1% (AOAC official methods of analysis, 
# 2019, appendix K, guidelines for dietary supplements and botanicals)



######################################################
################## IP6 CALCULATIONS #################
######################################################

### Calculate amt of IP6 in samples ###
gels<-gels.ext1and2.imputed %>%
  select(sample_id, crop, sample.weight, ext.rep, gel.rep, imputed.ip6.ug) %>% #do not include notes column
  mutate(ip6.mg.g = imputed.ip6.ug*(1/0.01)*(1/((gels.ext1and2$sample.weight)*1000))) %>% #to get mg/g which is often reported
  mutate(ip6.mg.kg = ip6.mg.g*1000) # to get mg/kg (ppm) to help assess RSD%

### Add variety back in as a factor ###
# grab the .csv file generated from the Ecuador Sample Collection Analysis RStudio project
var<-read.csv("crop samples and variety names.csv", header=T)
var<- var %>%
  select(id.1, var.clean)%>%
  rename(sample_id = id.1)
# combine data sets
gels<- left_join(gels, var, by="sample_id")


#Random effects = ext.rep, gel.rep, residual
gels$ext.rep<-as.factor(gels$ext.rep)
gels$gel.rep<-as.factor(gels$gel.rep)
#Fixed effects = sample_id, crop
gels$sample_id<-as.factor(gels$sample_id)
gels$crop<-as.factor(gels$crop)
gels$var.clean<-as.factor(gels$var.clean)


### GET SOME SUMMARY STATS ###
stats.gels<-gels %>%  # ok to use to get SD and RSD, don't use this for means and comparisons if data is unbalanced
  group_by(sample_id) %>% # grouping #can choose sample_id or ext.rep
  summarise( # summarize by group
    mean = mean(ip6.mg.g, na.rm = TRUE),
    sd = sd(ip6.mg.g, na.rm = TRUE)) %>%
  mutate(rsd = (sd/mean)*100) #RSD %
stats.gels



#### OUTLIERS ####
# Dixon's test for outliers (must have a minimum of 3 observations)
# Identify samples that you are concerned have outlier values; select samples with high RSD
# Look at stats.gels dataframe 
subset.gels<-gels %>%
  filter(sample_id == 904) #working from highest to lowest RSD, stopped at 
head(subset.gels)
dixon.high <- dixon.test(subset.gels$ip6.mg.g)
dixon.high # if p = 0.05 or is < 0.05, then the value indicated is an outlier
dixon.low <- dixon.test(subset.gels$ip6.mg.g, opposite = TRUE)
dixon.low
# look at extreme values that aren't in agreement with AT LEAST ONE OTHER VALUE
# Outlier detection has been completed! (1/31/2022)
# 218 highest value 20.716142 is an outlier (Row 95)
# 308 highest value 22.728088 is an outlier (Row 140)
# 309 highest value 13.3482528301887 is an outlier (Row 144)
# 404 highest value 8.621 is an outlier (Row 199)
# 601 highest value 13.400558490566 is an outlier (235)

# Remove outliers from data set
# work with "gels" data set
# Remove rows identified above
gels_new <- gels[-c(95, 140, 144, 199, 235), ]
  
#### RSD ####
# Redo summary stats for RSD values
stats.gels.new<-gels_new %>%  # ok to use to get SD and RSD, don't use this for means and comparisons if data is unbalanced
  group_by(sample_id) %>% # grouping #can choose sample_id or ext.rep or 
  summarise( # summarize by group
    mean = mean(ip6.mg.g, na.rm = TRUE),
    sd = sd(ip6.mg.g, na.rm = TRUE),
    se = std.error(ip6.mg.g, na.rm = TRUE)) %>%
  mutate(rsd = (sd/mean)*100) #RSD %
stats.gels.new
### Review RSD data
hist(stats.gels.new$rsd)
## if you want RSD by crop, then average RSD within crop
crops <- gels %>%
  select (sample_id, crop)
stats.gels.new<- left_join(stats.gels.new, crops, by="sample_id") # add crop names to stats data
var$sample_id <- as.factor(var$sample_id)
stats.gels.new<- left_join(stats.gels.new, var, by="sample_id") # add variety names to stats data
#by crop
stats.gels.avg.crop<-stats.gels.new %>%
  group_by(crop) %>%
  summarise(
    mean.rsd = mean(rsd, na.rm = TRUE),
    median = median(mean,na.rm = TRUE),
    mean = mean(mean,na.rm = TRUE)
    )
stats.gels.avg.crop
# by variety
stats.gels.avg.var<-stats.gels.new %>%
  group_by(var.clean) %>%
  summarise(
    mean.rsd = mean(rsd, na.rm = TRUE),
    median = median(mean,na.rm = TRUE),
    mean = mean(mean,na.rm = TRUE)
  )
stats.gels.avg.var
# figure out SE for var
stats.gels.se.var<-gels_new %>%  # ok to use to get SD and RSD, don't use this for means and comparisons if data is unbalanced
  group_by(var.clean) %>% # grouping
  summarise( # summarize by group
    median = median(ip6.mg.g,na.rm = TRUE),
    mean = mean(ip6.mg.g, na.rm = TRUE),
    sd = sd(ip6.mg.g, na.rm = TRUE),
    se = std.error(ip6.mg.g, na.rm = TRUE)) %>%
  mutate(rsd = (sd/mean)*100) #RSD %
stats.gels.se.var



#### SUBSET CROP SPECIES FOR SEPARATE DATA FRAMES ####
amaranth <- gels_new %>%
  filter(crop == "amaranth") # Do not analyze. Only one sample for the var Alegria
stats.am<-amaranth %>%  
  group_by(var.clean) %>%  
  summarise( 
    median = median(ip6.mg.g,na.rm = TRUE),
    mean = mean(ip6.mg.g, na.rm = TRUE),
    sd = sd(ip6.mg.g, na.rm = TRUE),
    se = std.error(ip6.mg.g, na.rm = TRUE)) %>%
  mutate(rsd = (sd/mean)*100) 
stats.am

barley <- gels_new %>%
  filter(crop == "barley") 
stats.bar<-barley %>%  
  group_by(var.clean) %>%  
  summarise( 
    median = median(ip6.mg.g,na.rm = TRUE),
    mean = mean(ip6.mg.g, na.rm = TRUE),
    sd = sd(ip6.mg.g, na.rm = TRUE),
    se = std.error(ip6.mg.g, na.rm = TRUE)) %>%
  mutate(rsd = (sd/mean)*100) 
stats.bar
barley <- gels_new %>%
  filter(crop == "barley") %>% # do not include perla or VNS
  filter(var.clean == "Cañicapac" | var.clean == "Pelada" | var.clean == "Pelada local" | var.clean == "Shiry")

fava <- gels_new %>%
  filter(crop == "fava") # do not analyze
stats.fav<-fava %>%  
  group_by(var.clean) %>%  
  summarise( 
    median = median(ip6.mg.g,na.rm = TRUE),
    mean = mean(ip6.mg.g, na.rm = TRUE),
    sd = sd(ip6.mg.g, na.rm = TRUE),
    se = std.error(ip6.mg.g, na.rm = TRUE)) %>%
  mutate(rsd = (sd/mean)*100) 
stats.fav

lentil <- gels_new %>%
  filter(crop == "lentil") # do not analyze; remove from data set
stats.len<-lentil %>%  
  group_by(var.clean) %>%  
  summarise( 
    median = median(ip6.mg.g,na.rm = TRUE),
    mean = mean(ip6.mg.g, na.rm = TRUE),
    sd = sd(ip6.mg.g, na.rm = TRUE),
    se = std.error(ip6.mg.g, na.rm = TRUE)) %>%
  mutate(rsd = (sd/mean)*100) 
stats.len

lupin <- gels_new %>%
  filter(crop == "lupin") # good
stats.lup<-lupin %>%  
  group_by(var.clean) %>%  
  summarise( 
    median = median(ip6.mg.g,na.rm = TRUE),
    mean = mean(ip6.mg.g, na.rm = TRUE),
    sd = sd(ip6.mg.g, na.rm = TRUE),
    se = std.error(ip6.mg.g, na.rm = TRUE)) %>%
  mutate(rsd = (sd/mean)*100) 
stats.lup

maize <- gels_new %>%
  filter(crop == "maize")
stats.mai<-maize %>%  
  group_by(var.clean) %>%  
  summarise( 
    median = median(ip6.mg.g,na.rm = TRUE),
    mean = mean(ip6.mg.g, na.rm = TRUE),
    sd = sd(ip6.mg.g, na.rm = TRUE),
    se = std.error(ip6.mg.g, na.rm = TRUE)) %>%
  mutate(rsd = (sd/mean)*100) 
stats.mai
maize <- gels_new %>%
  filter(crop == "maize") %>% # do not include cuzco, runapa shunku, wagal
  filter(var.clean == "Morocho" | var.clean == "Zhima")

pea <- gels_new %>%
  filter(crop == "pea") # do not analyze
stats.pea<-pea %>%  
  group_by(var.clean) %>%  
  summarise( 
    median = median(ip6.mg.g,na.rm = TRUE),
    mean = mean(ip6.mg.g, na.rm = TRUE),
    sd = sd(ip6.mg.g, na.rm = TRUE),
    se = std.error(ip6.mg.g, na.rm = TRUE)) %>%
  mutate(rsd = (sd/mean)*100) 
stats.pea

quinoa <- gels_new %>%
  filter(crop == "quinoa") # do not analyze
stats.qui<-quinoa %>%  
  group_by(var.clean) %>%  
  summarise( 
    median = median(ip6.mg.g,na.rm = TRUE),
    mean = mean(ip6.mg.g, na.rm = TRUE),
    sd = sd(ip6.mg.g, na.rm = TRUE),
    se = std.error(ip6.mg.g, na.rm = TRUE)) %>%
  mutate(rsd = (sd/mean)*100) 
stats.qui

wheat <- gels_new %>%
  filter(crop == "wheat")
stats.wheat<-wheat %>%  
  group_by(var.clean) %>%  
  summarise( 
    median = median(ip6.mg.g,na.rm = TRUE),
    mean = mean(ip6.mg.g, na.rm = TRUE),
    sd = sd(ip6.mg.g, na.rm = TRUE),
    se = std.error(ip6.mg.g, na.rm = TRUE)) %>%
  mutate(rsd = (sd/mean)*100) 
stats.wheat
  wheat <- gels_new %>%
  filter(crop == "wheat") %>% #do not include vivar or VNS
  filter(var.clean == "Napo" | var.clean == "Oregon" )

#### FIGURE FOR VARIETIES ####
# Remove all VNS - Allows me to then visually compare across named varieties 
  gels_new_noVNS <- filter(gels_new, 
                                    !var.clean %in% 
                                      c('VNS')) # exclude 'VNS' samples
  gels_new_noVNS$crop <- as.factor(gels_new_noVNS$crop)
  ip6_varieties <- ggplot(data = gels_new_noVNS,
                         aes(x = var.clean,
                             y = ip6.mg.g)) + #add data points
    geom_point() +
  ylab(expression(paste("IP6 (mg/g)"))) +  # change the title of the y-axis
  theme_pubr() +  # minimalist theme
  theme(
    legend.position = "none", # remove legend
    axis.title.x = element_blank() # remove x-axis title
    ) +
    #scale_color_manual(values = pal_okabe_ito_blue) +   # use custom color palette for the points
    #expand_limits(y=c(50,2000)) + #expand_limits(y=c(50,5000)) + # set y-axis limits
    NULL
ip6_varieties 
  
  

### REMOVE LENTIL ###
#Lentil stats (use these for results)
# Will not include lentil in analysis because there was only one sample
# Lentil mean = 5.92 ip6 mg/g
# Lentil sd = 4.85
# Lentil rsd% = 86.18
gels_new <- gels[-c(229:232), ]


##########################################################
############## STUDY HYPOTHESES  #######################
##########################################################
##########################################################

## Hypothesis 1: What are the differences among crops?
# Pairwise comparisons among crops
lmer.crop<-lmer(log(ip6.mg.g) ~  crop + (1|sample_id/ext.rep), #gel rep nested within extraction rep, but gel.rep isn't specified because it's the lowest level of replication
            data=gels_new, #Use the updated data set with the outliers removed
            na.action = na.omit)
summary(lmer.crop)


## Hypothesis 2: What are the differences among varieties, within crops?
# t-test comparisons between varieties with more than one sample (within crops)
lmer.bar<-lmer((ip6.mg.g) ~  var.clean + (1|sample_id/ext.rep), 
                data=barley, 
                na.action = na.omit)
t.test(ip6.mg.g ~ var.clean, data = wheat)
t.test(ip6.mg.g ~ var.clean, data = maize)
t.test(ip6.mg.g ~ var.clean, data = lupin)

##########################################################
########### ASSUMPTIONS  #######################
##########################################################
##########################################################

### CHECKING STATISTICAL ASSUMPTIONS USING LMER MODELS ###
plot(log(ip6.mg.g) ~ crop, data = gels_new)
plot((ip6.mg.g) ~ var.clean, data = barley)
#Good resource: https://ademos.people.uic.edu/Chapter18.html
plot(resid(lmer.crop)) # check linearity
plot(resid(lmer.bar))
### NORMALITY ####
qqnorm(resid(lmer.crop)) # check normality of residuals
qqmath(lmer.crop, id=0.05) # look for normality and extreme values
qqnorm(resid(lmer.bar)) # 
qqmath(lmer.bar, id=0.05)
### HOMOGENEITY OF VARIANCE ####
plot(lmer.crop) # check heteroscedacity
plot(lmer.bar)
#Bartlett’s test is sensitive to deviations from normality. 
#If you’ve verified that your data is normally distributed then Bartlett’s test 
#is a good first test but if your data is non-normal than you should always 
#verify results from Bartlett’s test to results from Levene’s and Flinger-Killeen’s 
#tests as there is a good chance you will get a false positive.
#source: https://uc-r.github.io/assumptions_homogeneity#levene
bartlett.test(log(ip6.mg.g) ~ crop, data=gels_new)
bartlett.test((ip6.mg.g) ~ var.clean, data=barley)
#The Levene’s test is slightly more robust to departures from normality than the Bartlett’s test. 
#Levene’s performs a one-way ANOVA conducted on the deviation scores; that is, the absolute difference 
#between each score and the mean of the group from which it came.1 To test, we use leveneTest() from 
#the car package. By default leveneTest() will test variance around the median but you can override this 
#by using the center = mean argument. 
#source: https://uc-r.github.io/assumptions_homogeneity#levene
library(car)
leveneTest(log(ip6.mg.g) ~ crop, data=gels_new)
leveneTest((ip6.mg.g) ~ var.clean, data=barley)
#Fligner-Killeen's test of homogeneity of variance, one per independent variable in each favored model ..:
#Less sensitive (and thus more reliable) for outliers than Levene's test
#Can also handle continuous independent variables
#Anything above p=0.05 is ok
#source: https://rpubs.com/loveb/mm
fligner.test(log(ip6.mg.g) ~ crop, data=gels_new)
fligner.test((ip6.mg.g) ~ var.clean,data=barley)
# Rejected the assumption of equal variance before log transformation

###### TRANSFORMATIONS ########
# Crop data needed to be transformed
# (log()) transformation applied to crop data,  rechecked assumptions with log(),  all assumptions met
# Barley data meets assumptions


##########################################################
############### CROPS ANALYSIS ##########################
##########################################################
##########################################################


### MIXED EFFECTS LINEAR MODEL ###
### ANOVA ###
### CROPS DATA LOG TRANSFORMED ###
### BARLEY DATA RAW ###
#fit the one-level mixed model with REML using the lmer function from the lme4 library
# Info on crossed vs nested random effects: https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified
library(lmerTest)
lmer.crop<-lmer(log(ip6.mg.g) ~  crop + (1|sample_id/ext.rep), #gel rep nested within extraction rep, but gel.rep isn't specified because it's the lowest level of replication
                data=gels_new, #Use the updated data set with the outliers removed
                na.action = na.omit)
summary(lmer.crop)
anova(lmer.crop) # Satterthwaite's method
anovalmer(lmer.crop) # Kenward-Roger approximation for degrees of freedom 
# Satterthwaite's method: Rejected the null hypothesis of no differences among crops; p = 0.0001413

lmer.bar<-lmer((ip6.mg.g) ~  var.clean + (1|sample_id/ext.rep), 
               data=barley, 
               na.action = na.omit)
anova(lmer.bar) # Satterthwaite's method


### ESTIMATED MARGINAL MEANS ###
### Pairwise comparisons among crops ###
# adapted from code source: https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/plotting-models.html
#Learn about EMMEANS and data transformations here: https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html
## Using this because the response is a variable that is the log of some other variable
emm.crop<-emmeans(lmer.crop, ~ crop, type = "response") #estimated marginal means because data is unbalanced
# pairwise with different commands
ip6_simple <- contrast(emm.crop,
                       method = "pairwise",
                       simple = "each",
                       combine = TRUE,
                       adjust = "none") %>%
  summary(infer = TRUE)
ip6_simple


##########################################################
########### CROPS ANALYSIS FIGURES #######################
##########################################################
##########################################################


### BASE PLOT USING GGPLOT2
ip6_response <- ggplot(data = gels_new,
                       aes(x = crop,
                           y = ip6.mg.g)) + #add data points
  #geom_jitter(width=0.2, shape=1) + #random position of data points
  geom_violin() +
  ylab(expression(paste("IP6 (mg/g)"))) +  # change the title of the y-axis
  theme_pubr() +  # minimalist theme
  theme(
    legend.position = "none", # remove legend
    axis.title.x = element_blank() # remove x-axis title
  ) +
  #scale_color_manual(values = pal_okabe_ito_blue) +   # use custom color palette for the points
  #expand_limits(y=c(50,2000)) + #expand_limits(y=c(50,5000)) + # set y-axis limits
  NULL
ip6_response

# convert bbg_emm to a data.table
emm_dt_crop <- summary(emm.crop) %>%
  data.table()

#Add modeled means and CIs
ip6_response <- ip6_response + 
  geom_errorbar(data = emm_dt_crop, 
                aes(y = response,# This differs from the non-transformed emm data # previously: y = emmean
                    ymin = lower.CL, 
                    ymax = upper.CL), 
                width = 0.1) +  # add layer containing error bars
  geom_point(data = emm_dt_crop,
             aes(y = response),# This differs from the non-transformed emm data # previously: y = emmean
             color = "#D55E00",
             size = 1) + # add layer containing means
  NULL
ip6_response

#Adding p-value brackets to response plot
#This function needs a column of p-values, and a pair of columns that define the left and right 
ip6_simple_dt <- data.table(ip6_simple)
# Group1: column containing x-position of the left side of the p-value bracket
# Use ip6_simple_dt
ip6_simple_dt[, group1 := c("amaranth",
                            "amaranth",
                            "amaranth",
                            "amaranth",
                            "amaranth",
                            "amaranth",
                            "amaranth",
                            "barley",
                            "barley",
                            "barley",
                            "barley",
                            "barley",
                            "barley",
                            "fava",
                            "fava",
                            "fava",
                            "fava",
                            "fava",
                            "lupin",
                            "lupin",
                            "lupin",
                            "lupin",
                            "maize",
                            "maize",
                            "maize",
                            "pea",
                            "pea",
                            "quinoa"
                            )]

# Group2: column containing x-position of the right side of the bracket
ip6_simple_dt[, group2 := c("barley",
                            "fava",
                            "lupin",
                            "maize",
                            "pea",
                            "quinoa",
                            "wheat",
                            "fava",
                            "lupin",
                            "maize",
                            "pea",
                            "quinoa",
                            "wheat",
                            "lupin",
                            "maize",
                            "pea",
                            "quinoa",
                            "wheat",
                            "maize",
                            "pea",
                            "quinoa",
                            "wheat",
                            "pea",
                            "quinoa",
                            "wheat",
                            "quinoa",
                            "wheat",
                            "wheat")]

ip6_simple_dt[, p_rounded := p_round(p.value,
                                     digits = 3)]
ip6_simple_dt[, p_pretty := p_format(p_rounded,
                                     digits = 3,
                                     accuracy = 1e-03,
                                     add.p = TRUE)]

# filter out p-values > 0.05
ip6_simple_dt_sig <- ip6_simple_dt %>%
  filter(p_rounded < 0.055)

# add the p-values
ip6_response_p <- ip6_response +
  stat_pvalue_manual(
    data = ip6_simple_dt_sig,
    label = "p_pretty",
    y.position = c(31, 33, 35, 37,39, 41, 43, 45, #position of p-values  
              47, 49, 51, 53, 55),
 size=2.4,
    tip.length = 0.01)
ip6_response_p

# reset working drive for figures
setwd("/Users/juliannekellogg/Documents/2. Ecuador Research/6_Data analysis_all outcomes/Ecuador Phytic Acid Analysis/Output")

# improve look
ip6_response_p_pretty <- ip6_response_p + #change axis labels and fonts/sizes
  scale_x_discrete(labels=c("amaranth" = "Amaranth", "barley" = "Barley","fava" = "Fava", 
                            "lupin" = "Lupin","maize" = "Maize", "pea" = "Pea",
                            "quinoa" = "Quinoa", "wheat" = "Wheat")) +
  theme(strip.text.x = element_text(size = 9, family ="Arial", color = "black" ),
        strip.text.y = element_text(size = 9, family ="Arial", color = "black"),
        axis.title.y = element_text(family = "Arial", size = 9, color = "black"),
        axis.text = element_text(family = "Arial", size = 8))+
  annotate("text", x=1, y=4, label="n = 3", size = 6/.pt) +
  annotate("text", x=2, y=-0.5, label="n = 22", size = 6/.pt) +
  annotate("text", x=3, y=-0.5, label="n = 11", size = 6/.pt) +
  annotate("text", x=4, y=-0.5, label="n = 5", size = 6/.pt) +
  annotate("text", x=5, y=-0.5, label="n = 17", size = 6/.pt) +
  annotate("text", x=6, y=-0.5, label="n = 5", size = 6/.pt) +
  annotate("text", x=7, y=3, label="n = 6", size = 6/.pt) +
  annotate("text", x=8, y=-0.5, label="n = 7", size = 6/.pt) 

ip6_response_p_pretty

ggsave(file = "crop_ip6_differences.jpeg", width = 110, height = 90, units = "mm")




##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################










