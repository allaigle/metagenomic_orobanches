############################################
##### Complete pipeline of the project #####
############################################

### Exception : UMAP using python (GoogleColab)


################################   README   ####################################
# In our case, we copied each raw data to conserve them
# Created a new folder, called data_modified 

# Note 1: All spaces/"_" have been replaced by dots
#   "Sample.ID"	"Week"	"Replicate"	"Soil.Infestation"	"Total.Fixation"
#   "Root.DM"	"Shoot.DM"	"Rel.abs.Elev" columns have been kept 
#         -> new csv file : metadata_16S.csv

# Note 2: in the new "metadata_16S.csv" : 
#         for "Week" column : int changed into strings (e.g. : 0 -> S0) 

# Note 3: not all "save the file" codes are written here, only the more useful ones

# Note 4: all phyloseq objects are clustered into the phyloseq_objects_created.R 
#           script          

# Note 5: libraries are called only once here, but they are called in each 
#           script if they are necessary

# Note 6: from abundance_table_16s_Mursay.csv : otu_table_16S.csv 
#         and taxa_16S_Mursay.csv as taxa_16S.csv
################################################################################


##### I. Metadata analysis #####
#   A. Trim metadata

# trim_metadata.R script

# Load library
library(readr) #to read excel files

## Upload metadata table 
# Warning : delim are semicolons
Tmetadata <- read_delim("data_modified/metadata_16S.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(Tmetadata)

################
### TRIM DATA ##
################
## Trim depending on soil
# keep rows that contains the conditions "Conducive"
metadata_CN <- subset(Tmetadata, Soil.Infestation == "Conducive") 
View(metadata_CN)

# keep rows that contains the condition "Suppressive"
metadata_SN <- subset(Tmetadata, Soil.Infestation == "Suppressive")

# keep rows depending on weeks 
metadata_W2 <- subset(Tmetadata, Week == "S2") # Week 2

metadata_W5 <- subset(Tmetadata, Week == "S5") # Week 5

metadata_W7 <- subset(Tmetadata, Week == "S7") # Week 7

## Depending on soils AND weeks
# based on metadata_CN
metadata_CN_W2 <- subset(metadata_CN, Week == "S2") # CN and Week 2
View(metadata_CN_W2)

metadata_CN_W5 <- subset(metadata_CN, Week == "S5") # CN and Week 5

metadata_CN_W7 <- subset(metadata_CN, Week == "S7") # CN and Week 7

# based on metadata_SN
metadata_SN_W2 <- subset(metadata_SN, Week == "S2") # SN and Week 2

metadata_SN_W5 <- subset(metadata_SN, Week == "S5") # SN and Week 5

metadata_SN_W7 <- subset(metadata_SN, Week == "S7") # SN and Week 7

# Also, need to remove rows where the week is S0
# For future boxplots
# if necessary : "library(dplyr)"
trim_meta = Tmetadata[- grep("S0", Tmetadata$Week),] 
View(trim_meta) 
# trim_meta saved as "./results/metadata_analysis/metadata_wo_S0.csv"


#   B. t-test on total fixation
# t_test_meta_TotalFixation.R script 

#############
## T-tests ##
#############
# depending on soils only
t.test(metadata_SN$Total.Fixation,metadata_CN$Total.Fixation)
# p-value = 0.001442

# depending on weeks only
t.test(metadata_W2$Total.Fixation,metadata_W5$Total.Fixation) #W2 versus W5
# p-value = 6.506e-06

t.test(metadata_W2$Total.Fixation,metadata_W7$Total.Fixation) #W2 versus W7
# p-value = p-value = 0.1033

t.test(metadata_W5$Total.Fixation,metadata_W7$Total.Fixation) #W5 versus W7
# p-value = p-value = 0.000375

# depending on soils AND weeks
t.test(metadata_CN_W2$Total.Fixation,metadata_SN_W2$Total.Fixation) #CN_W2 versus SN_W2
# p-value = 0.001718

t.test(metadata_CN_W5$Total.Fixation,metadata_SN_W5$Total.Fixation) #CN_W5 versus SN_W5
# p-value = 0.004538

t.test(metadata_CN_W7$Total.Fixation,metadata_SN_W7$Total.Fixation) #CN_W7 versus SN_W7
# p-value = 0.1521


#   C. Plot results on metadata
# plots_metadata.R script

## Load libraries
#library(readr) 
library(ggplot2)
theme_set(theme_bw()) # Define theme for plots

# Upload data : "./results/metadata_analysis/metadata_wo_S0.csv" as trim_meta

###########
## STATS ##
###########
## Function to produce summary statistics (mean and +/- sd) #### 
# to have stats incorporated to next plots 
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

## Plot using qplot function (test)
# trim_meta from I.A (without S0 neither Smoins1)
qplot(Week, Total.Fixation, data = trim_meta) + facet_wrap(~Soil.Infestation) 

##############
## BOXPLOTS ##
##############
## Boxplots of the totalFixation depending only on soils

# Change box plot colors by groups
boxplot_soil <- ggplot(trim_meta, aes(x = "Soil.Infestation", y = Total.Fixation, fill = Soil.Infestation)) + geom_boxplot()
boxplot_soil

boxplot_week <- ggplot(trim_meta, aes(x = "Week", y = Total.Fixation, fill = Week)) + geom_boxplot()
boxplot_week

boxplot_soil_week <- ggplot(trim_meta, aes(x = "Week", y = Total.Fixation, fill = Week)) + geom_boxplot() + facet_wrap(~Soil.Infestation)
boxplot_soil_week



##### II. Alpha-diversity #####

#   A. For soils    
# shannon_index_soils.R script

####################
## LOAD LIBRARIES ##
####################
# If necessary : install phyloseq 
# source(' http://bioconductor.org/biocLite.R')
# biocLite('phyloseq')
library(phyloseq)
# packageVersion('phyloseq') # Check version
library(dplyr) # to use the %>% symbol and "bind_rows" function

## Load custom phyloseq functions (feat. Dr. Mahendra Mariadasou INRA)
#WARNING : check paths inside the Custom_Functions.R file for the other scripts
source("./scripts/phyloseq_customs/Custom_Functions.R")
source("./scripts/phyloseq_customs/phyloseq_object_to_df.R")#export phylo object

###############
## LOAD DATA ##
###############
# better to do it by your own
otu_table <- read_delim("data_modified/otu_table_16S.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(otu_table) # Check how they look

taxa_table<- read_delim("data_modified/taxa_16S.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(taxa_table)

###############
## TRIM DATA ##
###############
## Trim OTU table --> Some data that we don't want
# Search columns where the name sample contains "Smoins1" (week-1)
test_grep_Smoins1 <- grep("Smoins1", names(otu_table), value = TRUE)
test_grep_Smoins1

# Remove columns containing "Smoins1"
otu_tab <- otu_table %>% select(-contains("Smoins1")) # do not forget the "-"

# Remove first column (numbers, not OTU.ID) that was automatically created
otu_tab <- otu_table %>% select(-contains("...1"))

# Save otu_tab 
write.csv(otu_tab,'./data_modified/otu_tab.csv')

# Verify : 
verif_search_Smoins1 <- grep("Smoins1", names(otu_tab), value = TRUE)
verif_search_Smoins1    # Worked : character(0)

## Trim depending on soil
# create a table for only OTUs and columns with the "conducive" soil
# by removing all columns with the conducive soil ('CN')
otu_CN <- otu_tab %>% select(contains("CN"))

# create a table for only OTUs and columns with the "suppressive" soil
otu_SN <- otu_tab %>% select(contains("SN"))

########################
## CREATE A PHYLOSEQ  ##   
##      OBJECT        ##
########################
# without metadata 
soil_CN <- phyloseq(otu_table(as.matrix(otu_CN), taxa_are_rows = TRUE),
                    tax_table(as.matrix(taxa_table))) # data.frame, not matrix
soil_CN
head(otu_table(soil_CN))

soil_SN <- phyloseq(otu_table(as.matrix(otu_SN), taxa_are_rows = TRUE),
                    tax_table(as.matrix(taxa_table))) # data.frame, not matrix
soil_SN

#####################
## ALPHA DIVERSITY ##
##    RICHNESS     ##
#####################
## Numeric values of alpha diversity indices for CN and SN
shannon.CN <- estimate_richness(soil_CN, measures = "Shannon")
head(shannon.CN)

shannon.SN <- estimate_richness(soil_SN,measures = "Shannon")
head(shannon.SN)

#################
## SAVE TABLES ##
#################
write.csv(shannon.CN,'./results/alpha_diversity/shannon_CN.csv')
write.csv(shannon.SN,'./results/alpha_diversity/shannon_SN.csv')

#####################
## CALCULATE MEAN  ##
##   of SHANNON    ##
#####################
Mean.Shannon.CN <- mean(shannon.CN$Shannon)
print(Mean.Shannon.CN)
# result : 5.639646

Mean.Shannon.SN <-  mean(shannon.SN$Shannon)
print(Mean.Shannon.SN) # result : 5.693977

############
## T-TEST ##
############
#To compare both Shannon indexes of CN and SN soils 
t.test(shannon.CN$Shannon,shannon.SN$Shannon)
# not significative : p-value = 0.3859

# conclusion : both CN and SN soils have a similar structure (alpha diversity)

## Boxplot 
# Add a column in each data frame of Soil Infestation
shannon.CN$Soil.Infestation <- "Conducive"
shannon.SN$Soil.Infestation <- "Suppressive"
View(shannon.CN) # two  columns except headers : Shannon and Soil.Infestation
View(shannon.SN)

alpha_div <- bind_rows(shannon.CN, shannon.SN) 
View(alpha_div) # worked
class(alpha_div) # data.frame 

## Boxplot
boxplot_shannon <- ggplot(alpha_div, aes(x = "Soil.Infestation", y = Shannon, fill = Soil.Infestation)) + geom_boxplot()
boxplot_shannon 


#   B. For weeks
# shannon_index_weeks.R script
# as for soils, used otu_table and taxa_table

###############
## TRIM DATA ##
###############

## Trim depending on weeks
# create a table for only OTUs and columns with : 
otu_Smoins1 <- otu_table %>% select(contains("Smoins1")) # Smoins1
View(otu_Smoins1)

otu_S0 <- otu_table %>% select(contains("S0")) # S0

otu_S2 <- otu_table %>% select(contains("S2")) # S2

otu_S5 <- otu_table %>% select(contains("S5")) # S5

otu_S7 <- otu_table %>% select(contains("S7")) # S7

##############################
## CREATE A PHYLOSEQ OBJECT ##
##############################
# without metadata 
week_Smoins1 <- phyloseq(otu_table(as.matrix(otu_Smoins1), taxa_are_rows = TRUE),
                         tax_table(as.matrix(taxa_table))) # data.frame, not matrix

week_S0 <- phyloseq(otu_table(as.matrix(otu_S0), taxa_are_rows = TRUE),
                    tax_table(as.matrix(taxa_table))) 

week_S2 <- phyloseq(otu_table(as.matrix(otu_S2), taxa_are_rows = TRUE),
                    tax_table(as.matrix(taxa_table))) 

week_S5 <- phyloseq(otu_table(as.matrix(otu_S5), taxa_are_rows = TRUE),
                    tax_table(as.matrix(taxa_table)))

week_S7 <- phyloseq(otu_table(as.matrix(otu_S7), taxa_are_rows = TRUE),
                    tax_table(as.matrix(taxa_table))) 

###################
## SHANNON INDEX ##
###################
## Numeric values of alpha diversity indices for CN and SN
shannon.Smoins1 <- estimate_richness(week_Smoins1, measures = "Shannon")
head(shannon.Smoins1)

shannon.S0 <- estimate_richness(week_S0, measures = "Shannon")
head(shannon.S0)

shannon.S2 <- estimate_richness(week_S2, measures = "Shannon")

shannon.S5 <- estimate_richness(week_S5, measures = "Shannon")

shannon.S7 <- estimate_richness(week_S7, measures = "Shannon")

## Calculate Means of Shannon
Mean.Shannon.Smoins1 <- mean(shannon.Smoins1$Shannon)
print(Mean.Shannon.Smoins1) # result : 5.38979

Mean.Shannon.S0 <- mean(shannon.S0$Shannon)
print(Mean.Shannon.S0) # result : 5.392502

Mean.Shannon.S2 <- mean(shannon.S2$Shannon)
print(Mean.Shannon.S2) # result : 5.780667

Mean.Shannon.S5 <- mean(shannon.S5$Shannon)
print(Mean.Shannon.S5) # result : 5.693753

Mean.Shannon.S7 <- mean(shannon.S7$Shannon)
print(Mean.Shannon.S7) # result : 5.708888

#############
## T-TESTS ##
#############
#To compare Shannon indices depending on weeks
t.test(shannon.S0$Shannon,shannon.S2$Shannon) #S0 vs S2
# result : p-value = 0.0006257 ***

t.test(shannon.S0$Shannon,shannon.S5$Shannon) #S0 vs S5
# result : p-value = 0.00544 **

t.test(shannon.S0$Shannon,shannon.S7$Shannon) #S0 vs S7
# result : p-value = 0.004573 **

t.test(shannon.S2$Shannon,shannon.S5$Shannon) #S2 vs S5
# result : p-value = 0.1636

t.test(shannon.S2$Shannon,shannon.S7$Shannon) #S2 vs S7
# result : p-value = 0.2897

t.test(shannon.S5$Shannon,shannon.S7$Shannon) #S5 vs S7
# result : p-value = 0.8358

## new t-tests with "Smoins1"
t.test(shannon.Smoins1$Shannon,shannon.S0$Shannon) #Smoins1 vs S0
# result : p-value = 0.9855

t.test(shannon.Smoins1$Shannon,shannon.S2$Shannon) #Smoins1 vs S2
# result : p-value = 0.02058 *

t.test(shannon.Smoins1$Shannon,shannon.S5$Shannon) #Smoins1 vs S5
# result : p-value = 0.05216

t.test(shannon.Smoins1$Shannon,shannon.S7$Shannon) #Smoins1 vs S7
# result : p-value = 0.04474 *

#############
## Boxplot ##
#############
## Preparation for boxplots
# Add a column in each data frame of Soil Infestation
shannon.Smoins1$Week <- "Smoins1"
shannon.S0$Week <- "S0"
shannon.S2$Week <- "S2"
shannon.S5$Week <- "S5"
shannon.S7$Week <- "S7"
View(shannon.S0) # two  columns except headers : Shannon index and Week (S0)

# Concatenate tables
alpha_div_week <- bind_rows(shannon.S0, shannon.S2, shannon.S5, shannon.S7) 
View(alpha_div_week) # worked
class(alpha_div_week) # data.frame 

# without S0
alpha_div_week2 <- bind_rows(shannon.S2, shannon.S5, shannon.S7) 

# with all (Smoins1, S0, S2, S5, S7)
alpha_div_week3 <- bind_rows(shannon.Smoins1, shannon.S0, shannon.S2, shannon.S5, shannon.S7) 
View(alpha_div_week3)

## Boxplot
# with S0, S2, S5, S7
boxplot_shannon_week <- ggplot(alpha_div_week, aes(x = "Week", y = Shannon, fill = Week)) + geom_boxplot()
boxplot_shannon_week 

# without S0
boxplot_shannon_week2 <- ggplot(alpha_div_week2, aes(x = "Week", y = Shannon, fill = Week)) + geom_boxplot()
boxplot_shannon_week2 

# with all
boxplot_shannon_week3 <- ggplot(alpha_div_week3, aes(x = "Week", y = Shannon, fill = Week)) + geom_boxplot()
boxplot_shannon_week3



##### III. Beta-diversity #####
#   A. Centered-Log Ratio (CLR) Transformation (normalization)

# CLR_centered_log_ratio.R script

library(compositions)

# Remove OTU.ID (if not : OTU.ID will be included in the clr transformation)
otu_tab2 <- otu_table %>% select(-contains("OTU.ID")) 
otu_clr <- clr(otu_tab2) # on the entire otu table (except the Smoins1)
View(otu_clr)

# Add the OTU.ID column as before 
otu_column <- otu_table %>% select(contains("OTU.ID")) # View(otu_column)
otu_clr <- cbind(otu_column, otu_clr)
View(otu_clr)

## Save file
write.csv(otu_clr,'./data_modified/otu_clr.csv')
## This file was then used for UMAP (Python script)

#   B. PCoA before rarefaction

## Ordination from the soil phyloseq object
soil.ord <- ordinate(soil, method = "PCoA", distance = "bray")

## plot_ordination
# PCoA depending on Soil Infestation
plot_soil <- plot_ordination(soil,soil.ord)
plot_soil2 <- plot_ordination(soil, soil.ord, color = "Soil.Infestation", title = "PCoA + BC, from otu_table")
plot_soil2

# Add ellipses
plot_soil3 <- plot_soil2 + stat_ellipse()  # basic ones
plot_soil3
plot_soil4 <- plot_soil2 + geom_point() + stat_ellipse(type = "norm", linetype = 2) + stat_ellipse(type = "t")
plot_soil4

# PCoA depending on Weeks
plot_soil_week <- plot_ordination(soil, soil.ord, color = "Week", title = "PCoA + BC, from otu_table")
plot_soil_week

# Add ellipses
plot_soil_week2 <- plot_soil_week + stat_ellipse()
plot_soil_week2
# Not necessary in this case
# plot_soil_week3 <- plot_soil_week + geom_point() + stat_ellipse(type = "norm", linetype = 2) + stat_ellipse(type = "t")


#   C. PCoA (with rarefaction : didn't work with otu_clr table)

## using phyloseq object and all data : otu, taxa, meta
# phyloseq_rarefaction_PCoA.R script

# Create the phyloseq object with all data
meta_table <- sample_data(Tmetadata) # here : only solution to be able to include metadata into "soil"
soil <- phyloseq (otu_table(as.matrix(otu_table), taxa_are_rows = TRUE), tax_table(as.matrix(taxa_table)), sample_data(data.frame(meta_table))) # data.frame, not matrix
soil

# If the creation doesn't work from the first time, do :
################################################################################
## LOAD DATA 
# better to do it by your own
otu_table <-read.table('./data_modified/otu_table_16S.csv', sep=';', dec='.', header=T,row.names=1)
meta_table <-read.table("./data_modified/metadata_16S.csv", sep=";", dec=".", header=T,row.names=1)
taxa_table <-read.table("./data_modified/taxa_16S.csv", sep=";", dec=".", header=T,row.names=1)

# Remove columns containing "Smoins1"
otu_tab <- otu_table %>% select(-contains("Smoins1")) # do not forget the "-"
verif_search_Smoins1 <- grep("Smoins1", names(otu_tab), value = TRUE)
verif_search_Smoins1    # Worked : character(0)

# and reload the creation step 
################################################################################

#################
## RAREFACTION ##
#################
rare_soil <- rarefy_even_depth(soil, sample.size = min(sample_sums(soil)), rngseed = 711, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Plot rarefaction curves (this step takes some time)
rare_soil_plot <- ggrare(rare_soil, step = 100, color = "Soil.Infestation", se = TRUE)
# Note : "se = TRUE" to have standard errors

# depending on weeks : 0/2/5/7
rare_week_plot <- ggrare(rare_soil, step = 100, color = "Week", se = TRUE)

# Custom in panels
rare <- rare_soil_plot + facet_wrap(~Soil.Infestation, ncol = 2) + theme_bw()
plot(rare)
rare2 <- rare_week_plot + facet_wrap(~Week, ncol = 4) + theme_bw()
plot(rare2)


#### Ordination and PCoA
## Ordination 
rare_soil.ord <- ordinate(rare_soil, method = "PCoA", distance = "bray")

## Plot depending on Soil Infestation 
plot_rare_soil <- plot_ordination(rare_soil, rare_soil.ord, color = "Soil.Infestation", title = "PCoA + BC, after rarefaction")
plot_rare_soil

## Plot depending on Weeks
plot_rare_week <- plot_ordination(rare_soil, rare_soil.ord, color = "Week", title = "PCoA + BC, after rarefaction")
plot_rare_week

# Add ellipses
plot_rare_soil2 <- plot_rare_soil + stat_ellipse() # only one ellipse
plot_rare_soil2
plot_rare_week2 <- plot_rare_week + stat_ellipse() 
plot_rare_week2 # one is enough here

plot_rare_soil3 <- plot_rare_soil + geom_point() + stat_ellipse(type = "norm", linetype = 2) + stat_ellipse(type = "t")
plot_rare_soil3 # better with 2 here



##### IV. Prevalence ##### 
# prevalence_otu.R script

### Gaol : Apply a prevalence filter before to perform FlashWeave ###

#   A. On soil phyloseq object

## Trim data 
# First : condition is true if a taxa has at least 5 positive counts (across samples) 
condition <-function(x) { sum(x > 0) >= 5 } 
taxaToKeep <-filter_taxa(soil, condition) 
taxa_kept <- prune_taxa(taxaToKeep, soil)
taxa_kept # result : 1199 taxa and 66 samples by 9 taxonomic ranks

## SUMMARY OF CUT-OFF TESTS
# Condition at least 20 -> result : 438 taxa and 66 samples by 9 taxonomic ranks
# Condition at least 10 -> result : 777 taxa and 66 samples by 9 taxonomic ranks
# Condition at least 3 -> result : 1667 taxa and 66 samples by 9 taxonomic ranks
# Condition at least 2 -> result : 2203 taxa and 66 samples by 9 taxonomic ranks
# Condition at least 1 -> result : 14700 taxa and 66 samples by 9 taxonomic ranks
# Last test : 
# Condition at least 66 -> result : 18 taxa and 66 samples by 9 taxonomic ranks


#     B. On otu_clr table  ("soil_clr_tax_met" phyloseq object)

## Create phyloseq object
# otu_clr + taxa + meta
soil_clr_tax_met <- phyloseq (otu_table(as.matrix(otu_clr), taxa_are_rows = TRUE), tax_table(as.matrix(taxa_table)), sample_names(meta_table)) 
soil_clr_tax_met

# If the creation doesn't work from the first time, do :
################################################################################
## Go back to the clr step : 
# Remove OTU.ID (if not : OTU.ID will be included in the clr transformation)
otu_tab2 <- otu_table %>% select(-contains("OTU.ID")) 
otu_clr <- clr(otu_tab2) # on the entire otu table (except the Smoins1)
View(otu_clr)

# Add the OTU.ID column as before 
otu_column <- otu_table %>% select(contains("OTU.ID")) # View(otu_column)
otu_clr <- cbind(otu_column, otu_clr)
View(otu_clr)

# Remove columns containing "Smoins1"
otu_tab <- otu_table %>% select(-contains("Smoins1")) # do not forget the "-"
verif_search_Smoins1 <- grep("Smoins1", names(otu_tab), value = TRUE)
verif_search_Smoins1    # Worked : character(0)

# and reload the creation step 
################################################################################

## Create the cut-off function
# want a data frame where all OTUs are found at least in 10% of all samples
condition3 <-function(x) { sum(x > 0) >= 6.6 } # as we have 66 samples 
taxaToKeep3 <-filter_taxa(soil_clr_tax_met, condition) 
taxa_kept3 <- prune_taxa(taxaToKeep3, soil_clr_tax_met)
taxa_kept3  # results : 847 taxa

# Extract these 847 taxa
OTU3 = as(otu_table(taxa_kept3), "matrix")
if(taxa_are_rows(taxa_kept3)){OTU3 <- t(OTU3)}
# Coerce to data.frame
df_847 = as.data.frame(OTU3)
View(df_847)
df_847_t <- t(df_847) # transpose

taxa3 = as(tax_table(taxa_kept3), "matrix")
if(taxa_are_rows(taxa_kept3)){taxa3 <- t(taxa3)}
# Coerce to data.frame
df_847_taxa = as.data.frame(taxa3)
df_847_taxa_t <- t(df_847_taxa) # transpose
View(df_847_taxa_t)

# Save this file 
write.csv(df_847_t,'./results/testFW/prevalence_clr_847otu.csv')
write.csv(df_847_taxa_t,'./results/testFW/prevalence_clr_847taxa.csv')
