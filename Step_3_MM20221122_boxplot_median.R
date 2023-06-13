# load libraries
library(ggplot2)
library(dplyr)
library(plotly)
library("DescTools")
library("openxlsx")
library(multcompView)

# clean environment
rm(list=ls())

# define wd
setwd("C:/Users/Mende012/Documents/Bioinformatics/ROS assays/DR_final_pipeline")

# import file
df <- read.csv("./Step_2_MM20230106_2_auc_data.csv")

# Output Prefix
# generate prefix string including date/ Initials
current_date <- gsub("-", "", as.character(Sys.Date()))# Get the current date as character
pre <- paste("MM", current_date, sep = "")

# Accidentally named HopN1 as AvrN1 - change to correct name
#df$treatment_id[df$treatment_id == 'D36E+AvrN1_Bacterial PAMPs'] <- 'D36E+HopN1_Bacterial PAMPs'

# select one assays per treatment combination
df <- df[!(df$assay_id == "plate029"),]
df <- df[!(df$assay_id == "plate030"),]
df <- df[!(df$assay_id == "plate032"),]
df <- df[!(df$assay_id == "plate033"),]
df <- df[!(df$assay_id == "plate034"),]
df <- df[!(df$assay_id == "plate037"),]

# df <- df[!((df$assay_id == "plate028") & (grepl(("AvrPto"), df$treatment_id))),]

#################################################################################
#################################################################################

# ANOVA after normalisation
# 1. ANOVA to analyse wheather there are differences between the groups
# H0 All medians are the same 
# H1 at least one median is different
#anova_ROS <- aov(area_under_curve ~  treatment_id, data = df)
#{summary(anova_ROS)
#  summary.lm(anova_ROS)}

# Dunnett's test - summarised analysis

#Dunnet <- DunnettTest(area_under_curve ~ treatment_id, data = df, control = "D36E_Bacterial PAMPs", conf.level = 0.95)
#print(Dunnet)

################################################################################
###############################    Statistics    ###############################
################################################################################

# Dunnett's test - per plate --> makes more sense due to inter-assay variation
# tbl_df <- df %>% as_tibble() %>% group_by(assay_id)

dunn_res <- lapply(unique(df$assay_id), function(x) as.data.frame(DunnettTest(area_under_curve ~ treatment_id, 
                                                    data = df[df$assay_id==x,], 
                                                    control = "D36E_Bacterial PAMPs", 
                                                    conf.level = 0.95)$`D36E_Bacterial PAMPs`)) %>% bind_rows()

# add column with asterix for significant results
dunn_res$ast <- ''
dunn_res$ast[dunn_res$pval<0.1] <- '.'
dunn_res$ast[dunn_res$pval<0.05] <- '*'
dunn_res$ast[dunn_res$pval<0.01] <- '**'
dunn_res$ast[dunn_res$pval<0.001] <- '***'


dunn_res_1 <- dunn_res
dunn_res$sample <- row.names(dunn_res_1)
rm(dunn_res_1)

# remove comparisons present in every assay 
remove <- c("D36E_MQ-D36E_Bacterial PAMPs","DC3000_Bacterial PAMPs-D36E_Bacterial PAMPs", "DC3000_MQ-D36E_Bacterial PAMPs")
dunn_res <- dunn_res[!grepl(paste(remove, collapse = "|"), dunn_res$sample),]

# sub string sample name to match names of treatments
dunn_res$sample <- sub("-D36E.*", "", dunn_res$sample)

# simple table 
simple <- data.frame(dunn_res$sample, dunn_res$ast, row.names = NULL)

# # safe res
sink(paste(pre, sep = "","ROS_supression_DC3000_library_Dunnet_prescale.txt"))
print(dunn_res)
sink()

###################################################################################
######################## scale by positive control - ##############################
###################################################################################

### STEP1: write function to show  ratio between median/median positive control #####
### and sample measurement ###                              c

## write a function to identify the median/average of the positive control D36+PAMPs

# identify unique assayIDs and sort them to ensure order
assayIDs <- unique(df$assay_id)
assayIDs <- sort(assayIDs)

# creat/test function to find median/mean
med_pos <- function(df, assayID, positive_control){
  a <- as.data.frame(df[df$assay_id == assayID,]) 
  b <- as.data.frame(a[a$treatment_id == positive_control,])
  contr_median <- median(b$area_under_curve) # change here from median to median
  return(contr_median)
}

# apply the function over all the assays 
contr_median <-   lapply(assayIDs, function(x){
  a <- as.data.frame(df[df$assay_id == x,]) 
  b <- as.data.frame(a[a$treatment_id == "D36E_Bacterial PAMPs",])
  contr_median <- median(b$area_under_curve)
  return(contr_median)
})     

## make small df to store median/median results per assay
df_ctrmed <- data.frame(assayIDs, unlist(contr_median))
colnames(df_ctrmed) <- c("assay_id", "median_positive_control") #change median/median
#rm(ls = contr_median)

## write function (function works with single assay) 
mult <- function(assayID){
    rows <- nrow(df[df$assay_id == assayID,]) # count rows for of each assay
  a <- df_ctrmed[df_ctrmed$assay_id == assayID,] # access the median/median of a respective control
  re_assay_id <- unlist(rep(a$assay_id, times = rows)) # times assayID as many times as there are rows for the assayID
  re_median <- unlist(rep(a$median_positive_control, times = rows)) #update median/median, times control/median as many times as there are rows for the assayID
  b <- data.frame(unlist(re_assay_id), unlist(re_median)) # make a df of length = assay rows
  return(b)
}

# use the function with lapply to repeat the same for all plates
rows <-   lapply(assayIDs, function(x){
  
  rows <- nrow(df[df$assay_id == x,])
  a <- df_ctrmed[df_ctrmed$assay_id == x,]
  re_assay_id <- unlist(rep(a$assay_id, times = rows))
  re_median <- unlist(rep(a$median_positive_control, times = rows))
  b <- data.frame(unlist(re_assay_id), unlist(re_median))
  return(b)
  }
  ) %>% bind_rows()

## reorganise and merge df
# rename columns in a medianingful way 
colnames(rows) <- c("assay_id", "median_positive_control") #change name median/median

# order df so that it is also ordered by plateID 
df <- df[order(df$assay_id, decreasing = FALSE),]

# add additional column to original df
df["median_positive_control"] <- rows$median_positive_control

# clean up 
rm(ls=df_ctrmed)
rm(ls=rows)


### use the median/median to scale the measurements
## write a function to calculate how the measurments compare to the positive control 
## calculate this as a % of the positive_control_median

prc_pos_avg <- (df$area_under_curve/df$median_positive_control)*100
df["prc_pos_avg"] <- prc_pos_avg


#################################################################################
#########################  export for summary heatmap ###########################
#################################################################################

e_line <- unique(df$treatment_id)
#e_line <- lapply(e_line,sort,decreasing=TRUE)

# apply contr_median function to find the median of each effector line percentage
sample_median <-   lapply(e_line, function(x){
  a <- as.data.frame(df[df$treatment_id == x,]) 
  med <- median(a$prc_pos_avg)
  return(med)
})

e_line <- as.character(e_line)
d2 <- cbind.data.frame("Effector" = e_line, "Median" = unlist(sample_median))
d2$Effector <- sub("*_Bacterial PAMPs", "", d2$Effector)

write.csv(d2, "./summary_median_ROS_percentage.xlsx")

head(df)
#################################################################################
##########################  Statistics post-scaling   ###########################
#################################################################################

# Dunnett's test - per plate with scaled data ----

dunn_scal <- lapply(unique(df$assay_id), function(x) as.data.frame(DunnettTest(prc_pos_avg ~ treatment_id, 
                                                                               data = df[df$assay_id==x,], 
                                                                               control = "D36E_Bacterial PAMPs", 
                                                                               conf.level = 0.95)$`D36E_Bacterial PAMPs`)) %>% bind_rows()

# add column with asterix for significant results
dunn_scal$ast <- ''
dunn_scal$ast[dunn_scal$pval<0.1] <- '.'
dunn_scal$ast[dunn_scal$pval<0.05] <- '*'
dunn_scal$ast[dunn_scal$pval<0.01] <- '**'
dunn_scal$ast[dunn_scal$pval<0.001] <- '***'


dunn_scal_1 <- dunn_scal
dunn_scal$sample <- row.names(dunn_scal_1)
rm(dunn_scal_1)

# remove comparisons present in every assay 
remove <- c("D36E_MQ-D36E_Bacterial PAMPs","DC3000_Bacterial PAMPs-D36E_Bacterial PAMPs", "DC3000_MQ-D36E_Bacterial PAMPs")
dunn_scal <- dunn_scal[!grepl(paste(remove, collapse = "|"), dunn_scal$sample),]

# sub string sample name to match names of treatments
dunn_scal$sample <- sub("-D36E.*", "", dunn_scal$sample)

# simple table 
simple_scal <- data.frame(dunn_scal$sample, dunn_scal$ast, row.names = NULL)

# safe res
sink(paste(pre, sep = "","ROS_supression_DC3000_library_Dunnet_scaled.txt"))
print(dunn_scal)
sink()

# ANOVA + TUKEY - per plate with scaled data ----
# The following analysis requires that ther is no "-" present in the sample name
# Replace "-" with "_"

df$treatment_id <- gsub("-", "_", df$treatment_id)

# 4. Statistical analysis across samples using ANOVA+TUKEY ----
# 4.1. ANOVA to analyse whether there are differences between the groups
#      H0 All means are the same 
#      H1 at least one mean is different

anova_prc <- aov(prc_pos_avg ~ treatment_id, data = df)

{summary(anova_prc)
  summary.lm(anova_prc)}

# 4.2 TUKEY Posthoc to see which groups are significantly different 
TUKEY <- TukeyHSD(anova_prc)

# 4.3 Summarize TUKEY test results for saving & plotting

# 4.3.1 write a funktion to generate labels for the test results
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc results
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  # order labels according to boxplot:
  Tukey.labels$Pst_strain = rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$Pst_strain) , ]
  return(Tukey.labels)
}

# 4.3.2 Apply the function on the df
LABELS <- generate_label_df(TUKEY , "treatment_id")
names(LABELS) <- c("Letters","treatment_id")

# 4.3.3 safe TUKEY test results
sink(paste(pre, sep = "","ROS_supression_DC3000_library_TUKEY_scaled.txt"))
print(LABELS)
sink()

# 4.3.4 obtain letter position for y axis using means
yvalue <- aggregate(prc_pos_avg ~ treatment_id, data = df, mean) 
LABELS$treatment_id <- (unlist(lapply(strsplit(LABELS$treatment_id,'[_]'),function(x){x[[1]]})))

# 4.3.5 generate df with results TUKEY test
final <- merge(yvalue,LABELS) 

# 4.3.6 For ploting extract label vector
stat_labels <- as.vector(LABELS$Letters)


# 5. Box Plots comparing Pst strains, 2 dpi, TUKEY results, save as pdf (default) ----

g_x <-  ggplot(df, aes(x=as.factor(treatment_id), y=prc_pos_avg)) + 
  
  geom_boxplot(fill="white",                 # box colour
               outlier.colour = "white",     # Outliers color, 
               alpha=0.9) +                  # Box color transparency
  
  theme_bw() +  # make the bg white
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), # remove background, frame
        axis.line = element_line(colour = "black")) +
  
  xlab("Pseudomonas strain") +                   # axis label
  ylab("log10 (CFU/cm2)") +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # turn the tick marks on the x axis 45 degree +
  #expand_limits( y = c(0, 10)) +
  
  geom_jitter(shape=16, position=position_jitter(0.1),
              aes(colour = factor(treatment_id))) +
  stat_boxplot(geom = "errorbar", # Error bars
               width = 0.2) +
  geom_text(data = final, aes(y= treatment_id, label = Letters), 
            vjust= -7 ,hjust= 0.5)


g_x



#################################################################################
###########################    boxplot pre-scaling    ###########################
################################################################################# ----

# simplify names in df
#df$treatment_id <- sub("D36E", "", df$treatment_id)
df$treatment_id <- sub("*_Bacterial PAMPs", "", df$treatment_id)
#df$treatment_id <- sub("_MQ", "mock", df$treatment_id)

# colour scheme
cbp1 <- c("#999999", "#44AA99", "#E69F00", "#661100", "#56B4E9", "#000000", "#117733", "#332288", "#F0E442", "#0072B2", "#D55E00", "#999933", "#CC79A7", "#E69F00")

#based on: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

# basic boxplot
g2 <-  ggplot(
  
  #input data
  df, aes(x=as.factor(treatment_id), y=area_under_curve)) +
  
  # generate basic boxplot
  geom_boxplot(fill="white",                 # box colour
               outlier.colour = "white",     # Outliers color, 
               alpha=0) +                    # Box color transparency
  
  # overlay with jitter
  geom_jitter(shape=16, position=position_jitter(0.1),
              aes(colour = factor(assay_id))) +
  scale_colour_manual(values = cbp1) +
  
  # define the theme of the boxplot
  theme_bw() +  # make the bg white
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), # remove background, frame
        axis.line = element_line(colour = "black")) +
  
  # label the axises 
  xlab("Pseudomonas strain") +                
  ylab("log2 (normalised area under the curve)") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) + # turn the tick marks on the x axis 45 degree 
  
  
  # define axis limits if needed
  expand_limits( y = c(0, 5)) +
  
  # add statistical information: error bars, statistical results
  stat_boxplot(geom = "errorbar", # Error bars
               width = 0.2) 

ggsave("MM20230202_ROS_boxpl_single_median.svg",width = 10, height = 6)
g2
dev.off()
g2

ggplotly(g2)
ggplotly(g2)


###############################################################################
################################### plot graph ################################
###############################################################################
length(unique(df$assay_id))

cbp1 <- rep("#E69F00", each = 20)
          
#cbp1 <- c("#999999", "#44AA99", "#E69F00", "#661100", "#56B4E9", "#000000", "#117733", "#332288", "#F0E442", "#0072B2", "#D55E00", "#999933", "#CC79A7", "#E69F00")
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

# scaled boxplot
g3 <-  ggplot(
  
  #input data
  df, aes(x=as.factor(treatment_id), y=prc_pos_avg)) +
  
  # generate basic boxplot
  geom_boxplot(fill="white",                 # box colour
               outlier.colour = "white",     # Outliers color, 
               alpha=0) +                    # Box color transparency
  
  
  # overlay with jitter
  geom_jitter(shape=16, position=position_jitter(0.1),
              aes(colour = factor(assay_id))) +
  scale_colour_manual(values = cbp1) +
  
  # define the theme of the boxplot
  theme_bw() +  # make the bg white
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), # remove background, frame
        axis.line = element_line(colour = "black")) +
  
  # label the axises 
  xlab("Pseudomonas strain") +                
  ylab("log2 (normalised area under the curve)") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) + # turn the tick marks on the x axis 45 degree 
  
  
  # define axis limits if needed
  expand_limits( y = c(0, 5)) +
  
  # add statistical information: error bars, statistical results
  stat_boxplot(geom = "errorbar", # Error bars
               width = 0.2) 


ggsave("MM20230202_ROS_boxpl_single_median.svg",width = 10, height = 6)
g3
dev.off()

ggplotly(g3)





# In progress LAST: Safe stat results for better comparison in one df ----
# Sort df
LABELS_sorted <- LABELS[order(LABELS$treatment_id), ]
simple_scal_sorted <- simple_scal[order(simple_scal$dunn_scal.sample), ]

stat_comp <- data.frame(  "Treatment" = LABELS_sorted$treatment_id ,
                          "Not_scaled_Dunnet's" = ,
                          "Scaled_Dunnet's" = simple_scal$dunn_scal.ast  ,
                          "Scaled_TUKEY" = LABELS$Letters)


sink(paste(pre, sep = "","CFU_count_DC3000_library_summary_stats.txt"))
print(stat_comp)
sink()