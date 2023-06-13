# load libraries
library(ggplot2)
library(dplyr)
library(plotly)
library("DescTools")

# clean environment
rm(list=ls())

# define wd
setwd("C:/Users/Mende012/Documents/Bioinformatics/ROS assays/DR_final_pipeline")

# import file
df <- read.csv("./Step_2_MM20230106_2_auc_data.csv")

# Accidentally named HopN1 as AvrN1 - change to correct name
df$treatment_id[df$treatment_id == 'D36E+AvrN1_Bacterial PAMPs'] <- 'D36E+HopN1_Bacterial PAMPs'

# select only one assays per treatment combination
df <- df[!(df$assay_id == "plate025"),]
df <- df[!(df$assay_id == "plate020"),]
df <- df[!(df$assay_id == "plate029"),]
df <- df[!(df$assay_id == "plate032"),]
df <- df[!(df$assay_id == "plate033"),]
df <- df[!(df$assay_id == "plate034"),]
df <- df[!(df$assay_id == "plate037"),]

df <- df[!((df$assay_id == "plate028") & (grepl(("AvrPto"), df$treatment_id))),]

#################################################################################
#################################################################################

# ANOVA after normalisation
# 1. ANOVA to analyse wheather there are differences between the groups
# H0 All means are the same 
# H1 at least one mean is different
#anova_ROS <- aov(area_under_curve ~  treatment_id, data = df)
#{summary(anova_ROS)
#  summary.lm(anova_ROS)}

# Dunnett's test - summarised analysis

#Dunnet <- DunnettTest(area_under_curve ~ treatment_id, data = df, control = "D36E_Bacterial PAMPs", conf.level = 0.95)
#print(Dunnet)

################################################################################
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

#################################################################################
###############################boxplot pre-scaling###############################
#################################################################################

# colour scheme
cbp1 <- c("#999999", "#44AA99", "#E69F00", "#661100", "#56B4E9", "#000000", "#117733", "#332288", "#F0E442", "#0072B2", "#D55E00", "#999933", "#CC79A7")
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

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
  
  # ANOVA/TUKEY test labels
  #geom_text(data = simple, aes(y= dunn_res.sample, x = dunn_res.ast)) 
  #vjust= -10 ,hjust= 0.5) 


  # significance for single comparisons
  #geom_signif(comparisons = list(c("D36", "D36+AvrE1"),
  #                              c("D36", "D36+HopM1"),
  #                             c("D36", "D36+HopU1"),
  #                            c("D36", "DC3000")),
  #        map_signif_level = TRUE,
#       y_position = c(3, 3.5, 4, 4.5))+


  #ggsignif::geom_signif(
    
  #comparisons = list(c("D36E_Bacterial PAMPs", "D36E+AvrE1_Bacterial PAMPs"),
  #                   c("D36E_Bacterial PAMPs", "D36E+HopAD1_Bacterial PAMPs"),
  #                   c("D36E_Bacterial PAMPs", "D36E+HopAF1_Bacterial PAMPs"),
   #                  c("D36E_Bacterial PAMPs", "D36E+HopAM1_Bacterial PAMPs"),
    #                 c("D36E_Bacterial PAMPs", "D36E+HopF2_Bacterial PAMPs"),
     #                c("D36E_Bacterial PAMPs", "D36E+HopM1_Bacterial PAMPs"),
      #               c("D36E_Bacterial PAMPs", "dCEL_Bacterial PAMPs")
       #               ),
  #annotations = c("***", "***", "***", "***", "***", "*", '***'),
  #y_position = c(17, 19, 21, 23, 25, 27, 29) 
  #)


ggsave("MM20230106_ROS_boxpl_single_mean.svg",width = 10, height = 6)
g2
dev.off()
g2
ggplotly(g2)


###################################################################################
######################## scale by positive control - ##############################
###################################################################################

### STEP1: write function to show  ratio between mean/mean positive control #####
### and sample measurement ###                              c

## write a function to identify the mean/average of the positive control D36+PAMPs

# identify unique assayIDs and sort them to ensure order
assayIDs <- unique(df$assay_id)
assayIDs <- sort(assayIDs)

# creat/test function
med_pos <- function(df, assayID, positive_control){
  a <- as.data.frame(df[df$assay_id == assayID,]) 
  b <- as.data.frame(a[a$treatment_id == positive_control,])
  contr_mean <- mean(b$area_under_curve) # change here from mean to mean
  return(contr_mean)
}

# apply the function over all the assays 
contr_mean <-   lapply(assayIDs, function(x){
  a <- as.data.frame(df[df$assay_id == x,]) 
  b <- as.data.frame(a[a$treatment_id == "D36E_Bacterial PAMPs",])
  contr_mean <- mean(b$area_under_curve)
  return(contr_mean)
})     

## make small df to store mean/mean results per assay
df_ctrmed <- data.frame(assayIDs, unlist(contr_mean))
colnames(df_ctrmed) <- c("assay_id", "mean_positive_control") #change mean/mean
rm(ls = contr_mean)

## write function (function works with single assay)
mult <- function(assayID){
    rows <- nrow(df[df$assay_id == assayID,]) # count rows for of each assay
  a <- df_ctrmed[df_ctrmed$assay_id == assayID,] # access the mean/mean of a respective control
  re_assay_id <- unlist(rep(a$assay_id, times = rows)) # times assayID as many times as there are rows for the assayID
  re_mean <- unlist(rep(a$mean_positive_control, times = rows)) #update mean/mean, times control/mean as many times as there are rows for the assayID
  b <- data.frame(unlist(re_assay_id), unlist(re_mean)) # make a df of length = assay rows
  return(b)
}

# use the function with lapply to repeat the same for all plates
rows <-   lapply(assayIDs, function(x){
  
  rows <- nrow(df[df$assay_id == x,])
  a <- df_ctrmed[df_ctrmed$assay_id == x,]
  re_assay_id <- unlist(rep(a$assay_id, times = rows))
  re_mean <- unlist(rep(a$mean_positive_control, times = rows))
  b <- data.frame(unlist(re_assay_id), unlist(re_mean))
  return(b)
  }
  ) %>% bind_rows()

## reorganise and merge df
# rename columns in a meaningful way 
colnames(rows) <- c("assay_id", "mean_positive_control") #change name mean/mean

# order df so that it is also ordered by plateID 
df <- df[order(df$assay_id, decreasing = FALSE),]

# add additional column to original df
df["mean_positive_control"] <- rows$mean_positive_control

# clean up 
rm(ls=df_ctrmed)
rm(ls=rows)


### use the mean/mean to scale the measurements
## write a function to calculate how the measurments compare to the positive control 
## calculate this as a % of the positive_control_mean

prc_pos_avg <- (df$area_under_curve/df$mean_positive_control)*100
df["prc_pos_avg"] <- prc_pos_avg

cbp1 <- cbp1 <- c("#999999", "#44AA99", "#E69F00", "#661100", "#56B4E9", "#000000", "#117733", "#332288", "#F0E442", "#0072B2", "#D55E00", "#999933", "#CC79A7")
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
  ylab("% log2 (normalised area under the curve)") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) + # turn the tick marks on the x axis 45 degree 
  
  
  # define axis limits if needed
  expand_limits( y = c(0, 5)) +
  
  # add statistical information: error bars, statistical results
  stat_boxplot(geom = "errorbar", # Error bars
               width = 0.2) 


ggsave("MM20230106_ROS_boxpl_single_scaled_mean.svg",width = 10, height = 6)
g3
dev.off()

ggplotly(g3)

#################################################################################
#################################################################################


# Dunnett's test - per plate with scaled data

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

