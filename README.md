# Analysis_ROS_burst_D36E
Analysing how ROS burst (in spinach) are altered by single effectors using a D36E effector-delivery-platform

- The repository contains scripts required to analyse the output of a luminometer using a 96-well-plate format (format is shown in raw data files)

- The scripts are to be run in order from Step_1 to Step_3

- There are two alternative scripts for every step. One pipeline normalises the data based on the control mean, the other pipeline based on the control median. 

- data wrangling and normalisation uses python

- data visualisation and statistical analysis is carried out in R

- Raw data files: are named by data_ROS_plate_processed - files essentially contain raw data, however, the first rows of the file were removed in order to input into the downstream analyisis pipeline (hence _processed)
Note: not all plates are ultimately included into the visualisation due to assay repeats/failed assays (information in .R scripts)

- Meta data:      meta data are collected in the MM20221209_summary_ROS_DC3000_ef_library.xlsx file. This file contains all information required to analyse the raw data. 
