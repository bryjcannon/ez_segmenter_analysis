# TBMIBIconcatenatecsv.R
# Author: Erin McCaffrey
# Date created: 190116
# Overview: This script reads in the csv for raw, cell size normalized, and asinh-transformed 
# cell size normalized expression data for all points in the study, appends a column with the 
# sample ID number and tissue type, then concatenates dataframes based on sample groupings of 
# data type. Exports as csv files for later use.

##..Set wd and load packages/libraries..##

#setwd("~/Desktop/MIBIProjects/Human_ATB_paper-cohort/181029_HumanATB_cohort/all-fields/denoised/dataPerCell")
setwd("~/Desktop/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")
require(gtools)
library(dplyr)

##..Assign tissue type to sample IDs..##

imm_tons<-c(2,10,17,24,31,38,45,51,56,78,83,86) 
imm_spl<-c(1,9,16,23,30,37,44,62)
imm_plac<-c(3,11,18,25,32,39,46,63)
granA_lung<-c(64,65,21,84,42,88,28,89,85,13,35,36,14,15,6,7) 
granA_MAC<-c(57,58,19,87)
granA_pleura<-c(33,34,26,27,40,61)
granA_endo<-c(47,48)
granA_LN<-c(54,55)
gran_BC<-c(79,80,81,82,52,53)

##..Read in all raw, size-normed, and size-normed+asinh-transformed csv..##

#raw
files_raw<-list.files(pattern='dataFCS_.*\\.csv') 
files_raw<-mixedsort(files_raw) 
data_raw<-lapply(files_raw,read.csv)

#cs normed
files_sizenorm<-list.files(pattern='dataScaleSizeFCS_.*\\.csv')
files_sizenorm<-mixedsort(files_sizenorm) 
data_sizenorm<-lapply(files_sizenorm,read.csv)

#asinh(cs normed)
files_sizenorm_asinh<-list.files(pattern='dataScaleSizeTransFCS_.*\\.csv')
files_sizenorm_asinh<-mixedsort(files_sizenorm_asinh) 
data_sizenorm_asinh<-lapply(files_sizenorm_asinh,read.csv)

##..Create concetenated matrices..##

#create dataframes of all with sample ID
all_raw<-bind_rows(data_raw, .id="SampleID") #raw

all_sizenorm<-bind_rows(data_sizenorm, .id="SampleID") #normalized by cell size

all_sizenorm_scale<-all_sizenorm #scaled by 1000 
all_sizenorm_scale[,4:50]<-all_sizenorm_scale[,4:50]*1000

all_sizenorm_asinh<-bind_rows(data_sizenorm_asinh, .id="SampleID")

#add tissue type as label

all_raw<-all_raw %>% mutate(Tissue=case_when(all_raw$SampleID %in% imm_tons ~ "tonsil",
                                          all_raw$SampleID %in% imm_spl ~ "spleen",
                                          all_raw$SampleID %in% imm_plac ~ "placenta",
                                          all_raw$SampleID %in% granA_lung ~ "gran_lung",
                                          all_raw$SampleID %in% granA_MAC ~ "gran_MAC",
                                          all_raw$SampleID %in% granA_pleura ~ "gran_pleura",
                                          all_raw$SampleID %in% granA_endo ~ "gran_endo",
                                          all_raw$SampleID %in% granA_LN ~ "gran_LN",
                                          all_raw$SampleID %in% gran_BC ~ "gran_lung_BC"
                                          ))

all_sizenorm$Tissue<-all_raw$Tissue
all_sizenorm_scale$Tissue<-all_raw$Tissue
all_sizenorm_asinh$Tissue<-all_raw$Tissue

#add internal Sample ID to the mycobacterial granulomas

all_raw<-all_raw %>% mutate(PatientID=case_when(all_raw$SampleID %in% c(64,65) ~ 1,
                                                all_raw$SampleID %in% c(21,84) ~ 2,
                                                all_raw$SampleID %in% c(42,88) ~ 3,
                                                all_raw$SampleID %in% c(28,89) ~ 4,
                                                all_raw$SampleID %in% c(85,13) ~ 11,
                                                all_raw$SampleID %in% c(35,36) ~ 34,
                                                all_raw$SampleID %in% c(14,15) ~ 17,
                                                all_raw$SampleID %in% c(57,58) ~ 21,
                                                all_raw$SampleID %in% c(19,87) ~ 23,
                                                all_raw$SampleID %in% c(6,7) ~ 30,
                                                all_raw$SampleID %in% c(33,34) ~ 18,
                                                all_raw$SampleID %in% c(26,27) ~ 20,
                                                all_raw$SampleID %in% c(40,61) ~ 31,
                                                all_raw$SampleID %in% c(47,48) ~ 29,
                                                all_raw$SampleID %in% c(54,55) ~ 12,
                                                all_raw$SampleID %in% c(79,80) ~ 8,
                                                all_raw$SampleID %in% c(81,82) ~ 6,
                                                all_raw$SampleID %in% c(52,53) ~ 41
))

all_sizenorm$PatientID<-all_raw$PatientID
all_sizenorm_scale$PatientID<-all_raw$PatientID
all_sizenorm_asinh$PatientID<-all_raw$PatientID

##..Save concatenated dataframes..##
write.csv(all_raw, file="allsamples_data3px_annotated.csv",row.names = FALSE)
write.csv(all_sizenorm, file="allsamples_dataCS3px_annotated.csv",row.names = FALSE)
write.csv(all_sizenorm_asinh, file="allsamples_dataCSasinh3px_annotated.csv",row.names = FALSE)
write.csv(all_sizenorm_scale, file="allsamples_dataCS3px_scale.csv",row.names = FALSE)
