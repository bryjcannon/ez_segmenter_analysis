if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("R.matlab", "digest", "rlang", "flowCore", "ks", "flowVS", "flowViz", "RColorBrewer", "gtools", "gplots", 
                       "ggplot2", "openxlsx", "samr", "lattice", "flowStats", "gdata", "Rtsne", "umap",
                       "FlowSOM", "dplyr", "plyr", "pryr", "doBy", "scales", "mixOmics", "reshape2", 
                       "plotly", "Rmisc", "Hmisc", "EBImage", "magick", "phonTools"))

library("R.matlab")
library("digest")
library("rlang")
library("flowCore")
library("ks")
library("flowVS") #
library("flowViz")
library("RColorBrewer")
library("gtools")
library("gplots")
library("ggplot2")
library("openxlsx") #
library("samr") #
library("lattice")
library("flowStats") #
library("gdata")
library("Rtsne")
library("umap")
library("FlowSOM") #
library("dplyr")
library('plyr')
library("pryr")
library("doBy") #
library("scales")
library("mixOmics") #
library("reshape2")
library("plotly") #
library("Rmisc")
library("Hmisc") #
library("EBImage")
library("magick")
library("phonTools")


#https://support.bioconductor.org/p/109128/ --> explains why use Biobase::exprs
exprs = Biobase::exprs

# color palette aken fro stackOverflow <- https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", "black"
)

# function to read in object datasets in csv form and attach additional data attributes to the csv's
#read_ez_data  <- function(ez_data, ) {