libs <- c('ggplot2','RColorBrewer','reshape2','devtools',
          'pheatmap','flowCore','plyr','scales','coin','grid',
          'gridExtra','tidyr','Rtsne','FlowSOM','randomcoloR','viridis')

for (L in libs){
  if(!require(L, character.only = T)) {
    install.packages(L, dep = T, quiet = T)
    
    if (!require(L, character.only = T)) {
      if (L %in% c('flowCore', 'flowWorkspace')) {
        install_github("RGLab/flowCore", ref="trunk")
      } else {
        suppressWarnings(BiocManager::install(L))
      }
    }
  }
}

library(flowCore)
library(readr)
setwd("/Volumes/KausaliaHD/R4R")
#source('applyAsinh.R')
outputFolder = "/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellSeg_uci2712J/Plots_Dynamic_Expansion"
PaperADPanel <- read_csv("/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellSeg_uci2712J/190604HiResCA2/info/190503_TMAADPanel.csv")

applyAsinh <- function(value, cofactor = 5) {
  value <- value * 100 - 1
  for(i in 1:length(value)) {
    if((value[i] < 0) | is.na(value[i])) value[i] <- rnorm(1, mean = 0, sd = 0.01)
  }
  value <- value / cofactor
  value <- asinh(value)
  
  return(value)
}

concatFCStoDF <- function(file.names, region, channel.metals) {
  cat("Reading FCS files...")
  cell.data <- data.frame()
  files.FCS <- read.flowSet(file.names, transformation = FALSE, emptyValue = FALSE, truncate_max_range = FALSE)
  cat("Complete\n")
  
  colnames(files.FCS) <- gsub("-|_|\\s+", "", colnames(files.FCS), ignore.case = T)
  
  cat("Applying asinh transformation...")
  files.FCS_transform <- files.FCS[, colnames(files.FCS) %in% channel.metals]
  
  
  files.FCS_nontransform <- files.FCS[, !(colnames(files.FCS) %in% channel.metals)]
  cat("Complete\n")
 
  for (i in 1:length(files.FCS)) {
    cat(paste("Concatenating:", file.names[i], "\n", sep = " "))
    
    file.fcs_transform <- files.FCS_transform[[i]]
    fcs.expr_transform <- as.data.frame(exprs(file.fcs_transform))
    fcs.expr_transform <- apply(fcs.expr_transform,2,applyAsinh)
    
    file.fcs_nontransform <- files.FCS_nontransform[[i]]
    fcs.expr_nontransform <- as.data.frame(exprs(file.fcs_nontransform))
    
    fcs.expr <- cbind(fcs.expr_transform, fcs.expr_nontransform)
    
    #fcs.expr$Region = gsub('190604HiRes', '', region)
    #fcs.expr$Region = gsub('190505HiRes', '', region)
    fcs.expr$Region = region
    fcs.expr$Point = gsub('\\.fcs', '', gsub('dataScaleSizeFCS_p', '', basename(file.names[i])))
    
    cell.data <- rbind(cell.data, fcs.expr)
  }
  cat("Concatenation complete.\n")
  
  return(cell.data)
}

fcsPath = '/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellSeg_uci2712J/'
fcsRuns = c('190604HiResCA2', '190505HiResDG')
fcsPixel = 'fcs_single_cell_dynamic_expansion'
fcsNames = list.files(path = paste0(fcsPath,fcsRuns[1], '/', fcsPixel), pattern = paste0("^dataScaleSizeFCS.*fcs$"), full.names = T)
allFcsFileNames = fcsNames
allDataFcs = concatFCStoDF(fcsNames, fcsRuns[1], PaperADPanel$Label)

fcsNames = list.files(path = paste0(fcsPath,fcsRuns[2], '/', fcsPixel), pattern = paste0("^dataScaleSizeFCS.*fcs$"), full.names = T)
allFcsFileNames = c(allFcsFileNames, fcsNames)
allDataFcs = rbind(allDataFcs, concatFCStoDF(fcsNames, fcsRuns[2], PaperADPanel$Label))

library(FlowSOM)
library(Rtsne)
seed = 140214

numClusters = 7
#palette <- distinctColorPalette(numClusters)
#palette[1] = ""
palette = c("#FF0000", "#1ABC9C","#800080","#0000ff", "#FFA500","#2980b9","#c0392b", "#00ffff","#00BFFF","#ec7063", "#ff00ff","#ffff00","#0000FF","#DC143C","#FFFFFF")
palette = palette[1:numClusters]
palette = palette[c(1,2,3,5,4)]
exclude = c("C12", "Na23", "Si28", "Ca40", "Background", "empty139", "Ta181", "Au197")
# clusteringMarkers = c("HistoneH3Lyo","Calretinin", "CD31",
#                       "MAP2",
#                       "PanGAD6567",
#                       "VGAT",
#                       "Reelin",
#                       "EEA1",
#                       "pTDP43",
#                       "CD45", "8OHGuano",
#                       "TH",
#                       "MAG",
#                       "VGLUT2",
#                       "Calbindin",
#                       "Iba1",
#                       "MBP",
#                       "MCT1",
#                       "CD33Lyo",
#                       "GFAP",
#                       "CD105" )
clusteringMarkers = c("Calretinin", "Calbindin", "MAP2", "VGAT", "8OHGuano", "GFAP", "Iba1", "CD45", "MCT1", "CD31", "CD105")
markers = setdiff(PaperADPanel$Label, exclude)

fSOM <- FlowSOM(allFcsFileNames,
                compensate = F, transform = T, toTransform = clusteringMarkers, 
                transformFunction = applyAsinh, scale = F,
                colsToUse = clusteringMarkers, nClus = numClusters, seed = seed)

table(fSOM$metaclustering)
metaClustering = fSOM$metaclustering

fSOM_clustering = data.frame(fSOM$FlowSOM$map$mapping)
colnames(fSOM_clustering) = c("Cluster", "Value")
fSOM_clustering$Meta = as.numeric(fSOM$metaclustering[fSOM_clustering[,1]])

### HACK, REMOVE AFTER USE
fSOM_clustering[which(fSOM_clustering$Meta == 4), "Meta"] = 2
fSOM_clustering[which(fSOM_clustering$Meta == 6), "Meta"] = 2
fSOM_clustering[which(fSOM_clustering$Meta == 5), "Meta"] = 4
fSOM_clustering[which(fSOM_clustering$Meta == 7), "Meta"] = 5
numClusters = 5
### END HACK HERE

colnames(allDataFcs) = gsub(' ','', colnames(allDataFcs))

all_data = cbind(allDataFcs[, c(markers, 'Region','Point', 'cellLabelInImage')], fSOM_clustering)

fSOMcodes = fSOM$FlowSOM$map$codes
write.csv (fSOMcodes, paste0(outputFolder,'/fSOM_FlowSOM_map_codes.csv'))
fSOMcodesMean = apply(fSOMcodes,2,mean)
write.csv(fSOMcodesMean, paste0(outputFolder,'/mean_of_fSOM_FlowSOM_map_codes.csv'))

# Make a barplot of cluster sizes
plot_data = data.frame(table(fSOM_clustering$Meta))
ggplot(plot_data, aes(x = Var1, y = Freq)) +
  geom_bar(aes(fill = Var1), stat = "identity", colour = "black", fill = palette) +
  #geom_text(aes(x = Var1, y = (max(plot_data$Freq)/2), label = paste0("Cluster ", Var1, "\n", Freq)), col='black', size = 2.5) +
  xlab("FlowSOM metaclusters") + ylab("Cluster size") +
  theme_bw(15) +
  theme(legend.position="none")
ggsave(filename = paste0(outputFolder,"/FlowSOMClusterSize.pdf"), width = 12, height = 3)

all_data = all_data[, -which(names(all_data) %in% c("Cluster", "Value"))]
all_data_sampled <- data.frame()
unifSize = 10000000
for (i in 1:numClusters) {
  all_data_subset = all_data
  all_data_subset = subset(all_data_subset, Meta == i)
  
  if (unifSize < nrow(all_data_subset)) {
    unif = sample(1:nrow(all_data_subset), size = unifSize)
  } else {
    unif = 1:nrow(all_data_subset)
  }
  unif = all_data_subset[unif,]
  unif = melt(unif, id = c("Meta"))
  colnames(unif) = c("Meta", "Channel", "Values")
  
  all_data_sampled <<- rbind(all_data_sampled, unif)
}

# Violin plots per fSOM cluster for each channel
plot_list = list()
for (channel in markers) {
  plot_data = subset(all_data_sampled, Channel == channel)
  
  # temp_data = data.frame()
  # for (clust in unique(plot_data$Meta)) {
  #   temp_data = rbind(temp_data, subset(plot_data, Meta == clust & Values < quantile(Values, 0.95)))
  # }
  # plot_data = temp_data
  # rm(temp_data)
  
  g = ggplot(plot_data, aes(x = as.factor(Meta), y = as.numeric(Values))) +
    geom_violin(aes(fill = Meta), scale = "width", width = 0.6) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
    xlab("") + ylab("") +
    ggtitle(channel) +
    #ylim(c(0, 10)) +
    theme_bw(15) +
    theme(legend.position = "none")
  plot_list[[channel]] <- g
}

# Plot violin + boxplots for each channel separately. Show all clusters on each plot
picName = paste0(outputFolder, "/fSOM_sample_violins_perChannel.pdf")
cat(paste0("Plotting: ", picName, "\n"))
pdf(picName, width = 40, height = 60, onefile = T)
do.call(grid.arrange, list(grobs = plot_list, ncol = 3))
dev.off()

# Same as above, but all plots are shown as one long column on the figure
# picName = paste0("Results/", fcsName, "__", popName, "_fSOM_sample_singleColumn.pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 8, height = 100, onefile = FALSE)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 1))
# dev.off()

# Violin plots per each channel for all fSOM clusters
plot_list = list()
for (meta in unique(all_data_sampled$Meta)) {
  plot_data = subset(all_data_sampled, Meta == meta)
  plot_data = subset(plot_data, Channel %in% markers)
  
  g = ggplot(plot_data, aes(x = Channel, y = as.numeric(Values))) +
    geom_violin(aes(fill = Channel), scale = "width", width = 0.6) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
    xlab("") + ylab("") +
    ggtitle(paste0("Cluster: ", meta)) +
    #ylim(c(0, 10)) +
    theme_bw(12) +
    theme(legend.position = "none") +
    theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
          axis.title.x = ggplot2::element_text(face = "bold", size = 10),
          axis.title.y = ggplot2::element_text(face = "bold", size = 10),
          axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
  plot_list[[meta]] <- g
}

# Plot violin + boxplots for each channel separately. Show all clusters on each plot
picName = paste0(outputFolder, "/fSOM_cluster_violins_perCluster.pdf")
cat(paste0("Plotting: ", picName, "\n"))
pdf(picName, width = 40, height = 2.5*numClusters, onefile = T)
do.call(grid.arrange, list(grobs = plot_list, ncol = 3))
dev.off()

# Make a heatmap for channels vs clusters with median expression values from all_data. Scaled by marker (row) version of the above
plot_data = all_data[, c(markers, "Meta")]
plot_data = melt(plot_data, id = c("Meta"))
plot_data = dcast(plot_data, variable ~ Meta, median)
plot_data$variable = NULL
plot_data = t(plot_data)
plot_data = apply(plot_data, 2, scale)
plot_data = rescale(plot_data, to = c(0, 2))
colnames(plot_data) = markers
rownames(plot_data) = 1:numClusters
dev.off()
pdf(paste0(outputFolder, "/FlowSOMClustersMediansScaled.pdf"), width = 8, height = 5, onefile = F)
pheatmap(plot_data, cluster_cols = T, cluster_rows = T, scale = "none", main = "FlowSOM clusters, scaled channels",
         color = viridis(length(plot_data),alpha=1, begin=0, end=1, direction=1, option = "D"), cellwidth = 10, cellheight = 10)
dev.off()


