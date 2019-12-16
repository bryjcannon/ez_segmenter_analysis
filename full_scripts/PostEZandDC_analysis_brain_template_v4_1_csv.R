# First Brain MIBI <- [ez_segmenter + deepCell] Paper Data Analysis Script
# Working with post-ez_segmenter data (MATLAB generated) for single object analysis
# Author: Bryan Cannon 2019 (multiple code snippets taken or built from DM, FH, EFM, DT)

##### STEPS #####:

##### 0) LOAD PKGS #####
  
  # load necessary packages
  # go to packages folder, this line will be converted into a function later

  # set seed for downstream analysis
  seed <- 123 
  set.seed(seed)

##### 1a) LOAD ezSeg DATA #####
  # move to location of your data, load in csv and .mat files (spatial information stored in .mat)
  ez_folder <- '/Volumes/BryJC_Stanford/For_Ez_Segmentor/HiADCase_Hippocampus /denoisedfft_HiResADuci2717J'
  setwd(ez_folder)
  # specify run folders
  ezRun <- c('ezSegResults_CA2+', 'ezSegResults_DG+')
  ezRunsShort = c('CA2', 'DG')
  object_types <- c('amyloidopathy', 'tauopathy', 'microglia_process', 'vessel_CD31_CD105', 'vessel_MCT1')
  
# create data.frame containing single object data from csv files - for each run and each object type
  master_obj_data <- data.frame()
  
  for (run_index in 1:length(ezRun)){
    obj_data_raw_all <- data.frame()
    
    for (obj_type in object_types) {
      csv_names <- list.files(path = paste0(ezRun[run_index],'/','objects_points'), recursive = T, full.names = T, pattern = paste0(obj_type, "_dataScaleSize.csv")) # read in csv files for object type
      csv_names <- mixedsort(csv_names)
      
      obj_data_raw <- lapply(csv_names, read.csv) %>% bind_rows() # grab data from csv's then convert data to data.frame
      
      obj_type_id <- rep(obj_type, dim(obj_data_raw)[1]) # create column of length object number with obj_type info
      obj_data_raw <- cbind(obj_data_raw, obj_type_id) # add obj_type_id to data
      
      obj_data_raw_all <- rbind(obj_data_raw_all, obj_data_raw) # collate to run data.frame
    }
    
    run_type_id <- rep(ezRunsShort[run_index], dim(obj_data_raw_all)[1]) # create column of length object number with run_type info, name run names to reflect actual regions
    obj_data_raw_all <- cbind(obj_data_raw_all, run_type_id) # add run_type_id to data
    
    master_obj_data <- rbind(master_obj_data, obj_data_raw_all) # collate to master data.frame
    
    rm(obj_data_raw)
    rm(obj_data_raw_all)
  }

##### 1b) LOAD Deep Cell DATA #####
  # move to location of your data, load in fcs and .mat files (spatial information retained here)
  data_folder_DC <- '/Volumes/BryJC_Stanford/deepCell/'
  setwd(data_folder_DC)
  deepcellRun <- c('190505HiResDG', '190604HiResCA2')
  deepcellRunsShort = c('CA2', 'DG')
  
  # create FlowSet from fcs files
  
  master_cell_data <- data.frame()
  
  for (run_index in 1:length(deepcellRun)) {
    fcs_names <- c(list.files(path = paste0(deepcellRun[run_index],'/','fcs_single_cell_dynamic_expansion'), full.names = T, pattern = "(dataScaleSizeFCS.*fcs)"))
    cell_flowSet <- read.flowSet(files = fcs_names, alter.names = T, transformation = FALSE, emptyValue = FALSE, truncate_max_range = FALSE)
    cell_data_matrix <- fsApply(cell_flowSet, Biobase::exprs) #https://support.bioconductor.org/p/109128/ --> explains why use Biobase::exprs
    cell_data_raw <- as.data.frame(cell_data_matrix)
    
    # add columns to each file denoting point number, type of tissue, name of objects
    point_names <- sampleNames(cell_flowSet)
    point_numbers <- as.numeric(gsub("[^0-9]*", '', gsub(".*p[^0-9]*", '', point_names))) # grabs point numbers from filenames and converts to numeric
    point_id <- as.vector(x = NULL)  # create a list 'point_source' that assigns a point_id to each object
    
    dim = fsApply(cell_flowSet, dim)
    dim = as.numeric(dim[,1])
    for(i in 1:length(dim)) { #loop creates the actual vector with point id's for each object
      temp_point_id <- rep(point_numbers[i], dim[i])
      point_id <- c(point_id, temp_point_id)
    }
    cell_data_raw <- cbind(cell_data_raw, point_id) # add point_id to data
    
    cell_type_id <- rep('cell', dim(cell_data_raw)[1]) # create column of length object number with cell_type info (just 'cell' for now)
    cell_data_raw <- cbind(cell_data_raw, cell_type_id) # add cell_type_id to data
    
    run_type_id <- rep(deepcellRunsShort[run_index], dim(cell_data_raw)[1]) # create column of length object number with run_type info, name run names to reflect actual regions
    cell_data_raw <- cbind(cell_data_raw, run_type_id) # add run_type_id to data
    
    master_cell_data <- rbind(master_cell_data, cell_data_raw)
    
    rm(cell_data_matrix)
    rm(cell_flowSet)
    rm(cell_data_raw)
  }

##### 1c) Combine Deep Cell and ezSeg Data #####
  # rename deepCell labels to match ez_segmenter labels
  master_cell_data <- plyr::rename(master_cell_data, c('cellLabelInImage' = 'obj_id', 'cellSize' = 'obj_size', 'cell_type_id' = 'obj_type_id')) # rename deep cell columns to match ez names
  master_cell_data <- master_cell_data %>% select(point_id, everything()) # move point_id to start of deep cell data to match ez column order
  
  # row bind the deepCell dataFrame data with ezSeg dataFrame data
  master_data <- rbind(master_cell_data, master_obj_data[,-c(51:55)]) # excluding composite channels from  imported ez data
  
##### 2) TRANSFORMATION & SAMPLING #####
# linear, arcsinh, quantile normalization of data + setting up subsamples (need to revisit - do per object type)

  # assign and standardize panels if needed
  panel <- names(master_obj_data[ ,8:48][, -c(4,5)]) # remove metals, composites, other labels
    
  # linear transformation
  data_linear_transform <- master_data
  data_linear_transform[,4:50] <- data_linear_transform[,4:50]*100 # multiply all counts by 100 (linear transform)
  
  # do arcsinh transformation only for the clustering.channels
  data_lin_asinh_transf <- data_linear_transform
  asinh_scale <- 5
  data_lin_asinh_transf[,4:50] <- asinh(data_lin_asinh_transf[,4:50] / asinh_scale)
  
  # PERCENTILE normalize expression values from 0 to 1
  data_normalized <- data_lin_asinh_transf
  normalization_vector <- apply(data_lin_asinh_transf[,4:50], 2, function(x) quantile(x, 0.9999, names = F))
  data_normalized[,4:50] <- t(t(data_normalized[,4:50]) / as.numeric(normalization_vector))
  # check whether you adjusted the range approximately from 0 to 1
  apply(data_normalized[,4:50], 2, max)
  
  # sample (sample_n) for later use by object type
  n_sub_fraction <- 1500

  object_types <- c(object_types, 'cell')
  subsetted_norm_data <- data.frame()
  for (obj_type in object_types) {
    typed_obj <- filter(data_normalized, obj_type_id == obj_type)
    sampled_data <- sample_n(typed_obj, n_sub_fraction)
    subsetted_norm_data <- rbind(subsetted_norm_data, sampled_data)
  }

##### 3) DIAGNOSTIC ANALYSIS #####
# Initial plots to quickly observe the data

##### 4) CLUSTERING & PHENOTYPING #####
  
library(FlowSOM)

# initial cluster number  
num_clusters <- 5
# color palete of choice, here using c25 from brain_data_pkg
palette <- c25

clustering_markers <- c('HistoneH3Lyo', 'MAP2', 'VGAT', 'VGLUT1', 'VGLUT2', 'X8OHGuano', 'GFAP', 'Iba1', 'CD45', 'MCT1', 'CD31')
exclude <- c("C12", "Na23", "Si28", "Ca40", "Background", "Ta181", "Au197", "empty113" )

markers <- setdiff(panel, exclude)

# create a flowFrame from previously imported, transformed data to use in FlowSOM calculation. Currently set to only cluster on 'cell' objects from deepCell data.
flowFrame_cluster <-  new("flowFrame", exprs = as.matrix(subset(master_data, obj_type_id == 'cell', select = panel)))

# run FlowSOM on the above flowFrame which has already had transformation performed
fSOM <- FlowSOM(flowFrame_cluster,
                compensate = F, scale = F, colsToUse = clustering_markers, nClus = num_clusters, seed = 140214)

table(fSOM$metaclustering) 
metaClustering <- fSOM$metaclustering

fSOM_clustering <- data.frame(fSOM$FlowSOM$map$mapping)
colnames(fSOM_clustering) <- c("Cluster", "Value")
fSOM_clustering$Meta <- as.numeric(fSOM$metaclustering[fSOM_clustering[,1]])

fSOMcodes = fSOM$FlowSOM$map$codes
write.csv (fSOMcodes, paste0(outputFolder,'/fSOM_FlowSOM_map_codes.csv'))
fSOMcodesMean = apply(fSOMcodes,2,mean)
write.csv(fSOMcodesMean, paste0(outputFolder,'/mean_of_fSOM_FlowSOM_map_codes.csv'))

# VIZ ABOVE
plot_data = data.frame(table(fSOM_clustering$Meta))
plot_ly(x = plot_data$Var1, y = plot_data$Freq, type = 'bar')
  
##### 5) VISUALISE (DIM REDUC) #####
  
  # tSNE analysis
  # prepare object data - pick a transformed data.frame and revert back to matrix for RtSNE
  data_rtsne <- data_normalized[, panel]
  data_rtsne <- as.matrix(data_rtsne)
  head(data_rtsne)
  colnames(data_rtsne)
  dim(data_rtsne)
  
  # run RtSNE
  set.seed(123)
  out_data_rtsne <- Rtsne(data_rtsne, dims = 2, perplexity = 50, theta = 0.5, #Run Rtnse.multicore if using Ubuntu/S3IT
                        max_iter = 1000, verbose = T, pca = F, check_duplicates=F)
  
  # prepare for plotting (double check assignments are correct between categories and objects)
  data_tsne_plot <- as.data.frame(out_data_rtsne$Y)
  colnames(data_tsne_plot) <- c("tSNE1", "tSNE2")
  data_tsne_plot_all_data <- data_tsne_plot
  data_tsne_plot_all_data <- cbind(data_tsne_plot, data_normalized)
  
  # plot tSNE
  p1 <- ggplot(data_tsne_plot_all_data, aes(x = tSNE1, y = tSNE2, color = obj_type_id)) +
    geom_point(size = 1) + 
    coord_fixed(ratio = 1)
  p1
  
  # UMAP analysis
  data_umap <- data_normalized[, c(panel)]
  data_umap <- as.matrix(data_umap)
  out_data_umap <- umap(data_umap)
  
  data_umap_plot <- as.data.frame(out_data_umap$layout)
  colnames(data_umap_plot) <- c("UMAP1", "UMAP2")
  data_umap_plot_all_data <- data_umap_plot
  data_umap_plot_all_data <- cbind(data_umap_plot, data_normalized)
  
  p1_U <- ggplot(data_umap_plot_all_data, aes(x = UMAP1, y = UMAP2, color = obj_type_id)) +
    geom_point(size = 1) + 
    coord_fixed(ratio = 1)
  p1_U
  
  # UMAP analysis - SUBSETTED DATA
  data_sub_umap <- subsetted_norm_data[, c(panel)]
  data_sub_umap <- as.matrix(data_sub_umap)
  out_data_sub_umap <- umap(data_sub_umap)
  
  data_sub_umap_plot <- as.data.frame(out_data_sub_umap$layout)
  colnames(data_sub_umap_plot) <- c("UMAP1", "UMAP2")
  data_sub_umap_plot_all_data <- data_sub_umap_plot
  data_sub_umap_plot_all_data <- cbind(data_sub_umap_plot, subsetted_norm_data)
  
  p1_sub_U <- ggplot(data_sub_umap_plot_all_data, aes(x = UMAP1, y = UMAP2, color = obj_type_id)) +
    geom_point(size = 1) + 
    coord_fixed(ratio = 1)
    
  p1_sub_U
  
  #EXPRESSION OVERLAYS###################################################### NEEDS UPDATE
  # prepare the expression data
  data.x <- data.frame(data[,c(panel, "cell.id")])
  data.x$cell.id  <- as.factor(data.x$cell.id)
  data.ix.df <- data.frame(data.x[,])
  #library(reshape2)
  data.melt <- melt(data.ix.df, variable.name = "antigen", value.name = "expression")
  joined.expr <- merge(data.melt, data_plot, by = "cell.id")
  #joined.expr <- merge(joined.expr, md, by = "gate.source")
  head(joined.expr)
  
  # plot tSNEs with expression overlayed
  p3 <- ggplot(joined.expr, aes(x = bhSNE1, y = bhSNE2, color = expression)) +
    geom_point(size = 0.25) +
    coord_fixed(ratio = 1) +
    scale_colour_gradient2(low = "black", mid = "yellow", high = "red", midpoint = 0.5, limits = c(0,1), 
                           space = "rgb", guide = "colourbar") +
    facet_wrap(~ antigen, ncol = 4, scales = "free") 
  
  p3 + theme_tsne5
  
  # save the plots as png
  ggsave(filename = "tsne_all.markers.png", plot = p3 + theme_tsne5, 
         scale = 1, width = 10, height = 5, units = c("in"))
##### 6) FURTHER ANALYSIS #####