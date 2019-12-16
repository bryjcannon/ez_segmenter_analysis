# First Brain [MIBI - ez_segmenter] Paper Data Analysis Script
# Working with post-ez_segmenter data (MATLAB generated) for single object analysis
# Author: Bryan Cannon 2019 (multiple code snippets taken or built from from DM, FH, EFM, DT)

##### STEPS #####:

##### 0) LOAD PKGS #####
  
  # load necessary packages
  # go to packages file, this line will be converted into a function later

  # set seed for downstream analysis
  seed <- 123 
  set.seed(seed)

##### 1) LOAD ezSeg DATA #####
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
  
##### 2) TRANSFORMATION & SAMPLING #####
# linear, arcsinh, quantile normalization of data + setting up subsamples (need to revisit - do per object type)
  
  # assign and standardize panels if needed
  panel_start = 8
  panel_end = 48
  panel <- names(master_obj_data[,panel_start:panel_end][, -c(4,5)]) # remove metals, composites, other labels
  
  # linear transformation
  obj_linear_transform <- master_obj_data
  obj_linear_transform[,panel_start:panel_end] <- obj_linear_transform[,panel_start:panel_end]*100 # multiply all counts by 100 (linear transform)
  
  # do arcsinh transformation only for the clustering.channels
  obj_lin_asinh_transf <- obj_linear_transform
  asinh_scale <- 5
  obj_lin_asinh_transf[,panel_start:panel_end] <- asinh(obj_lin_asinh_transf[,panel_start:panel_end] / asinh_scale)
  
  # PERCENTILE normalize expression values from 0 to 1
  obj_normalized <- obj_lin_asinh_transf
  normalization_vector <- apply(obj_lin_asinh_transf[,panel_start:panel_end], 2, function(x) quantile(x, 0.9999, names = F))
  obj_normalized[,panel_start:panel_end] <- t(t(obj_normalized[,panel_start:panel_end]) / as.numeric(normalization_vector))
  # check whether you adjusted the range approximately from 0 to 1
  apply(obj_normalized[,panel_start:panel_end], 2, max)
  
  # sample (sample_n) for later use by object type
  n_sub_fraction <- 1500

  subsetted_norm_data <- data.frame()
  for (obj_type in object_types) {
    typed_obj <- filter(obj_normalized, obj_type_id == obj_type)
    sampled_obj <- sample_n(typed_obj, n_sub_fraction)
    subsetted_norm_data <- rbind(subsetted_norm_data, sampled_obj)
  }

##### 3) DIAGNOSTIC ANALYSIS #####
# Initial plots to quickly observe the data

##### 4) CLUSTERING & PHENOTYPING #####
##### 5) VISUALISE (DIM REDUC) #####
  
  # tSNE analysis
  # prepare object data - pick a transformed data.frame and revert back to matrix for RtSNE
  obj_rtsne <- obj_normalized[, panel]
  obj_rtsne <- as.matrix(obj_rtsne)
  head(obj_rtsne)
  colnames(obj_rtsne)
  dim(obj_rtsne)
  
  # run RtSNE
  set.seed(123)
  out_obj_rtsne <- Rtsne(obj_rtsne, dims = 2, perplexity = 50, theta = 0.5, #Run Rtnse.multicore if using Ubuntu/S3IT
                        max_iter = 1000, verbose = T, pca = F, check_duplicates=F)
  
  # prepare for plotting (double check assignments are correct between categories and objects)
  obj_tsne_plot <- as.data.frame(out_obj_rtsne$Y)
  colnames(obj_tsne_plot) <- c("tSNE1", "tSNE2")
  obj_tsne_plot_all_data <- obj_tsne_plot
  obj_tsne_plot_all_data <- cbind(obj_tsne_plot, obj_normalized)
  
  # plot tSNE
  p1 <- ggplot(obj_tsne_plot_all_data, aes(x = tSNE1, y = tSNE2, color = obj_type_id)) +
    geom_point(size = 1) + 
    coord_fixed(ratio = 1)
  p1
  
  # UMAP analysis
  obj_umap <- obj_normalized[, c(panel)]
  obj_umap <- as.matrix(obj_umap)
  out_obj_umap <- umap(obj_umap)
  
  obj_umap_plot <- as.data.frame(out_obj_umap$layout)
  colnames(obj_umap_plot) <- c("UMAP1", "UMAP2")
  obj_umap_plot_all_data <- obj_umap_plot
  obj_umap_plot_all_data <- cbind(obj_umap_plot, obj_normalized)
  
  p1_U <- ggplot(obj_umap_plot_all_data, aes(x = UMAP1, y = UMAP2, color = ApoE4)) +
    geom_point(size = 1) + 
    coord_fixed(ratio = 1)
  p1_U
  
  # UMAP analysis - SUBSETTED DATA
  obj_sub_umap <- subsetted_norm_data[, c(panel)]
  obj_sub_umap <- as.matrix(obj_sub_umap)
  out_obj_sub_umap <- umap(obj_sub_umap)
  
  obj_sub_umap_plot <- as.data.frame(out_obj_sub_umap$layout)
  colnames(obj_sub_umap_plot) <- c("UMAP1", "UMAP2")
  obj_sub_umap_plot_all_data <- obj_sub_umap_plot
  obj_sub_umap_plot_all_data <- cbind(obj_sub_umap_plot, subsetted_norm_data)
  
  p1_sub_U <- ggplot(obj_sub_umap_plot_all_data, aes(x = UMAP1, y = UMAP2, color = CD47)) +
    geom_point(size = 1) + 
    coord_fixed(ratio = 1)
    
  p1_sub_U

##### 6) FURTHER ANALYSIS #####