# First Brain [MIBI - ez_segmenter] Paper Data Analysis Script
# Working with post-ez_segmenter data (MATLAB generated) for single object analysis
# Author: Bryan Cannon 2019 (multiple code snippets taken or built from from DM, FH, EFM, DT)

##### STEPS #####:

##### 0) LOAD PKGS #####
  
  # load necessary packages


  # set seed for downstream analysis
  set.seed(123)

##### 1a) LOAD ezSeg DATA #####
  # move to location of your data, load in csv and .mat files (spatial information stored in .mat)
  ez_folder <- '/Volumes/BryJC_Stanford/For_Ez_Segmentor/HiADCase_Hippocampus /denoisedfft_HiResADuci2717J'
  setwd(ez_folder)
  ezRun <- c('ezSegResults_CA2+', 'ezSegResults_DG+')
  object_types <- c('amyloidopathy', 'tauopathy', 'microglia_process', 'vessel_CD31_CD105', 'vessel_MCT1')
  
# create data.frame containing single object data from csv files - for each run and each object type
  master_obj_data <- data.frame()
  
  for (run in ezRun) {
    obj_data_raw_all <- data.frame()
    
    for (obj_type in object_types) {
      csv_names <- list.files(path = paste0(run,'/','objects_points'), recursive = T, full.names = T, pattern = paste0(obj_type, "_dataScaleSize.csv")) # read in csv files for object type
      csv_names <- mixedsort(csv_names)
      
      obj_data_raw <- lapply(csv_names, read.csv) %>% bind_rows() # grab data from csv's then convert data to data.frame
      
      obj_type_id <- rep(obj_type, dim(obj_data_raw)[1]) # create column of length object number with obj_type info
      obj_data_raw <- cbind(obj_data_raw, obj_type_id) # add obj_type_id to data
      
      obj_data_raw_all <- rbind(obj_data_raw_all, obj_data_raw) # collate to run data.frame
    }
    
    run_type_id <- rep(run, dim(obj_data_raw_all)[1]) # create column of length object number with run_type info
    obj_data_raw_all <- cbind(obj_data_raw_all, run_type_id) # add run_type_id to data
    
    master_obj_data <- rbind(master_obj_data, obj_data_raw_all) # collate to master data.frame
    
    rm(obj_data_raw)
    rm(obj_data_raw_all)
  }

##### 1b) LOAD Deep Cell DATA #####
  # move to location of your data, load in fcs and .mat files (spatial information retained here)
  data_folder_DC <- '/Volumes/BryJC_Stanford/'
  setwd(data_folder_DC)
  deepcellRun <- c('190505HiResDG', '190604HiResCA2')
  
  # create FlowSet from fcs files
  
  master_cell_data <- data.frame()
  
  for (run in deepcellRun) {
    fcs_names <- c(fcs_names, list.files(path = paste0(run,'/','fcs_single_cell_dynamic_expansion'), full.names = T, pattern = "(dataScaleSizeFCS.*fcs)"))
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
    
    run_type_id <- rep(run, dim(cell_data_raw)[1]) # create column of length object number with run_type info
    cell_data_raw <- cbind(cell_data_raw, run_type_id) # add run_type_id to data
    
    master_cell_data <- rbind(master_cell_data, cell_data_raw)
    
    rm(cell_data_matrix)
    rm(cell_flowSet)
    rm(cell_data_raw)
  }

##### 1b) Combine Deep Cell and ezSeg Data #####
  
  master_cell_data <- plyr::rename(master_cell_data, c('cellLabelInImage' = 'obj_id', 'cellSize' = 'obj_size', 'cell_type_id' = 'obj_type_id')) # rename deep cell clolumns to match ez names
  master_cell_data <- master_cell_data %>% select(point_id, everything()) # move point_id to start of deep cell data to match ez column order
  
  master_data <- rbind(master_cell_data, master_obj_data[,-c(51:55)])
      
##### 2) TRANSFORMATION & SAMPLING #####
# linear, arcsinh, quantile normalization of data + setting up subsamples (need to revisit - do per object type)

  # assign and standardize panels if needed
  panel <- names(master_obj_data[ ,8:48][, -c(4,5)]) # remove metals, composites, other labels
    
  # linear transformation
  obj_linear_transform <- master_obj_data
  obj_linear_transform[,4:55] <- obj_linear_transform[,4:55]*100 # multiply all counts by 100 (linear transform)
  
  # do arcsinh transformation only for the clustering.channels
  obj_lin_asinh_transf <- obj_linear_transform
  asinh_scale <- 5
  obj_lin_asinh_transf[,4:55] <- asinh(obj_lin_asinh_transf[,4:55] / asinh_scale)
  
  # PERCENTILE normalize expression values from 0 to 1
  obj_normalized <- obj_lin_asinh_transf
  normalization_vector <- apply(obj_lin_asinh_transf[,4:55], 2, function(x) quantile(x, 0.9999, names = F))
  obj_normalized[,4:55] <- t(t(obj_normalized[,4:55]) / as.numeric(normalization_vector))
  # check whether you adjusted the range approximately from 0 to 1
  apply(obj_normalized[,4:55], 2, max)
  
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
  
  p1_U <- ggplot(obj_umap_plot_all_data, aes(x = UMAP1, y = UMAP2, color = obj_type_id)) +
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
  
  p1_sub_U <- ggplot(obj_sub_umap_plot_all_data, aes(x = UMAP1, y = UMAP2, color = obj_type_id)) +
    geom_point(size = 1) + 
    coord_fixed(ratio = 1)
    
  p1_sub_U

##### 6) FURTHER ANALYSIS #####