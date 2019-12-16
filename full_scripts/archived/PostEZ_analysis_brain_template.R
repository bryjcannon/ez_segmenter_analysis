# First Brain [MIBI - ez_segmenter] Paper Data Analysis Script
# Working with post-ez_segmenter data (MATLAB generated) for single object analysis
# Author: Bryan Cannon 2019 (multiple code snippets taken from DM and FH)

##### STEPS #####:

##### 0) LOAD PKGS #####
  
  # load and install necessary packages

  # run pkgs_needed_for_brain_ez
  
  # set seed for downstream analysis
  set.seed(123)

##### 1) LOAD DATA #####
  # move to location of your data, load in fcs and .mat files (spatial information retained here)
  data_folder = '/Volumes/BryJC_Stanford/For_Ez_Segmentor/MedADCase_Hippocampus/'
  setwd(data_folder)
  ezRun = 'ezSegResults_MedAD_MedRes_HC'
  
# create FlowSet from fcs files
  init_fcsNames = list.files(path = paste0(ezRun,'/','fcs_all'), full.names = T, pattern = ".fcs")
  init_flowSet = read.flowSet(files = init_fcsNames, alter.names = T, transformation = FALSE, column.pattern = '($P*R)', invert.pattern = TRUE, emptyValue = FALSE, truncate_max_range = FALSE)
  

##### 2) PRE-PROCESSING #####
# add columns to each file denoting point number, type of tissue, name of objects
  
  # create a list 'point_source' that assigns a point_id to each object
  dim = fsApply(init_flowSet, dim)
  dim = as.numeric(dim[,1])
  
  # set up lists for adding point_ids for each object
  point_names = sampleNames(init_flowSet)
  point_numbers = as.numeric(gsub("[^0-9]*", '', gsub(".*Point[^0-9]*", '', point_names))) # grabs point numbers from filenames and converts to numeric
  point_id = as.vector(x = NULL)
  # set up lists for adding object_type_ids for each object
  object_types = gsub("_dataScaleSize_Point[0-9]*.fcs", "", point_names)
  object_type_id = as.vector(x = NULL)
  # set up lists for adding region_ids for each object
  regions = rep('CA2', length(dim))
  region_id = as.vector(x = NULL)
  
  for(i in 1:length(dim)) { #loop creates the actual vector with point id's for each object
    temp_point_id = rep(point_numbers[i], dim[i])
    point_id = c(point_id, temp_point_id)
    temp_object_type_id = rep(object_types[i], dim[i])
    object_type_id = c(object_type_id, temp_object_type_id)
    temp_region_id = rep(regions[i], dim[i])
    region_id = c(region_id, temp_region_id)
  }
  
  # convert init_data format to matrix (numeric info), then dataFrame (characters + factors), then bind labeling info to dataFrame (keeps string labels)
  objects_data_matrix = fsApply(init_flowSet, Biobase::exprs) #https://support.bioconductor.org/p/109128/ --> explains why use Biobase::exprs
  objects_data_frame = as.data.frame(objects_data_matrix)
  objects_data_frame = cbind(objects_data_frame, point_id, object_type_id, region_id)
  
  # assign and standardize panels if needed
  panel = colnames(objects_data_frame[,7:47][, -c(4,5)])

##### 3) TRANSFORMATION & SAMPLING #####
# linear, arcsinh, quantile normalization of data + setting up subsamples (need to revisit - do per object type)
  
  # sample for later use
  n_sub = 2000
  n = nrow(objects_data_frame)
  set.seed(123)
  subsetted = sample(1:n, n_sub)
  
  # linear transformation
  obj_linear_transform = objects_data_frame
  obj_linear_transform[,3:54] = obj_linear_transform[,3:54]*100 # multiply all counts by 100 (linear transform)
  
  # do arcsinh transformation only for the clustering.channels
  obj_lin_asinh_transf = obj_linear_transform
  asinh_scale = 5
  obj_lin_asinh_transf[,3:54] = asinh(obj_lin_asinh_transf[,3:54] / asinh_scale)
  
  #PERCENTILE normalize expression values from 0 to 1
  obj_normalized = obj_lin_asinh_transf
  normalization_vector = apply(obj_lin_asinh_transf[,3:54], 2, function(x) quantile(x, 0.9999, names = F))
  obj_normalized[,3:54] = t(t(obj_normalized[,3:54]) / as.numeric(normalization_vector))
  
  # check whether you adjusted the range approximately from 0 to 1
  apply(obj_normalized[,3:54], 2, max)

##### 4) DIAGNOSTIC ANALYSIS #####
# Initial plots to quickly observe the data
  
# tSNE analysis
  # prepare object data - pick a transformed data.frame and revert back to matrix for RtSNE
  obj_rtsne = obj_normalized[, c(panel)]
  obj_rtsne = as.matrix(obj_rtsne)
  head(obj_rtsne)
  colnames(obj_rtsne)
  dim(obj_rtsne)
  
  # run RtSNE
  set.seed(123)
  out_obj_rtsne = Rtsne(obj_rtsne, dims = 2, perplexity = 50, theta = 0.5, #Run Rtnse.multicore if using Ubuntu/S3IT
                     max_iter = 1000, verbose = T, pca = F, check_duplicates=F)
  
  # prepare for plotting (double check assignments are correct between categories and objects)
  obj_tsne_plot = as.data.frame(out_obj_rtsne$Y)
  colnames(obj_tsne_plot) = c("tSNE1", "tSNE2")
  obj_tsne_plot_all_data = obj_tsne_plot
  obj_tsne_plot_all_data = cbind(obj_tsne_plot, obj_normalized)
  
  # plot tSNE
  p1 = ggplot(obj_tsne_plot_all_data, aes(x = tSNE1, y = tSNE2, color = object_type_id)) +
    geom_point(size = 1) + 
    coord_fixed(ratio = 1)
  p1 + theme_tsne # from DM / FH script
  
# UMAP analysis
  # use tSNE prepped data for UMAP
  obj_umap = obj_normalized[, c(panel)]
  obj_umap = as.matrix(obj_umap)
  out_obj_umap = umap(obj_umap)
  
  obj_umap_plot = as.data.frame(out_obj_umap$layout)
  colnames(obj_umap_plot) = c("UMAP1", "UMAP2")
  obj_umap_plot_all_data = obj_umap_plot
  obj_umap_plot_all_data = cbind(obj_umap_plot, obj_normalized)
  
  p1_U = ggplot(obj_umap_plot_all_data, aes(x = UMAP1, y = UMAP2, color = object_type_id)) +
    geom_point(size = 1) + 
    coord_fixed(ratio = 1)
  p1_U + theme_tsne

##### 5) CLUSTERING & PHENOTYPING #####
##### 6) VISUALISE (DIM REDUC) #####

##### 7) FURTHER ANALYSIS #####