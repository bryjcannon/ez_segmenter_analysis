##### ez_Overlay_working #####
# Working with post-ez_segmenter data (MATLAB GUI generated) for single object overlay, use this script only for testing new code #
# Author: Bryan Cannon 2019 (multiple code snippets taken or built from DM, FH, EFM, DT)

# install and / or load ez_pkgs if you haven't already done so (will also install color scheme)
install_ez_packages(T)
load_ez_packages(T)

##### FILE & DATA SETTINGS - user must adjust these! #####
# run on data that has been already read and transformed into R using ezAnalysis script (script will output all_data)
data_for_overlays <- master_obj_data
# enter head folder where your ez data is stored
dataPath = '/Volumes/BryJC_Stanford/For_Ez_Segmentor/HiADCase_Hippocampus /denoisedfft_HiResADuci2717J'
setwd(dataPath)
# runs, i.e. regions scanned
dataRuns = c('ezSegResults_CA2+', 'ezSegResults_DG+')
# shorthand for regions or runs, used in run_type_id
dataRunsShort = c('CA2', 'DG')
# where the objects / data / csv / mat files are stored
dataContainer = 'objects_points'
# object types you plan to cluster)
object_types <- c('amyloidopathy', 'tauopathy', 'microglia_process', 'vessel_CD31_CD105', 'vessel_MCT1', 'cell')
# boolean denoting if clustering was performed on the objects (T) or not (F)
objects_clustered = F
# color palete of choice
palette <- c25
# the resolution of your images, e.g. 512 x 512, 1024 x 1024, etc.
resolution_dim = c(1024, 1024)
# list of actual point numbers (i.e. which points you actually want to overlay, e.g. c(4,6:9) for Points 4 and 6 through 9)
point_list <- c(1:3)

##### OVERLAY CONSTRCTION PROCESS - run entire code block #####
# for each run, for each point, for each object, create and save colored overlays
for (i in 1:length(dataRunsShort)) {
  region = dataRunsShort[i]
  
  #############NEEDS TO BE CHANGED - WILL MISS POINTS IF NOT IN POURE SEQUENTIAL ORDER
  for (p in point_list) {
    print(paste0("Processing region: ", region, ", Point: ", p))
    
    for (obj_index in 1:length(object_types)) {
      # read matlab matrix of specific object type from specific point into R, if no entry create empty mask
      tryCatch({
        obj_properties <- readMat(paste0(dataPath, '/', dataRuns[i], '/', dataContainer, '/Point', p, '/', object_types[obj_index], '_objData.mat'))
        # pull out mapped object ids as a matrix (defaults to image dimensions)
        obj_mask <- obj_properties$mapped.obj.ids
      }, error = function(err) {
        obj_mask <- zeros(resolution_dim[1], resolution_dim[2])
      }) # end of tryCatch loop
      
      # create three channels (red, green, blue) for the mask to use in later Image creation - inits with empty matrix size of the image
      obj_mask_r = zeros(resolution_dim[1], resolution_dim[2])
      obj_mask_g = zeros(resolution_dim[1], resolution_dim[2])
      obj_mask_b = zeros(resolution_dim[1], resolution_dim[2])
      
      #for each object in point p we replace object ID with colour based on meta number.
      pointdata =  subset(data_for_overlays, run_type_id == region & point_id == p & obj_type_id == object_types[obj_index])
      for (uniq_object in unique(pointdata$obj_id)){
        # try to set rgb colors according to palette - object id match, otherwise set it all to 0
        tryCatch({
          if (objects_clustered == F) {
            obj_mask_r[obj_mask == uniq_object] = col2rgb(palette)[1, obj_index]
            obj_mask_g[obj_mask == uniq_object] = col2rgb(palette)[2, obj_index]
            obj_mask_b[obj_mask == uniq_object] = col2rgb(palette)[3, obj_index]
          }
          else if (objects_clustered == T) {
            obj_mask_r[obj_mask == uniq_object] = col2rgb(palette)[1, subset(pointdata, obj_id == uniq_object)$cluster_id]
            obj_mask_g[obj_mask == uniq_object] = col2rgb(palette)[2, subset(pointdata, obj_id == uniq_object)$cluster_id]
            obj_mask_b[obj_mask == uniq_object] = col2rgb(palette)[3, subset(pointdata, obj_id == uniq_object)$cluster_id]
          }
        }, error = function(err) {
          obj_mask_r[obj_mask == uniq_object] = 0
          obj_mask_g[obj_mask == uniq_object] = 0
          obj_mask_b[obj_mask == uniq_object] = 0
        })
      }
      obj_mask_r = obj_mask_r / 255
      obj_mask_g = obj_mask_g / 255
      obj_mask_b = obj_mask_b / 255
      
      img <- transpose(rgbImage(Image(obj_mask_r), Image(obj_mask_g), Image(obj_mask_b)))
      dir.create(paste0(dataPath, "/", dataRuns[i], "/data_overlays_", dataRunsShort[i]), showWarnings = FALSE)
      dir.create(paste0(dataPath, "/", dataRuns[i], "/data_overlays_", dataRunsShort[i], "/Point", p), showWarnings = FALSE)
      dir.create(paste0(dataPath, "/", dataRuns[i], "/data_overlays_", dataRunsShort[i], "/Point", p, "/TIFs/"), showWarnings = FALSE)
      writeImage(img, paste0(dataPath, "/", dataRuns[i], "/data_overlays_", dataRunsShort[i], "/Point", p, "/TIFs/", object_types[obj_index], ".tif"))
      
      # composite image = new vector of images
      # for each object type in object types
      # import the newModL csv or matrix from .mat file for this point, this object
      # if clusters in the data: TBD----------------------------------------------------> still need to think about this, more relevant for DM
      # if no clusters in the data, then assign a color to this object type
      # convert newModL csv (now a matrix or data frame) into an RGB image using color type
      # hold image in composite queue.
      # combo imgae = use magick to merge these images into one
      # save image in the point for later tiling in matlab (eventually R)
      
      # to incorporate deep cell data, use code from deep cell overlay to produce 
      # cluster annotation to RGB assignment and produce an image (ensure no clash with ez objects)
      # add image to composite image and proceed with merge as planned
    }
  }
}
