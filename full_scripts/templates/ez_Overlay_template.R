##### ez_Overlay_current #####
  # Working with post-ez_segmenter data (MATLAB GUI generated) for single object overlay
  # Author: Bryan Cannon 2019 (multiple code snippets taken or built from DM, FH, EFM, DT)

# install and / or load ez_pkgs if you haven't already done so (will also install color scheme)
# to do this make sure you go to ez_lib script and source it. below functions will then work
install_ez_packages(T)
load_ez_packages(T)

##### FILE & DATA SETTINGS - user must adjust these! #####
  # run on data that has been already read and transformed into R using ezAnalysis script (script will output all_data)
  data_for_overlays <- master_obj_data
  # enter head folder where your ez data is stored
  data_folder = '/Volumes/BryJC_Stanford/For_Ez_Segmentor/HiADCase_Hippocampus /denoisedfft_HiResADuci2717J'
  setwd(data_folder)
  # runs, i.e. regions scanned
  ez_runs = c('ezSegResults_CA2+', 'ezSegResults_DG+')
  # shorthand for regions or runs, used in run_type_id
  ez_runs_short = c('CA2', 'DG')
  # where the objects / data / csv / mat files are stored
  ez_data_container = 'objects_points'
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
  for (i in 1:length(ez_runs_short)) {
    region = ez_runs_short[i]
    
    #############NEEDS TO BE CHANGED - WILL MISS POINTS IF NOT IN POURE SEQUENTIAL ORDER
    for (p in point_list) {
      print(paste0("Processing region: ", region, ", Point: ", p))
      
      for (obj_index in 1:length(object_types)) {
        # read matlab matrix of specific object type from specific point into R, if no entry create empty mask
        tryCatch({
          obj_properties <- readMat(paste0(data_folder, '/', ez_runs[i], '/', ez_data_container, '/Point', p, '/', object_types[obj_index], '_objData.mat'))
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
        dir.create(paste0(data_folder, "/", ez_runs[i], "/data_overlays_", ez_runs_short[i]), showWarnings = FALSE)
        dir.create(paste0(data_folder, "/", ez_runs[i], "/data_overlays_", ez_runs_short[i], "/Point", p), showWarnings = FALSE)
        dir.create(paste0(data_folder, "/", ez_runs[i], "/data_overlays_", ez_runs_short[i], "/Point", p, "/TIFs/"), showWarnings = FALSE)
        writeImage(img, paste0(data_folder, "/", ez_runs[i], "/data_overlays_", ez_runs_short[i], "/Point", p, "/TIFs/", object_types[obj_index], ".tif"))
      }
    }
  }
