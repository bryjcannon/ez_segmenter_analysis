# ezSegmenter data prep for DeepCell Analysis
library("R.matlab")

ez_folder <- '/Volumes/BryJC_Stanford/For_Ez_Segmentor/HiADCase_Hippocampus /denoisedfft_HiResADuci2717J'
setwd(ez_folder)
ezRun <- c('ezSegResults_CA2+', 'ezSegResults_DG+')

for (run in ezRun) {
  path = paste0(run,'/',"objects_points")
  
  pointList <- list.files(path)
  pointList <- gsub("Point", "", pointList)
  
  for (point in pointList) {
    pointFiles <- list.files(paste0(path, "/Point", point))
    
    pointFiles <- pointFiles[grepl("_cellData", pointFiles, ignore.case = T)]
    
    for (cellDataMat in pointFiles) {
      matlabData <- readMat(paste0(path, "/Point", point, "/", cellDataMat))
      newLmod <- matrix(unname(matlabData["mapped.obj.ids"]))[[1]]
      
      cellDataMat <- gsub("_cellData.mat", "", cellDataMat, ignore.case = T)
      write.csv(newLmod, paste0(path, "/Point", point, "/", cellDataMat, "_newLmod.csv"))
    }
  }
}
