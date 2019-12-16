library(EBImage)
library(nnet)
library(reshape2)
setwd("/Volumes/KausaliaHD/R4R")

#Run DeepCellAnalysis.R before using this Rscript to get "gated_UMapData" or "all_data" clusters!!!!!!

fcsPath = '/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellSeg_uci2712J/'
fcsRuns = c('190604HiResCA2', '190505HiResDG')
fcsRunsShort = c('CA2', 'DG')
fcsPixel = 'fcs_single_cell_dynamic_expansion'
fcsNames = list.files(path = paste0(fcsPath,fcsRuns[1], '/', fcsPixel), pattern = paste0("^dataScaleSizeFCS.*fcs$"), full.names = T)
allFcsFileNames = fcsNames
allDataFcs = concatFCStoDF(fcsNames, fcsRuns[1], PaperADPanel$Label)

fcsNames = list.files(path = paste0(fcsPath,fcsRuns[2], '/', fcsPixel), pattern = paste0("^dataScaleSizeFCS.*fcs$"), full.names = T)
allFcsFileNames = c(allFcsFileNames, fcsNames)
allDataFcs = rbind(allDataFcs, concatFCStoDF(fcsNames, fcsRuns[2], PaperADPanel$Label))

#palette = brewer.pal(n = numClusters, name = "Paired") # Paired has only 12 colours, choose differetn palette if numClusters > 12

for (i in 1:length(fcsRuns)) {
  region = fcsRuns[i]
  for (p in 1:length(unique(subset(all_data, Region == region)$Point)) ) {
    print(paste0("Processing region: ", region, ". Point: ", p))
    newLmod <- read.csv(paste0("/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellSeg_uci2712J/",region, "/single_cell_dynamic_expansion/Point", p, "/newLmod.csv"), header=FALSE)
    overlayImage = newLmod
    
    # Set background to zero from whatever index it is. Assume: background is the most prevalent object (most pixels) in the image
    bg_finder = melt(newLmod, id.vars=NULL)
    overlayImage[overlayImage == (which.is.max(table(bg_finder$value)))] = 0
    rm(bg_finder)
    
    overlayImage_r = overlayImage
    overlayImage_g = overlayImage
    overlayImage_b = overlayImage
    
    #for each cell in Point p we replace cell ID with colour based on meta number. Change between all_data or gated_UMapData.
    pointdata =  subset(gated_UMapData, Region == region & Point == p)
    for (cell in unique(pointdata$cellLabelInImage)){
      overlayImage_r[overlayImage == cell] = col2rgb(palette)[1, subset(pointdata, cellLabelInImage == cell)$Meta]
      overlayImage_g[overlayImage == cell] = col2rgb(palette)[2, subset(pointdata, cellLabelInImage == cell)$Meta]
      overlayImage_b[overlayImage == cell] = col2rgb(palette)[3, subset(pointdata, cellLabelInImage == cell)$Meta]
    }
    
    overlayImage_r = overlayImage_r / 255
    overlayImage_g = overlayImage_g / 255
    overlayImage_b = overlayImage_b / 255
    
    img <- transpose(rgbImage(Image(as.matrix(overlayImage_r)), Image(as.matrix(overlayImage_g)), Image(as.matrix(overlayImage_b))))
    dir.create(paste0(fcsPath, region, "/expansion_figures_", fcsRunsShort[i], "_Gated"), showWarnings = FALSE)
    dir.create(paste0(fcsPath, region, "/expansion_figures_", fcsRunsShort[i], "_Gated/Point", p), showWarnings = FALSE)
    writeImage(img, paste0(fcsPath, region, "/expansion_figures_", fcsRunsShort[i], "_Gated/Point", p, "/overlay", numClusters, ".tif"))
  }
}

