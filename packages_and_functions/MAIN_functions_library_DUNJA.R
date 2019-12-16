###############################################
#collections of functions by Felix J. Hartmann#
###############################################

#load and if nececcary install libraries
source("http://bioconductor.org/biocLite.R")
if (!require(flowCore)) {biocLite("flowCore")} 
if (!require(flowVS)) {biocLite("flowVS")} 
if (!require(flowViz)) {biocLite("flowViz")} 
if (!require(RColorBrewer)) {biocLite("RColorBrewer")} 
if (!require(gplots)) {biocLite("gplots")} 
if (!require(ggplot2)) {biocLite("ggplot2")} 
if (!require(samr)) {biocLite("samr")} 
if (!require(lattice)) {biocLite("lattice")} 
if (!require(flowStats)) {biocLite("flowStats")} 
if (!require(gdata)) {biocLite("gdata")}
if (!require(Rtsne)) {biocLite("Rtsne")}
if (!require(FlowSOM)) {biocLite("FlowSOM")}
if (!require(plyr)) {biocLite("plyr")}
if (!require(pryr)) {biocLite("pryr")}
if (!require(doBy)) {biocLite("doBy")}
if (!require(scales)) {biocLite("scales")}
if (!require(mixOmics)) {biocLite("mixOmics")}
if (!require(reshape2)) {biocLite("reshape2")}
if (!require(plotly)) {biocLite("plotly")}
if (!require(Rmisc)) {biocLite("Rmisc")}
if (!require(Hmisc)) {biocLite("Hmisc")}
#if (!require(Rtsne.multicore)) {biocLite("Rtsne.multicore")} #don't run on a Mac


`%notin%` <- function(x,y) !(x %in% y) 

# sem function
b.median <- function(data, num) {
  resamples <- lapply(1:num, function(i) sample(na.omit(data), replace = T))
  r.median <- sapply(resamples, median)
  std.err <- sqrt(var(r.median))
  std.err
  #list(std.err=std.err, resamples=resamples, medians=r.median)   
}


# median function
median_fun  <- function(data, median.column, minimum.cells){
  
  #make empty results matrix
  max <- max(data[,median.column])
  min <- min(data[,median.column])
  length <- max-min+1
  results <- matrix(data = NA, nrow = length, ncol = ncol(data))
  colnames(results) <- colnames(data)
  rownames(results) <- rownames(md)
  
  #loop through all gate_sources and calculate median for all samples (with minimum cells)
  for(i in 1:length) {
    temp <- data[data[,median.column] == (i+min-1),]
    med_temp <- apply(temp, 2, function(x){median(x, na.rm = TRUE)})
      results[i,] <- as.vector(med_temp)
  }
  results
}



# bin function
bin_fun  <- function(data, median.column, minimum.cells){
  
  #make empty results matrix
  cluster_min <- min(data[,median.column])
  cluster_max <- max(data[,median.column])
  cluster_num <- (cluster_min-cluster_max-1)*(-1)

  #move it to positive space
  data[,median.column] <- data[,median.column]+(cluster_min*(-1))+1
  
  results <- matrix(data = NA, nrow = cluster_num, ncol = ncol(data))
  colnames(results) <- colnames(data)
  
  #loop through all gate_sources and calculate median for all samples (with minimum cells)
  for(i in 1:cluster_num) {
    temp <- data[data[,median.column] == i,]
    if(nrow(temp) >= minimum.cells) {
      med_temp <- apply(temp, 2, function(x) {median(x, na.rm = T)})
      results[i,] <- as.vector(med_temp)
    } else {  }
  }
  results[,median.column]  <-  (results[,median.column]+((cluster_min)-1))
  results
}



#gating function
cutoff_fun  <- function(ff, cutoff, jit_num){
  exp_ff <- exprs(ff)
  red_mat <- exp_ff[,names(cutoff)]
  comp_mat <- cutoff < t(red_mat)
  comp_mat[comp_mat == F] <- 0
  comp_mat[comp_mat == T] <- 1
  comp_mat <- (t(comp_mat))
  if (jit_num!=0) {
    comp_mat <- jitter(comp_mat, jit_num)
  } 
  colnames(comp_mat) <- paste("cut", colnames(comp_mat), sep = "_")
  return_ff  <- cbind2(ff,comp_mat)
  mergedesc <- c(colnames(ff),colnames(comp_mat))
  return_ff@parameters$desc <- mergedesc
  return_ff
}


#concatenation function
concat_fun  <- function(fs_loops){
  exp_loops <- fsApply(fs_loops, function(x) exprs(x))
  concat_ff <- fs_loops[[2]]
  exprs(concat_ff)  <- exp_loops
  concat_ff
}



#median frequency function
freq_fun  <- function(ff,fm,cc,cluster_num){
  vec  <- exprs(ff)[,cc]
  f_if <- vector(mode="integer", length=0)
  f_else <- as.integer(rep(NA, cluster_num))
  
  if(length(vec)>mn)
  {
    for(i in 1:cluster_num)
    {
      f_if[i]  <- 100*(length(which(vec==i))/length(vec))
    }
    fm  <- rbind(fm, f_if)
  } else {
    fm  <- rbind(fm, f_else)
  }
  fm #returns fm
}


#subsampling function
sample_fun  <- function(ff, sub_vec){
  mat  <- exprs(ff)
  mat_1 <- mat[sample(nrow(mat), replace=T, size=sub_vec), ]
  colnames(mat_1) <- as.vector(colnames(mat))
  exprs(ff)  <- mat_1
  ff
}


#split merged files function
split_fun  <- function(exp_ff, names){
  #exp_ff <- exprs(ff)
  for (i in min(exp_ff[,"gate_source"]):max(exp_ff[,"gate_source"])) {
    temp_fcs <- as.matrix(exp_ff[exp_ff[,"gate_source"]==i,])
    ff_split  <- flowFrame(exprs = temp_fcs)
    if (i==min(exp_ff[,"gate_source"])) {
      fs_new <- flowSet(ff_split)
    } else {
      fs_new <- rbind2(fs_new, ff_split)
    }
    sampleNames(fs_new)  <- names[1:length(fs_new)]
  }
  fs_new
}

# frequency function for combined dataset
freq_fun_mat  <- function(data, source.column, cluster.column, minimum.cells, names.vec){
  
  #make empty results matrix 
  gate_max  <- max(data[,source.column])
  gate_min <- min(data[,source.column])
  gate_length <- gate_max-gate_min+1
  cluster_num <- max(data[,cluster.column])
  results <- matrix(data = NA, nrow = gate_length, ncol = cluster_num)
  colnames(results) <- paste("cluster_", 1:cluster_num, sep="")
 # rownames(results) <- names.vec
  
  
  for(ib in 1:gate_length) {
    
    vec <- data[data[,source.column] == (ib+gate_min-1), cluster.column]
    f_if <- vector(mode="integer", length=0)
    f_else <- as.integer(rep(NA, cluster_num))
    
    if(length(vec)>minimum.cells) {
      
      for(i in 1:cluster_num) {
        f_if[i]  <- 100*(length(which(vec==i))/length(vec)) #put total cell number instead of length(vec)?
      }
      results[ib,] <- f_if 
    } else {}
    
  }
  results
}




# percentile function for  a whole flowset
percentile_fun  <- function(ff, quantiles_vec){
  exp_ff  <- exprs(ff)[,names(quantiles_vec)]
  sc <- t(t(exp_ff) / as.numeric(quantiles_vec))
  exprs(ff)[,names(quantiles_vec)]  <- as.matrix(sc)
  ff
}


# this is for ff specific scaling
ff_percentile_fun  <- function(ff, lc, low_ms, q_vec){
  #for low
  lc_low <- lc[lc %in% low_ms]
  exp_ff_low  <- exprs(ff)[,lc_low]
  quantiles_fcs_low <- apply(exp_ff_low, 2, function(x) quantile(x, max(q_vec), names = FALSE))
  sc_low <- t(t(exp_ff_low) / as.numeric(quantiles_fcs_low))
  exprs(ff)[,lc_low]  <- as.matrix(sc_low)
  #for high
  lc_high <- lc[!lc %in% low_ms]
  exp_ff_high  <- exprs(ff)[,lc_high]
  quantiles_fcs_high <- apply(exp_ff_high, 2, function(x) quantile(x, min(q_vec), names = FALSE))
  sc_high <- t(t(exp_ff_high) / as.numeric(quantiles_fcs_high))
  exprs(ff)[,lc_high]  <- as.matrix(sc_high)
  #return ff
  ff
}



heat_fun  <- function(ce, cc, clustering, minimum_cells, color) {
  #initialize matrix
  #rm(return_mat)
  mn <- minimum_cells
  my_color <- color
  return_mat <- matrix(,nrow=max(ce[,cc]), ncol=ncol(ce))
  colnames(return_mat) <- colnames(ce)
  
  #go through every cluster and calculate median for all channels
  cluster_num <- max(ce[,cc])
  for(i in 1:cluster_num)
  {
    temp_mat <- ce[which(ce[,cc]==i),]
    return_mat[i,] <- as.matrix(apply(temp_mat, 2, median))
  }
  return_mat
  
  #plot heatmap
  heatmap.2(return_mat[,colnames(return_mat) %in% clustering], 
            scale="none",
            Colv=T,
            dendrogram="both",
            trace="none",
            #hclustfun=method,
            col=my_color)
}

##### colour palettes
# display.brewer.pal(8, "Dark2")

#For palette choices:
library("RColorBrewer")
display.brewer.all()

# colours for clusters
db1 <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))
db2 <- c((dichromat_pal("DarkRedtoBlue.12")(12)))
db3 <- c(brewer.pal(8, "Dark2"))
db4 <- c(brewer.pal(12, "Paired"))
db5 <- c(dichromat_pal("BluetoDarkOrange.18")(18), "black", "grey")
db6 <- c(dichromat_pal("DarkRedtoBlue.18")(18), "black", "grey")
db7 <- c(db1, db6)
db8 <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))
db9 <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Set3"), brewer.pal(8, "Set1"))
db10 <- c(brewer.pal(8, "Dark2"), brewer.pal(7, "Set3"), "grey")
db11 <- c ("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
           "#1F78B4", "#666666", "#FB9A99", "blue")

# Custom colors
# dbr <- sample(db9, 28, replace = F)

db19 <- c("#7570B3", "#D95F02", "#D95F02", "#E7298A", "#66A61E" , "#E6AB02",
          "#A6761D", "#666666", "#A6CEE3", "#A6CEE3", "indianred", "gold", 
          "lightseagreen", "brown" , "hotpink", "#7570B3", "grey" , "firebrick1", "lightsalmon")

db13 <- c("#1B9E77","gold" , "#7570B3", "#E7298A", "#D95F02" ,  "lightpink1", "#666666",
          "#57B0FF", "#66CD00", "#1F78B4", "#B2DF8A", "firebrick1", 
          "seagreen1", "brown" , "#6A3D9A", "mediumorchid", "grey" , "gold", "gold",
          "lightpink1", "#BC80BD", "gold", "#80B1D3", "gold","#80B1D3", "mediumorchid")

db13.3 <- c("gold", "#7570B3", "#E7298A", "hotpink4", "hotpink", 
            "#D95F02", "#D95F02", "#D95F02", "#D95F02", "#D95F02", "lightpink1",
            "#666666", "#666666", "#666666", "#666666",
            "#57B0FF", "dodgerblue2", "darkturquoise", "lightblue", 
            "#66CD00", "mediumorchid1")

db13.4 <- c("#1B9E77", "gold", "#7570B3", "#E7298A", "hotpink4", "hotpink", 
            "#D95F02", "lightpink1",
            "#666666",
            "#57B0FF", "dodgerblue2", "darkturquoise", "lightblue", 
            "#66CD00", "mediumorchid1", "#1B9E77")

db13.5 <- c("gold", "#7570B3", "#E7298A", "hotpink4", "hotpink", 
            "#D95F02", "#D95F02", "#D95F02", "#D95F02", "#D95F02", "lightpink1",
            "#666666", "#666666", "#666666", "#666666",
            "#57B0FF", "dodgerblue2", "darkturquoise", "lightblue", 
            "#66CD00", "mediumorchid1", "#1B9E77")

db14 <- c("#1B9E77", "gold", "#7570B3", "#DCDCDC")

db15 <- c("gold" , "#7570B3", "#E7298A", "#D95F02" ,  "lightpink1", "#666666",
          "#57B0FF", "#66CD00", "#1F78B4", "#B2DF8A", "firebrick1", 
          "seagreen1", "brown" , "#6A3D9A", "mediumorchid", "grey" , "gold", "gold",
          "lightpink1", "#BC80BD", "gold", "#80B1D3", "gold","#80B1D3", "mediumorchid")


e#Theme for tsne  
theme_tsne <-  theme (panel.spacing = unit(1.3, "lines"), 
                       strip.text = element_text(size = rel(1.5)), 
                       axis.ticks.length = unit(0.3, "lines"),
                       axis.text = element_blank(),
                       axis.title.x = element_text(size = 35),
                       axis.title.y = element_text(size = 35),
                       plot.background = element_rect(color="white"),
                       strip.background = element_blank(), 
                       panel.background = element_blank(),
                       axis.ticks = element_blank(),
                       aspect.ratio = 1,
                       legend.text=element_text(size=15), 
                       legend.title=element_text(size=15))

theme_tsne1 <-  theme (panel.spacing = unit(1.3, "lines"), 
                      strip.text = element_text(size = rel(1.5)), 
                      axis.ticks.length = unit(0.3, "lines"),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      plot.background = element_rect(color="white"),
                      strip.background = element_blank(), 
                      panel.background = element_blank(),
                      aspect.ratio = 1,
                      legend.text=element_text(size=15), 
                      legend.title=element_text(size=15))

theme_tsne2 <-  theme (panel.spacing = unit(1.3, "lines"), 
                       strip.text = element_text(size = rel(1.5)), 
                       axis.ticks.length = unit(0.3, "lines"),
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       plot.background = element_rect(color="white"),
                       strip.background = element_blank(), 
                       panel.background = element_blank(),
                       aspect.ratio = 1,
                       legend.text=element_text(size=15), 
                       legend.title=element_blank())

theme_tsne3 <-  theme (panel.spacing = unit(1.3, "lines"), 
                       strip.text = element_text(size = rel(1.4)), 
                       axis.ticks.length = unit(0.3, "lines"),
                       axis.text = element_blank(),
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       plot.background = element_blank(),
                       strip.background = element_blank(), 
                       panel.background = element_rect(fill = "black"),
                       panel.grid = element_blank(),
                       axis.ticks = element_blank(),
                       aspect.ratio = 1,
                       legend.text=element_text(size=25), 
                       legend.title=element_text(size=25),
                       legend.key = element_rect(fill = "black"))

theme_tsne4 <- theme(axis.text.x = element_blank(),
                    axis.text.y = element_blank(),  
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.ticks = element_blank(),
                    legend.position = "none",
                    aspect.ratio =1,
                    panel.spacing = unit(1.3, "lines"), 
                    strip.text = element_text(size = rel(1.4)))

theme_tsne5 <-  theme (panel.spacing = unit(1.3, "lines"), 
                       strip.text = element_text(size = rel(1.5)), 
                       axis.ticks.length = unit(0.3, "lines"),
                       axis.text = element_blank(),
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       plot.background = element_rect(color="white"),
                       strip.background = element_blank(), 
                       panel.background = element_blank(),
                       axis.ticks = element_blank(),
                       aspect.ratio = 1,
                       legend.position = "none")
                       
theme_biaxial <- theme(axis.title.x = element_text(size=20),
                       axis.title.y = element_text(size=20),
                       axis.text = element_text(size=20))

theme_biaxial2 <- theme_classic()  + 
  theme(
    axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
    axis.text.y = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=0,face="plain"),
    axis.ticks.length = unit(1, "lines"),
    axis.title.x = element_text(colour="black",size=30,angle=0,hjust=.5,vjust=0,face="plain"),
    axis.title.y = element_text(colour="black",size=30,angle=90,hjust=.5,vjust=.5,face="plain"),
    legend.position = "none")


theme_pub <-  theme(strip.text = element_text(size = rel(1.25)), 
                    strip.background = element_blank(), 
                    axis.text = element_text(size = rel(1.25), color = "black"),
                    axis.text.x = element_text(angle = 90, hjust = 1, color = "black", size = rel(1)),
                    panel.spacing = unit(1.3, "lines"), 
                    axis.ticks.length = unit(0.3, "lines"),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.line = element_blank(),
                    plot.background = element_rect(color="white"),
                    panel.background = element_rect(fill="white", color=NA))

theme_pub2 <-  theme(strip.text = element_text(size = rel(1.25)), 
                    strip.background = element_blank(), 
                    axis.text = element_text(size = rel(1.25), color = "black"),
                    axis.text.x = element_text(angle = 0, hjust = 1, color = "black", size = rel(0.5)),
                    panel.spacing = unit(1.3, "lines"), 
                    axis.ticks.length = unit(0.3, "lines"),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.line = element_blank(),
                    plot.background = element_rect(color="white"),
                    panel.background = element_rect(fill="white", color=NA))

theme_legend <- theme(legend.title=element_blank(), legend.text=element_text(size=15))

theme_facs <- theme(strip.background = element_blank(),
                    #aspect.ratio = 1,
                    panel.background = element_rect(colour = "black", size = 1, fill = "white"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    #axis.text = element_text(size = rel(2), color = "black"),
                    #axis.title = element_text(size = rel(3.3)),
                    legend.position = "right",
                    #axis.ticks.length = unit(0.5, "lines"),
                    #title = element_text(size = rel(2)),
                    #plot.title = element_text(size = rel(2)),
                    #strip.text = element_text(size = rel(2))
)

##### metaclustering/cluster combination function #####


meta_fun <- function(ff, cc) {
  test <- exprs(ff)[,cc]
  manual_metaclusters <- as.vector(rep(0, length(test)))
  test2 <- cbind(test, manual_metaclusters)
  #head(test2)
  
  #defines the clusters
  cd4_naive <- c(4,2,13,10,18)
  cd4_memory <- c(1,12,7)
  cd8_naive <- c(3,17)
  cd8_memory  <- c(5,14)
  gd_cells <- c(15)
  nk_cells <- c(6)
  mono <- c(16)
  b_cells <- c(9,8,11)
  
  #give them the right metacluster number
  test2[test2[,"test"] %in% cd4_naive, "manual_metaclusters"]  <- 1
  test2[test2[,"test"] %in% cd4_memory, "manual_metaclusters"]  <- 2
  test2[test2[,"test"] %in% cd8_naive, "manual_metaclusters"]  <- 3
  test2[test2[,"test"] %in% cd8_memory, "manual_metaclusters"]  <- 4
  test2[test2[,"test"] %in% gd_cells, "manual_metaclusters"]  <- 5
  test2[test2[,"test"] %in% nk_cells, "manual_metaclusters"]  <- 6
  test2[test2[,"test"] %in% mono, "manual_metaclusters"]  <- 7
  test2[test2[,"test"] %in% b_cells, "manual_metaclusters"]  <- 8
  
  #test2[,"manual_metaclusters"]
  new_mat <- cbind(exprs(ff),as.vector(test2[,"manual_metaclusters"]))
  colnames(new_mat)[45]  <- "manual_metaclusters"
  #head(new_mat)
  
  ff_new  <- flowFrame(exprs = new_mat)
  ff_new
}

