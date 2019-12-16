##########################################################################################################################
######################### HIGH DIMENSIONAL CYTOMETRY ANALYSIS WORKFLOW DUNJA MRDJEN ######################################
################## PLEASE CITE THIS PAPER AS A REFERENCE FOR THIS ANALYSIS WORKFLOW AND SCRIPT:###########################
#### Hartmann F., et al., High-dimensional single-cell analysis reveals the immune signature of narcolepsy, JEM, 2016 ####
##########################################################################################################################

#color palette
db13 <- c("#1B9E77","gold" , "#7570B3", "#E7298A", "#D95F02" ,  "lightpink1", "#666666",
          "#57B0FF", "#66CD00", "#1F78B4", "#B2DF8A", "firebrick1", 
          "seagreen1", "brown" , "#6A3D9A", "mediumorchid", "grey" , "gold", "gold",
          "lightpink1", "#BC80BD", "gold", "#80B1D3", "gold","#80B1D3", "mediumorchid")


# Set the working directory (where outputs will be and where the files are)
setwd("/Volumes/Neuroimmunology/People/Mrdjen/BS 6-1-17 Astro FACS/Dunja Analysis_2017.09.27/")
load(".RData")

####################################################
##### reading the merged or single fcs file(s) #####
####################################################

# load the fcs-files into a flowSet (fs)
fs <- read.flowSet(path = ("transformed/"), pattern = "*.fcs", alter.names = T)
fs
sampleNames(fs)
colnames(fs)

# combine the flowset into a single flowframe with an additional gate.source channel
dim <- fsApply(fs, dim)
dim <- as.numeric(dim[,1])
dim
gate.source <- as.vector(x = NULL)
for(i in 1:length(dim)) {temp.source <- rep(i, dim[i])
gate.source <- c(gate.source, temp.source)}

# combine data into a matrix
data <- fsApply(fs, exprs)
data <- cbind(data, gate.source)
data <- data.frame(data)
table(data$gate.source)
colnames(data)

#change channel names
channel.names <- read.xls("channel.names.xlsx")
channel.names
new.names <- channel.names$new.names
colnames(data) <- new.names
colnames(data)
range(data)
range(data[,"gate.source"])

#OR
#read in a single fcs file
ffm <- read.FCS("final_merged_export_0.fcs", transformation = F)
colnames(ffm)
data <- flowCore::exprs(ffm)
params <- flowCore::parameters(ffm)
desc <- flowCore::description(ffm)

#########################################################################################
###### PRE-PROCESSING #### HIGHLY INFLUENCED BY THE OUTLIERS IN YOUR DATA ###############
#########################################################################################

data.original <- data
data <- data.original

# do arcsinh transformation only for the clustering.channels
asinh_scale <- 50
data[,"GS"] <- asinh(data[,"GS"] / asinh_scale)

# define clustering columns #(panel)
data <- data.frame(data)
dim(data)
colnames(data)
panel <- colnames(data)[1:10]
panel
length(panel)

#ONLY FOR FACS DATA: adjust data to start from (approximately) 0
q.vector <- apply(data[,panel], 2, function(x) quantile(x, 0.0001, names = F))
q.vector
data.shift <- data
data.shift[,panel] <- sweep(data.shift[,panel], 2, q.vector)

# check the min values
apply(data.shift[,panel], 2, min)
apply(data[,panel], 2, min)

#PERCENTILE normalize expression values from 0 to 1
per.vector <- apply(data.shift[,panel], 2, function(x) quantile(x, 0.9999, names = F))
per.vector
data.shift[,panel] <- t(t(data.shift[,panel]) / as.numeric(per.vector))

# check whether you adjusted the range approximately from 0 to 1
apply(data.shift[,panel], 2, max)

# show biaxial plot
ggplot(data = data.frame(data.shift[1:100000,]), aes(x = ACSA2, y = GS)) + geom_point(alpha=1, size=0.5) +
  theme_biaxial #+ xlim(0,1) + ylim(0,1)

# replace the un-processed data with data.shift
data <- data.shift

# export FCS file for double checking in FlowJo
data <- data.matrix(data)
write.FCS(x = flowFrame(data), filename = "data.normed.fcs")
#save norm'd data as robject
saveRDS(data, file = "data.normed.robject", compress = F)

#########################################################################################

##### read in excel with sample.if #####
md <- read.xls("sample.id.xlsx")
md

# merge data with sample.id tags
data <- merge(data, md, by = "gate.source")
head(data)
table(data$gate.source)

#### Important: give each cell an cell.id to keep track of #####
data <- as.data.frame(data)
data$cell.id <- 1:nrow(data)
str(data)
dim(data)
head(data)

#########################################################################################
########################### CALCULATE tSNE ##############################################
#########################################################################################

dim(data)
colnames(data)

# create subsample vector
n_sub <- 20000
n <- nrow(data)
set.seed(123)
ix <- sample(1:n, n_sub)

# check/redefine panel
panel

# prepare data for Rtsne
data_rtsne <- data[ix, c(panel)]
data_rtsne <- as.matrix(data_rtsne)
head(data_rtsne)
colnames(data_rtsne)
dim(data_rtsne)

# run bh SNE
set.seed(123)
out_rtsne <- Rtsne(data_rtsne, dims = 2, perplexity = 50, theta = 0.5, #Run Rtnse.multicore if using Ubuntu/S3IT
                   max_iter = 1000, verbose = T, pca = F, check_duplicates=F)


# prepare for plotting
data_plot <- as.data.frame(out_rtsne$Y)
colnames(data_plot) <- c("bhSNE1", "bhSNE2")
data_plot$cell.id <- data[ix,"cell.id"]
data_plot$gate.source <- data[ix,"gate.source"]
head(data_plot)


# keep track of data_plot for all cells
data_plot1 <- data_plot
data_plot2 <- data_plot

# plot tSNE
p1 <- ggplot(data_plot_wt, aes(x = bhSNE1, y = bhSNE2)) +
  geom_point(size = 1) + 
  coord_fixed(ratio = 1)
p1 + theme_tsne

# if it's good then save data_plot as an R Object to keep it safe
saveRDS(data_plot, file = "data_plot69.robject", compress = F)

ggsave(filename = "data_plot2.png", plot = p1 + theme_tsne, 
       scale = 1, width = 5, height = 5, units = c("in"))

#separate data_plot by group
data_plot_md <- merge(data_plot, md, by = "gate.source")
data_plot_wt <- droplevels(subset(data_plot_md, group == "wt"))
data_plot_ko <- droplevels(subset(data_plot_md, group == "ko"))

# combine with other expression values
tsne <- merge(data_plot[,1:3], data, by = "cell.id")

# plot tSNE
p2 <- ggplot(tsne, aes(x = bhSNE1, y = bhSNE2)) +
  geom_point(size = 0.5, aes(colour = S100b)) + 
  scale_colour_gradient2(low = "black", mid = "yellow", high = "red", midpoint = 0.5, limits = c(0,1), 
                         space = "rgb", guide = "colourbar") +      
  coord_fixed(ratio = 1)
p2 + theme_tsne

ggsave(filename = "tsne_S100b.png", plot = p2 + theme_tsne, 
       scale = 1, width = 6, height = 6, units = c("in"))

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



#########################################################################################
########################### FLOWSOME CLUSTERING #########################################
#########################################################################################
colnames(data)
data <- data.matrix(data)

# run FlowSOM (with set.seed for reproducibility)
set.seed(123)
out_fSOM <- FlowSOM::ReadInput(flowFrame(exprs = data, desc = list(FIL = 1)), transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = panel)
out_fSOM <- FlowSOM::BuildMST(out_fSOM)
labels <- out_fSOM$map$mapping[,1]

# set the max value for k, which is the maximum metaclusters you want
max  <- 20

# initialize matrix
gate_source <- as.factor(data[,"gate.source"])
meta_results <- data.frame(gate_source)

# do a manual metaclustering for all values up to max
for (i in 3:max) {
  set.seed(123)
  out_meta <- FlowSOM::metaClustering_consensus(out_fSOM$map$codes, k = i)
  meta_results <- cbind(meta_results, as.factor(out_meta[labels]))}

meta_results <- meta_results[,2:ncol(meta_results)]
colnames(meta_results) <- paste("k.", 3:max, sep = "")
meta_results <- cbind(meta_results, data[,"cell.id"])
colnames(meta_results)[19] <- "cell.id"
colnames(meta_results)

#keep track of different meta_results


# prepare the joined metaclustering data
library(reshape2)
meta.melt <- melt(meta_results, id.vars = "cell.id", variable.name = "k.value", value.name = "cluster.assigment")
meta.melt$cluster.assigment <- as.factor(meta.melt$cluster.assigment)
joined.meta <- merge(meta.melt, data_plot, by = "cell.id")
joined.meta <- merge(joined.meta, md, by = "gate.source")

# plot tSNEs with k.values overlayed
p4 <- ggplot(joined.meta, aes(x = bhSNE1, y = bhSNE2, color = cluster.assigment)) +
  geom_point(size = 0.25) +
  coord_fixed(ratio = 1) +
  scale_colour_manual(name = NULL,values = c(db13)) +
  facet_wrap(~ k.value, ncol =6, scales = "free") +
  guides(colour = guide_legend(override.aes = list(size=5), title="Cluster"))

p4 + theme_tsne

# save the plots as png
ggsave(filename = "clust_data.eae.png", plot = p3 + theme_tsne, 
       scale = 1, width = 20, height = 20, units = c("in"))

# plot a specific k.value for all cells
data_plot_clust <- merge(meta_results, data_plot, by = "cell.id")
data_plot_clust <- merge(data_plot_clust, md, by = "gate.source")
dim(data_plot_clust)

p5 <- ggplot(data_plot_clust, aes(x = bhSNE1, y = bhSNE2, color = as.factor("k.x"))) +
  geom_point(size = 1) +
  coord_fixed(ratio = 1) + 
  scale_colour_manual(name = NULL, 
                      values = c(db13.4)) +
  guides(colour = guide_legend(override.aes = list(size=2), title="Cluster"))
p5 + theme_tsne5

# save the plots as png
ggsave(filename = "Steady state CD45positive/data.ss.young_DCs.png", plot = p5 + theme_tsne5, 
       scale = 1, width = 8, height = 8, units = c("in"))


# plot with a contour backgroung
p6 <- ggplot(data = clust_old, aes(x = bhSNE1, y = bhSNE2)) +
  geom_density_2d(data = data_plot_clust, aes(x = bhSNE1, y = bhSNE2), colour = "lightgrey", size = 0.5, bins=20) +
  coord_fixed(ratio = 1) + 
  ylim(-40,40) + xlim(-40,40) +
  geom_point(data = clust_old, size = 0.5, aes(colour = as.factor(manual.metaclusters))) +
  scale_colour_manual(name = NULL, values = c(db13)) +
  guides(colour = guide_legend(override.aes = list(size=2), title="Cluster"))
p6 + theme_tsne5

# for an interactive t-sne map:
library(plotly)
ggplotly()

# save the plots as png
ggsave(filename = "Steady state CD45positive/data.ss.young_DCs.png", plot = p5 + theme_tsne5, 
       scale = 1, width = 8, height = 8, units = c("in"))

#########################################################################################
########################### MAKE HEATMAPS ###############################################
#########################################################################################

# define heatmap color
heat_col  <- "black" #for black and white
my_color <- colorRampPalette(c("white",heat_col))(n = 300) #black and white
my_color <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)) #blue and red

# make combined expression matrix data.ss
ce <- merge(data, meta_results, by = "cell.id")
ce <- data.matrix(ce)

# plot heatmap
cols.to.plot <- panel
heat_fun(ce = data.matrix(ce), cc = "k.x", 
         clustering = cols.to.plot, minimum_cells = 10, color = my_color)


#################################
##### manual metaclustering #####
#################################
colnames(meta_results)
colnames(data)

cluster.column <- "k.x"
mmc <- merge(data.ss, meta_results10[, c("k.x", "cell.id")], by = "cell.id")
colnames(mmc)
colnames(mmc)[38] <- cluster.column
manual.metaclusters <- as.vector(rep(0, nrow(mmc)))
mmc2 <- cbind(mmc, manual.metaclusters)
head(mmc2)

#for data.all
meta1 <- c(1) #Microglia
meta2 <- c(2,3) #Macrophages
meta3 <- c(10) #Neutrophils
meta4 <- c(4,6) #Monocytes
meta5 <- c(5,8) #cDCs
meta6 <- c(11) #pDCs
meta7 <- c(14,12) #B cells
meta8 <- c(13,15) #T cell, NK cells, ILCs
meta9 <- c(9) #Eosinophils
meta10 <- c(7) #unidentified

#give them the right metacluster number
mmc2[mmc2[,cluster.column] %in% meta1, "manual.metaclusters"]  <- 1
mmc2[mmc2[,cluster.column] %in% meta2, "manual.metaclusters"]  <- 2
mmc2[mmc2[,cluster.column] %in% meta3, "manual.metaclusters"]  <- 3
mmc2[mmc2[,cluster.column] %in% meta4, "manual.metaclusters"]  <- 4
mmc2[mmc2[,cluster.column] %in% meta5, "manual.metaclusters"]  <- 5
mmc2[mmc2[,cluster.column] %in% meta6, "manual.metaclusters"]  <- 6
mmc2[mmc2[,cluster.column] %in% meta7, "manual.metaclusters"]  <- 7
mmc2[mmc2[,cluster.column] %in% meta8, "manual.metaclusters"]  <- 8
mmc2[mmc2[,cluster.column] %in% meta9, "manual.metaclusters"]  <- 9
mmc2[mmc2[,cluster.column] %in% meta10, "manual.metaclusters"]  <- 10

colnames(mmc2)
table(mmc2$manual.metaclusters)

#keep track of mmc2 versions
mmc2.ss <- mmc2

###############################
###### subset clusters ########
###############################


################################
##### frequency analysis #######
################################
ce <- merge(data, meta_results, by = "cell.id")
cc = "k.x"

# make frequency table
fm <- freq_fun_mat(data = data.matrix(ce), source.column = "gate.source", cluster.column = cc, minimum.cells = 2)
rownames(fm) <- 1:x

# now we need to make a dataframe for plotting
dimnames <- list(gate.source = rownames(fm), cluster = colnames(fm))
mat <- matrix(fm, ncol = ncol(fm), nrow = nrow(fm), dimnames = dimnames)
df <- as.data.frame(as.table(mat))
df <- merge(df, md, by="gate.source")
df_plot <- df

#rename the clusters
colnames(df_plot) <- c("name1", "name2")

# apply t-test for all clusters but you can only do this if you have 3 or more samples per group
tres <- ddply(df_plot, "cluster", function(x) {
  w <- t.test(Freq ~ group, data = x)
  with(w,data.frame(statistic,method,p.value))
})
tres

# correct for multiple comparisons using the Benjamini-Hochberg test
tres$corrected.p.value <- p.adjust(tres$p.value, method="BH")
tres <- tres[with(tres, order(p.value)), ]
tres$rounded.p.value <- round(tres$corrected.p.value, 4)
tres

# write csv file
write.csv(tres, file="ttest_name.csv",eol="\r")

# make an overview plot
p7 <- ggplot(plot, aes(x = age, y = Freq, ymin = 0, ymax = Freq*1.3)) +
  facet_wrap(~ cluster, ncol=3, scales="free") +
  geom_point(size=1) +
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot()
p7 + theme_pub


# save the plots as pdf
ggsave(filename = "Macrophages/Young and Old/FREQ_MACs_Young_VS_Old_k.3_stats.pdf", plot = p4 + theme_pub, useDingbats = F, 
       scale = 1, width = 6*1.5, height = 3*3.5, units = c("in"))

# print summary statitics in a table
sum <-summaryBy(. ~ cluster,  
                data = df_plot, 
                FUN = function(x) { c(m = median(x, na.rm = TRUE), s = sd(x, na.rm = TRUE), n = length(x)) } )
sum



#########################################################################################
########################### SAVE AN FCS FILE ############################################
#########################################################################################


#you will need the data called data.out
data.out <- data.matrix(mmc2[,])
ff.out <- flowFrame(exprs = data.out)
write.FCS(ff.out, filename = "name.fcs")


# save the R environment
save.image(file = "data_astrocytes_BS.rData", compress = T)


