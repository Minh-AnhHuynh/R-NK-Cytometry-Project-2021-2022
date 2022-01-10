# 1. Setup -------------------------------------------------------------------
# Make sure to install librarian before launching the script.
# Make sure to have opened with an R project file first.
librarian::shelf(flowCore, ggplot2, ggrepel, stringr, uwot, pheatmap, reshape2)
set.seed(123)




# 2. FlowSet data loading ----------------------------------------------------
files <- list.files("./02_data/NK_FCS/", full.names = TRUE)
FlowSet <- read.flowSet(files)

# Selection of markers of interest
clustering_markers <- read.delim("./02_data/NK_panel.txt")
markers2norm <- clustering_markers[-c(1:7, 18, 24, 29, 39, 43:45), 1]

# Arcsinh markers transformation
ArcTrans <- arcsinhTransform()
TransList <- transformList(from = markers2norm, tfun = ArcTrans)
FlowSet <- TransList %on% FlowSet

# Creation of a data frame with row = cell and column = marker
data <- data.frame()
for (i in 1:length(FlowSet)) {
  print(paste0("fichier numero ", i))
  data_sub <- exprs(FlowSet[[i]])
  data_sub <- data.frame(data_sub)
  data_sub$sample <- basename(files)[i]
  data <- rbind(data, data_sub)
}

# Selection of markers of interest
markernames <- markernames(FlowSet)
colnames(data) <- c("Time", markernames, "sample")
markerUMAP <- clustering_markers[-c(1:7, 18, 24, 29, 39, 43:45), 2]

# It's possible to make a sample with the following command line to reduce
# time doing the UMAP (for example 5000) :
#
# data = data[sample(nrow(data), 5000), ]
# The following code takes data from every cell and puts the data into forUMAP
data <- data[!colnames(data) %in% "cells"]
forUMAP <- data[, colnames(data) %in% markerUMAP]




# 3. Creation of UMAP --------------------------------------------------------

UMAP <- umap(forUMAP,
  min_dist = 0.5,
  n_threads = 40,
  verbose = TRUE,
  n_sgd_threads = 40
)
forUMAP$umap1 <- UMAP[, 1]
forUMAP$umap2 <- UMAP[, 2]
forUMAP$sample <- data$sample

# Creation of a table of conditions with times, day, individual for each sample
assignments <- data.frame(
  bc = str_sub(data$sample, 1, 2),
  day = str_sub(data$sample, 3, 5),
  hour = str_sub(data$sample, 6, 8),
  ind = str_sub(data$sample, 10, 14),
  sample = data$sample
)

# Change the order of the conditions to before-prime then post-prime then post-boost
forUMAP$condition_order <- factor(assignments$bc, levels = c("BP", "PP", "PB"))

# Plot and save UMAP
plotUMAP <-
  ggplot(forUMAP, aes(x = umap1, y = umap2, color = condition_order)) +
  geom_point(size = 0.2) +
  xlab(label = "X") +
  ylab(label = "Y") +
  facet_wrap( ~ assignments$bc, ncol = 3) +
  facet_grid(. ~ condition_order) +
  theme_bw() +
  guides(color = guide_legend(title = "Conditions")) +
  theme(
    aspect.ratio = 1,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(
  filename = paste0("UMAP.pdf"),
  path = "./05_figures"
)
ggsave(
  filename = paste0("UMAP.png"),
  path = "./05_figures"
)

# 4.Perform k-means clustering -----------------------------------------------

# Code to detect optimal number of clusters, however vector size is too big.
# librarian::shelf(NbClust)
# NbClust(
# data = as.matrix(forUMAP),
# diss = NULL,
# distance = "euclidean",
# min.nc = 2,
# max.nc = 15,
# method = "kmeans",
# index = "all",
# alphaBeale = 0.1
# )


# Perform k-means clustering with k = 12 clusters according to SPADE heatmap
# metaclusters
k_centers <- 12
km <- kmeans(forUMAP[1:31], centers = k_centers, nstart = 20)

# View results
km

# Put results of k-means into
forUMAP$clusters <- km$cluster

# Plot results of final k-means model
plot_kmeans <- ggplot(forUMAP, aes(x = umap1, y = umap2, color = as.factor(clusters))) +
  geom_point(size = 0.2) +
  xlab(label = "X") +
  ylab(label = "Y") +
  facet_wrap(~ assignments$bc, ncol = 10) +
  facet_grid(. ~ condition_order) +
  theme_bw() +
  guides(color = guide_legend(title = "Cluster")) +
  theme(
    aspect.ratio = 1,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(
  filename = paste0("kmeans_centers", k_centers, ".pdf"),
  path = "./05_figures"
)
ggsave(
  filename = paste0("kmeans_centers", k_centers, ".png"),
  path = "./05_figures"
)



# 5. Represent the different markers according to the clusters ---------------

# Heatmap of clusters according to the expression matrix
# Average expression according to each marker

# With melt : we make a long format => for each line (= cell), we add its
# corresponding marker
melt_forUMAP <- reshape2::melt(forUMAP, id.vars = c("umap1", "umap2", "sample", "clusters", "condition_order"), variable.name = "markers", value.name = "expression")
heatmap <- reshape2::dcast(melt_forUMAP,
  melt_forUMAP$clusters ~ melt_forUMAP$markers,
  value.var = "expression",
  fun = mean
)
row.names(heatmap) <- heatmap$`melt_forUMAP$clusters`
heatmap[1] <- NULL
pheatmap(
  mat = as.matrix(heatmap),
  filename = paste0(("./05_figures/pheatmap"), k_centers, "clusters.pdf"),
  cluster_cols = TRUE,
  cluster_rows = TRUE,
)
pheatmap(
  mat = as.matrix(heatmap),
  filename = paste0(("./05_figures/pheatmap"), k_centers, "clusters.png"),
  cluster_cols = TRUE,
  cluster_rows = TRUE,
)
