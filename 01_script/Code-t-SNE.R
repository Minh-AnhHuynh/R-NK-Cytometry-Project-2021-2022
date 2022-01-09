# 1. Setup -------------------------------------------------------------------
# Make sure to install librarian before launching the script.
# Make sure to have opened with an R project file first.
#
# FCS_files are in "./02_data/NK_FCS"
librarian::shelf(flowCore, Rtsne, ggplot2, ggrepel, stringr, uwot)
set.seed(123)

# If help is needed : browseVignettes("flowCore")
# If needed : fs <- read.flowSet(path = "./02_data/NK_FCS/")
files <- list.files("./02_data/NK_FCS/", full.names = TRUE)
FlowSet <- read.flowSet(files)

# Selection of markers of interest
clustering_markers <- read.delim("./02_data/NK_panel.txt")
markers2norm <- clustering_markers[-c(1:7, 18, 24, 29, 39, 43:45), 1]


# Arcsinh markers transformation
ArcTrans <- arcsinhTransform()
TransList <- transformList(from = markers2norm, tfun = ArcTrans)
FlowSet <- TransList %on% FlowSet


# 2. Creation of a data frame with row = cell and column = marker ----------
data <- data.frame()
for (i in 1:length(FlowSet)) {
  print(paste0("fichier numero ", i))
  data_sub <- exprs(FlowSet[[i]])
  data_sub <- data.frame(data_sub)
  data_sub$sample <- basename(files)[i]
  data <- rbind(data, data_sub)
}

markernames <- markernames(FlowSet)
colnames(data) <- c("Time", markernames, "sample")
markerTSNE <- clustering_markers[-c(1:7, 18, 24, 29, 39, 43:45), 2]

# It's possible to make a sample with the following command line to reduce
# time doing the t-SNE (for example 100 000) :
sample_value = 100000
data <-  data[sample(nrow(data), sample_value), ]

# Only sampling will produce tSNE in a meaningful time
data <- data[!colnames(data) %in% "cells"]
forTSNE <- data[, colnames(data) %in% markerTSNE]


# threads = nombre de processeur virtuel
tSNE <-
  Rtsne(
    forTSNE,
    num_threads = 12,
    perplexity = 30,
    max_iter = 1000,
    verbose = TRUE
  )
data$tsne1 <- tSNE$Y[, 1]
data$tsne2 <- tSNE$Y[, 2]

# création d'un tableau des conditions, heures, jour, individu selon chaque échantillon
assignments <- data.frame(
  bc = str_sub(data$sample, 1, 2),
  day = str_sub(data$sample, 3, 5),
  hour = str_sub(data$sample, 6, 8),
  ind = str_sub(data$sample, 10, 14),
  sample = data$sample
)

# changement de l'ordre d'affichage des conditions pour avoir before-prime puis post-prime puis post-boost
data$condition_order <- factor(assignments$bc, levels = c("BP", "PP", "PB"))

# Création du viSNE
plotVISNE <- ggplot(data, aes(x = tsne1, y = tsne2)) +
  geom_point(size = 0.2) +
  scale_color_gradient(low = "yellow", high = "red") +
  xlab("viSNE1") +
  ylab("viSNE2") +
  facet_wrap(~ assignments$bc, ncol = 10) +
  facet_grid(. ~ condition_order) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

pdf(paste0("./05_figures/t-SNE_sample", sample_value, ".pdf"), height = 10, width = 17)
plot(plotVISNE)
dev.off()
