# FCS_files are in "./02_data/NK_FCS"
librarian::shelf(flowCore, Rtsne, ggplot2,ggrepel, stringr,uwot)
# browseVignettes("flowCore")


setwd("./R_cytometry_project")

# fs <- read.flowSet(path = "./02_data/NK_FCS/")
files = list.files("./02_data/NK_FCS/", full.names = TRUE)
FlowSet = read.flowSet(files)

# sélection des marqueurs d'intérêt
clustering_markers = read.delim('./02_data/NK_panel.txt')
markers2norm = clustering_markers[-c(1:7,18,24,29,39,43:45),1]


# transformation arcsinh des marqueurs
ArcTrans = arcsinhTransform()
TransList = transformList(from = markers2norm, tfun = ArcTrans)
FlowSet = TransList %on% FlowSet


# création d'une dataframe avec ligne = cellule et colonne = marqueur
data = data.frame()
for (i in 1:length(FlowSet)) {
  print(paste0("fichier numero ", i))
  data_sub <- exprs(FlowSet[[i]])
  data_sub <- data.frame(data_sub)
  data_sub$sample <- basename(files)[i]
  data <- rbind(data, data_sub)
}

markernames = markernames(FlowSet)
colnames(data) = c("Time",markernames,"sample")

markerTSNE <- clustering_markers[-c(1:7,18,24,29,39,43:45),2]

# sample : On prend 5000 parmi le nombre de ligne de la dataframe
data <- data[sample(nrow(data), 5000), ]
data <- data[!colnames(data) %in% "cells"]

forTSNE = data[, colnames(data) %in% markerTSNE]


# threads = nombre de processeur virtuel
tsne <- Rtsne(forTSNE, num_threads = 10)
data$tsne1 <- tsne$Y[, 1]
data$tsne2 <- tsne$Y[, 2]

# création d'un tableau des conditions, heures, jour, individu selon chaque échantillon
assignments = data.frame(bc = str_sub(data$sample,1,2),
                           day = str_sub(data$sample,3,5),
                           hour = str_sub(data$sample,6,8),
                           ind = str_sub(data$sample,10,14),
                           sample = data$sample)

# changement de l'ordre d'affichage des conditions pour avoir before-prime puis post-prime puis post-boost
data$condition_order = factor(assignments$bc, levels=c('BP','PP','PB'))

# Création du viSNE
plotVISNE = ggplot(data, aes(x = tsne1, y = tsne2)) +
                    geom_point(size = 0.2) +
                    scale_color_gradient(low = "yellow", high = "red") +
                    xlab("viSNE1") +
                    ylab("viSNE2") +
                    facet_wrap(~assignments$bc, ncol = 10) +
                    facet_grid(.~condition_order) +
                    theme_bw() +
                    theme(aspect.ratio = 1,
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()
                  )

plot(plotVISNE)