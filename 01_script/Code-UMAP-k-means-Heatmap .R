librarian::shelf(flowCore, ggplot2,ggrepel, stringr,uwot)


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


# sélection des marqueurs d'intérêt
markernames = markernames(FlowSet)
colnames(data) = c("Time",markernames,"sample")
markerUMAP = clustering_markers[-c(1:7,18,24,29,39,43:45),2]

# sample : On prend 5000 parmi le nombre de ligne de la dataframe
# data = data[sample(nrow(data), 100000), ]
data = data[!colnames(data) %in% "cells"]

forUMAP = data[, colnames(data) %in% markerUMAP]

set.seed(123)
### Création de UMAP
UMAP = umap(forUMAP,
            min_dist = 0.5,
            n_threads = 40,
            verbose = TRUE,
            n_sgd_threads = 40
)


forUMAP$umap1 = UMAP[, 1]
forUMAP$umap2 = UMAP[, 2]
forUMAP$sample = data$sample

# création d'un tableau des conditions, heures, jour, individu selon chaque échantillon
assignments = data.frame(bc = str_sub(data$sample,1,2),
                         day = str_sub(data$sample,3,5),
                         hour = str_sub(data$sample,6,8),
                         ind = str_sub(data$sample,10,14),
                         sample = data$sample)

# changement de l'ordre d'affichage des conditions pour avoir before-prime puis post-prime puis post-boost
forUMAP$condition_order = factor(assignments$bc, levels=c('BP','PP','PB'))

plotUMAP = ggplot(forUMAP, aes(x = umap1, y = umap2, color = condition_order)) +
  geom_point(size = 0.2) +
  xlab("umap1") +
  ylab("umap2") +
  facet_wrap(~assignments$bc, ncol = 3) +
  facet_grid(.~condition_order) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )
# geom_point(kmeans@centers)
# km[["centers"]]
plotUMAP
ggsave(filename = paste0("UMAP.pdf"),
       path = "./05_figures")




#########################################
#LOAD DATA
#########################################

##########################################
#PERFORM K-MEANS CLUSTERING
##########################################

#make this example reproducible
set.seed(123)
#perform k-means clustering with k = 4 clusters
k_centers = 10
km <- kmeans(forUMAP [1:31], centers = k_centers, nstart = 20)
#view results
km


forUMAP$clusters = km$cluster
#plot results of final k-means model

plotkmeans = ggplot(forUMAP, aes(x = umap1, y = umap2, color = as.factor(clusters))) +
  geom_point(size = 0.2) +
  xlab("umap1") +
  ylab("umap2") +
  facet_wrap(~assignments$bc, ncol = 10) +
  facet_grid(.~condition_order) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )
plotkmeans
ggsave(
  filename = paste0("kmeans_centers", k_centers, ".pdf"),
  path = "./05_figures"
)




####Représenter les différents marqueurs en fonction des clusters
librarian::shelf(pheatmap, reshape2)

#Heatmap de clusters en fonction de la matrice d'expression
#Moyenne d'expression en fonction de chacun des marqueurs

#Avec melt : on fait un format long => pour chaque ligne (=cellule)), on rajoute son marqueur

melt_forUMAP =  reshape2::melt (forUMAP, id.vars = c("umap1", "umap2", "sample","clusters", "condition_order"), variable.name = "markers", value.name = "expression")


heatmap = reshape2::dcast(melt_forUMAP,
                          melt_forUMAP$clusters ~ melt_forUMAP$markers,
                          value.var = "expression",
                          fun = mean)
row.names(heatmap) = heatmap$`melt_forUMAP$clusters`
heatmap[1]=NULL

pheatmap(
  mat = as.matrix(heatmap),
  filename = paste0(("./05_figures/pheatmap"), k_centers, "clusters.pdf"),
  cluster_cols = TRUE,
  cluster_rows = TRUE,
)


# METACLUSTERS :


# forUMAP$clusters == "1" <- "NKcluster"
# Ensuite on fait color = clusters
# ggline ; ggboxplot
