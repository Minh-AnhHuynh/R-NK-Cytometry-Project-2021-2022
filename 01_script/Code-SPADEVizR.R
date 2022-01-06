# 1. Setup -------------------------------------------------------------------
# Make sure to install librarian before launching the script.
# Make sure to have opened with an R project file first.
librarian::shelf(tchitchek - lab / SPADEVizR, data.table, ggplot2, stringr, dplyr, flowCore)
set.seed(123) # for reproducibility
seed_number <- 123





# 2. Data loading ------------------------------------------------------------
clustering_markers <- read.delim("./02_data/NK_panel.txt")
# Filtering markers
excluded_markers <- clustering_markers[c(1:7, 18, 24, 29, 39, 43, 44, 45), 1]
excluded_markers2 <- c("File Number", "density", "cells-(Ir191)Di", "cells-(Ir193)Di")
excluded_markers <- c(excluded_markers, excluded_markers2)

## Select k_value
k_value <- 100

## Load SPADE results
output_dir <- paste0("./03_spade_analysis/spade_k", k_value)
spade_results <- importResultsFromSPADE(output_dir,
  exclude.markers = excluded_markers
)
# To check if spade_results is correct :
# head(spade_results@cluster.phenotypes)





# 3. Quality Control ---------------------------------------------------------

# Generate a QC report to detect cell cluster phenotype with uniform marker distribution expressions
accuracyQC <- qcUniformClusters(spade_results,
  density.PDFfile = paste0("./04_SPADEVizR/QualityControls/UniformClusters_density_k", k_value, "_seed", seed_number, ".pdf"),
  heatmap.PDFfile = paste0("./04_SPADEVizR/QualityControls/UniformClusters_heatmap_k", k_value, "_seed", seed_number, ".pdf"),
  uniform.test = "unimodality"
)

# Display the percentage of accuracy clusters
print(accuracyQC$perc)

# Save the accuracy matrix in a text file
write.table(accuracyQC$accuracy.matrix, paste0("./04_SPADEVizR/QualityControls/UniformClusters-accuracy-matrix_k_", k_value, ".txt"), quote = FALSE, sep = "\t", col.names = NA)

# Generate a QC report to detect cell cluster phenotypes with low number of cell clusters
smallclusterQC <- qcSmallClusters(spade_results, PDFfile = paste0("./04_SPADEVizR/QualityControls/QCreport-SmallClusters_heatmap_k", k_value, ".pdf"), th = 500)

# Display the percentage of cell cluster phenotypes with low number of cell clusters
print(smallclusterQC$perc)

# Save the small cluster matrix in a text file
write.table(smallclusterQC$small.clusters, paste0("./04_SPADEVizR/QualityControls/UniformClusters-smallclusters-matrix_k", k_value, ".txt"), quote = FALSE, sep = "\t", col.names = NA)





# 4. Annotation of metadata information -----------------------------------


## 4.1.1 Data frame creation : row = cell ; col = markers =================
##  This step is done in order to get the assignment data.frame without selecting samples manually

# Load data using flowSet
files <- list.files("./02_data/NK_FCS/", full.names = TRUE)
FlowSet <- read.flowSet(files)

# Select markers of interest
clustering_markers <- read.delim("./02_data/NK_panel.txt")
markers2norm <- clustering_markers[-c(1:7, 18, 24, 29, 39, 43:45), 1]

# Arcsinh transformation of markers to be able to use Flowset
ArcTrans <- arcsinhTransform()
TransList <- transformList(from = markers2norm, tfun = ArcTrans)
FlowSet <- TransList %on% FlowSet

# Generation of files
data <- data.frame()
for (i in 1:length(FlowSet)) {
  print(paste0("fichier numero ", i))
  data_sub <- exprs(FlowSet[[i]])
  data_sub <- data.frame(data_sub)
  data_sub$sample <- basename(files)[i]
  data <- rbind(data, data_sub)
}




## 4.2.2 Assigning sample biological conditions (metadata information) ####
assignments <- data.frame(
  bc = str_sub(data$sample, 1, 8),
  tp = str_sub(data$sample, 1, 8),
  # hour = str_sub(data$sample,6,8),
  ind = str_sub(data$sample, 10, 14),
  sample_names = str_sub(data$sample, 1, 23)
)

# Filtering samples to obtain only 42 values corresponding to 42 samples
# Note : This uses the dplyr packages which converts a dataframe into a tibble to work with.
assignments <- assignments %>%
  group_by(sample_names) %>%
  slice(1)

# Filtering according to timepoints : BP, PP, BP
df1 <- as_tibble(assignments)
df2 <- slice(df1, 1:4, 22:42, 5:21)
assignments <- df2

# Put sample names in row names
sample <- assignments$sample_names
assignments[4] <- NULL
row.names(assignments) <- sample

# Put assignment data.frame tidied into spade_results
spade_results <- assignContext(spade_results, assignments = assignments)

# Transformation of the dataframe assignments from a tibble into a data.frame to make MDSViewer work
spade_results@assignments <- as.data.frame(spade_results@assignments)




## 4.2.3 Select samples according to their conditions respectively ========

# Before Prime
condition_BP <- c("BPD19H00_BB078_CD3-CD8+", "BPD19H00_BB231_CD3-CD8+", "BPD19H00_BC641_CD3-CD8+", "BPD19H00_BD620_CD3-CD8+")

test <- assignments %>% select(starts_with("BP"))
test <- assignments %>% tidyr::pivot_longer(sample_names, names_to = NULL, values_to = "sample_names", )


condition_BP = filter(assignments, grepl("BP", assignments$sample_names))
Mydata1 = filter(mydata, grepl(0,hp))

condition_BP = row.names(assignments)[,grepl("BP", assignments$sample_names) == TRUE))


test <- filter(assignments, grepl("BP", row.names(assignments)))
test <- na.omit (test)
grepl("BP", row.names(assignments))


condition_BP = sample[sample %in% grepl("BP", rownames(assignments))]


# Post Prime
condition_PP <- c("PPD00H00_BB078_CD3-CD8+", "PPD00H00_BB231_CD3-CD8+", "PPD00H00_BD620_CD3-CD8+", "PPD00H03_BB078_CD3-CD8+", "PPD00H03_BB231_CD3-CD8+", "PPD00H03_BC641_CD3-CD8+", "PPD00H03_BD620_CD3-CD8+", "PPD00H06_BB078_CD3-CD8+", "PPD00H06_BB231_CD3-CD8+", "PPD00H06_BC641_CD3-CD8+", "PPD01H00_BB078_CD3-CD8+", "PPD01H00_BB231_CD3-CD8+", "PPD01H00_BC641_CD3-CD8+", "PPD01H00_BD620_CD3-CD8+", "PPD03H00_BB078_CD3-CD8+", "PPD03H00_BB231_CD3-CD8+", "PPD03H00_BC641_CD3-CD8+", "PPD14H00_BB078_CD3-CD8+", "PPD14H00_BB231_CD3-CD8+", "PPD14H00_BC641_CD3-CD8+", "PPD14H00_BD620_CD3-CD8+")

# Post Boost
condition_PB <- c("PBD00H00_BB078_CD3-CD8+", "PBD00H00_BB231_CD3-CD8+", "PBD00H00_BD620_CD3-CD8+", "PBD00H03_BB078_CD3-CD8+", "PBD00H03_BB231_CD3-CD8+", "PBD00H03_BC641_CD3-CD8+", "PBD00H06_BB078_CD3-CD8+", "PBD00H06_BB231_CD3-CD8+", "PBD00H06_BC641_CD3-CD8+", "PBD00H06_BD620_CD3-CD8+", "PBD01H00_BB078_CD3-CD8+", "PBD01H00_BB231_CD3-CD8+", "PBD01H00_BC641_CD3-CD8+", "PBD01H00_BD620_CD3-CD8+", "PBD03H00_BB078_CD3-CD8+", "PBD03H00_BB231_CD3-CD8+", "PBD03H00_BD620_CD3-CD8+")




## 4.2.4 Cluster Manipulation  ============================================
# Merge the abundances and the phenotypes of clusters into a new cluster in a Results object
M.I = c("9", "24", "62", "36", "28", "35", "59", "77")
M.II = c("3", "31", "11", "68", "84", "69", "49", "96")
M.III = c("20", "74", "83", "33", "72", "22", "53")
M.IV = c("26", "94", "2", "91","38", "43", "58", "75", "89", "55", "87", "76", "65", "81")
M.V = c("27", "6", "82", "63", "79", "52", "54")
M.VI = c("42", "90", "88", "78", "34", "71", "70", "80", "8", "10", "1", "18")
M.VII = c("97", "44", "100", "41", "50", "95")
M.VIII = c("60", "25", "61", "17", "23")
M.IX = c("51", "30", "40", "98")
M.X = c("92", "47", "64", "57", "48", "66", "32", "85", "13", "99", "21", "37", "39", "46", "7", "19", "12", "29", "14", "16", "4", "56")
M.XI = c("73", "93", "86", "45", "67")
M.XII = c("5", "15")
cluster_results <- mergeClusters(spade_results, clusters = M.I, name = "M.I")
cluster_results <- mergeClusters(cluster_results, clusters = M.II, name = "M.II")
cluster_results <- mergeClusters(cluster_results, clusters = M.III, name = "M.III")
cluster_results <- mergeClusters(cluster_results, clusters = M.IV, name = "M.IV")
cluster_results <- mergeClusters(cluster_results, clusters = M.V, name = "M.V")
cluster_results <- mergeClusters(cluster_results, clusters = M.VI, name = "M.VI")
cluster_results <- mergeClusters(cluster_results, clusters = M.VII, name = "M.VII")
cluster_results <- mergeClusters(cluster_results, clusters = M.VIII, name = "M.VIII")
cluster_results <- mergeClusters(cluster_results, clusters = M.IX, name = "M.IX")
cluster_results <- mergeClusters(cluster_results, clusters = M.X, name = "M.X")
cluster_results <- mergeClusters(cluster_results, clusters = M.XI, name = "M.XI")
cluster_results <- mergeClusters(cluster_results, clusters = M.XII, name = "M.XII")
print(cluster_results@cluster.names)

# Delete clusters
cluster_results <- removeClusters(cluster_results, clusters = c("9", "24", "62", "36", "28", "35", "59", "77"))
cluster_results <- removeClusters(cluster_results, clusters = c("3", "31", "11", "68", "84", "69", "49", "96"))
cluster_results <- removeClusters(cluster_results, clusters = c("20", "74", "83", "33", "72", "22", "53"))
cluster_results <- removeClusters(cluster_results, clusters = c("26", "94", "2", "91","38", "43", "58", "75", "89", "55", "87", "76", "65", "81"))
cluster_results <- removeClusters(cluster_results, clusters = c("27", "6", "82", "63", "79", "52", "54"))
cluster_results <- removeClusters(cluster_results, clusters = c("42", "90", "88", "78", "34", "71", "70", "80", "8", "10", "1", "18"))
cluster_results <- removeClusters(cluster_results, clusters = c("97", "44", "100", "41", "50", "95"))
cluster_results <- removeClusters(cluster_results, clusters = c("60", "25", "61", "17", "23"))
cluster_results <- removeClusters(cluster_results, clusters = c("51", "30", "40", "98"))
cluster_results <- removeClusters(cluster_results, clusters = c("92", "47", "64", "57", "48", "66", "32", "85", "13", "99", "21", "37", "39", "46", "7", "19", "12", "29", "14", "16", "4", "56"))
cluster_results <- removeClusters(cluster_results, clusters = c("73", "93", "86", "45", "67"))
cluster_results <- removeClusters(cluster_results, clusters = c("5", "15"))
print(cluster_results@cluster.names)


# Keep only NK clusters

NK_results = removeClusters(spade_results, clusters = c("9", "24", "62", "36", "28", "35", "59", "77"))
NK_results <- removeClusters(NK_results, clusters = c("3", "31", "11", "68", "84", "69", "49", "96"))
NK_results <- removeClusters(NK_results, clusters = c("97", "44", "100", "41", "50", "95"))
NK_results <- removeClusters(NK_results, clusters = c("60", "25", "61", "17", "23"))
NK_results <- removeClusters(NK_results, clusters = c("51", "30", "40", "98"))
NK_results <- removeClusters(NK_results, clusters = c("92", "47", "64", "57", "48", "66", "32", "85", "13", "99", "21", "37", "39", "46", "7", "19", "12", "29", "14", "16", "4", "56"))
NK_results <- removeClusters(NK_results, clusters = c("73", "93", "86", "45", "67"))
NK_results <- removeClusters(NK_results, clusters = c("5", "15"))
print(cluster_results@cluster.names)



## 4.2.5 Cluster Annotations ===============================================
# defines an annotation dataframe
annotations <- data.frame()
annotations["M.I", "CD3"] <- "c(9, 24, 62, 36, 28, 35, 59, 77)"
annotations["M.II", "CD3"] <- "c(3, 31, 11, 68, 84, 69, 49, 96)"
annotations["M.III", "GranzymeB"] <- "c(20, 74, 83, 33, 72, 22, 53)"
annotations["M.IV", "GranzymeB"] <- "c(26, 94, 2, 91,38, 43, 58, 75, 89, 55, 87, 76, 65, 81)"
annotations["M.V", "GranzymeB"] <- "c(27, 6, 82, 63, 79, 52, 54)"
annotations["M.VI", "GranzymeB"] <- "c(42, 90, 88, 78, 34, 71, 70, 80, 8, 10, 1, 18)"
annotations["M.VII - Monocytes", "CD14"] <- "c(97, 44, 100, 41, 50, 95)"
annotations["M.VIII", "CD11c"] <- "c(60, 25, 61, 17, 23)"
annotations["M.IX", "CD66"] <- "c(51, 30, 40, 98)"
annotations["M.X", "CD66"] <- "c(92, 47, 64, 57, 48, 66, 32, 85, 13, 99, 21, 37, 39, 46, 7, 19, 12, 29, 14, 16, 4, 56)"
annotations["M.XI", "CD66"] <- "c(73, 93, 86, 45, 67)"
annotations["M.XII", "CD20"] <- "c(5, 15)"
print(annotations)

# annotates the cell clusters in a Results object#
# cell clusters are renamed according to the population names
cluster_results <- annotateClusters(spade_results, annotations = annotations)
print(cluster_results@cluster.names)




# 5. Analysis and Figures ----------------------------------------------------

## 4.2 Heatmap generation ====================================================
# General heatmap on 100 clusters
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/heatmap_spade_k", k_value, ".pdf"), height = 10, width = 17)
heatmapViewer(spade_results)
dev.off()

# Focused heatmap on metaclusters
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/heatmap_spade_k", k_value, "_NK-clusters.pdf"), height = 10, width = 17)
heatmapViewer(NK_results)
dev.off()




## 4.3 Identification of Differential Abundant Clusters and Volcano Plot ====
# Volcano plot Condition BP vs PP
resultsDAC_BPvsPP <- identifyDAC(
  spade_results,
  condition1 = condition_BP,
  condition2 = condition_PP,
  th.pvalue = 0.05,
  th.fc = 2,
  method.paired = FALSE
)
export(resultsDAC_BPvsPP, filename = paste0("./04_SPADEVizR/ResultsDAC/resultsDAC_k", k_value, "_BPvsPP.txt"))
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/VolcanoDAC_k", k_value, "_BPvsPP.pdf"), height = 10, width = 17)
volcanoViewer(resultsDAC_BPvsPP)
dev.off()

# Volcano plot Condition PP v PB
resultsDAC_PPvsPB <- identifyDAC(
  spade_results,
  condition1 = condition_PP,
  condition2 = condition_PB,
  th.pvalue = 0.05,
  th.fc = 2,
  method.paired = FALSE
)
export(resultsDAC_PPvsPB, filename = paste0("./04_SPADEVizR/ResultsDAC/resultsDAC_k", k_value, "_PPvsPB.txt"))
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/VolcanoDAC_k", k_value, "_PPvsPB.pdf"), height = 10, width = 17)
volcanoViewer(resultsDAC_PPvsPB)
dev.off()

# Volcano plot Condition BP vs PB
resultsDAC_BPvsPB <- identifyDAC(
  spade_results,
  condition1 = condition_BP,
  condition2 = condition_PB,
  th.pvalue = 0.05,
  th.fc = 2,
  method.paired = FALSE
)
export(resultsDAC_BPvsPB, filename = paste0("./04_SPADEVizR/ResultsDAC/resultsDAC_k", k_value, "_BPvsPB.txt"))
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/VolcanoDAC_k", k_value, "_BPvsPB.pdf"), height = 10, width = 17)
volcanoViewer(resultsDAC_BPvsPB)
dev.off()

# Volcano plot for NK clusters only, condition BP vs PP
resultsDAC_BPvsPP_NK <- identifyDAC(
  NK_results,
  condition1 = condition_BP,
  condition2 = condition_PP,
  th.pvalue = 0.05,
  th.fc = 2,
  method.paired = FALSE
)
export(resultsDAC_BPvsPP_NK, filename = paste0("./04_SPADEVizR/ResultsDAC/resultsDAC_k", k_value, "_BPvsPP_NK-clusters.txt"))
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/VolcanoDAC_k", k_value, "_BPvsPP_NK-clusters.pdf"), height = 10, width = 10)
volcanoViewer(resultsDAC_BPvsPP_NK)
dev.off()

# Volcano plot for NK clusters only, condition PP vs PB
resultsDAC_PPvsPB_NK <- identifyDAC(
  NK_results,
  condition1 = condition_PP,
  condition2 = condition_PB,
  th.pvalue = 0.05,
  th.fc = 2,
  method.paired = FALSE
)
export(resultsDAC_PPvsPB_NK, filename = paste0("./04_SPADEVizR/ResultsDAC/resultsDAC_k", k_value, "_PPvsPB_NK-clusters.txt"))
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/VolcanoDAC_k", k_value, "_PPvsPB_NK-clusters.pdf"), height = 10, width = 10)
volcanoViewer(resultsDAC_PPvsPB_NK)
dev.off()

# Volcano plot for NK clusters only, condition BP vs PB
resultsDAC_BPvsPB_NK <- identifyDAC(
  NK_results,
  condition1 = condition_BP,
  condition2 = condition_PB,
  th.pvalue = 0.05,
  th.fc = 2,
  method.paired = FALSE
)
export(resultsDAC_BPvsPB_NK, filename = paste0("./04_SPADEVizR/ResultsDAC/resultsDAC_k", k_value, "_PPvsPB_NK-clusters.txt"))
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/VolcanoDAC_k", k_value, "_BPvsPB_NK-clusters.pdf"), height = 10, width = 10)
volcanoViewer(resultsDAC_BPvsPB_NK)
dev.off()

## 4.4 Tree Viewer ============================================================

# displays in a pdf the SPADE tree for sample PBD28_BB078
# pdf("./05_figures/treeViewer-PBD28_BB078.pdf", width=15, height=15)
# treeViewer(spade_results,samples="PBD28_BB078")
# dev.off()

# displays in a pdf the aggregated SPADE tree for all samples
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/treeViewer_k", k_value, "_GranzymeB.pdf"), width = 15, height = 15)
treeViewer(spade_results, marker = "GranzymeB")
dev.off()

pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/treeViewer_k", k_value, "_CD8.pdf"), width = 15, height = 15)
treeViewer(spade_results, marker = "CD8")
dev.off()

# displays in a pdf the aggregated SPADE tree for some samples
pdf("./04_SPADEVizR/SPADEVizR-figures/treeViewer-condition-BP.pdf", width = 15, height = 15)
treeViewer(NK_results, highlight = resultsDAC_BPvsPP_NK)
dev.off()


# Condition PP vs PB
pdf("./04_SPADEVizR/SPADEVizR-figures/treeViewer-condition-PP.pdf", width = 15, height = 15)
treeViewer(spade_results, samples = condition_PP)
dev.off()

# Condition BP vs PB
pdf("./04_SPADEVizR/SPADEVizR-figures/treeViewer-condition-PB.pdf", width = 15, height = 15)
treeViewer(spade_results, samples = condition_PB)
dev.off()

# displays in a pdf the aggregated SPADE tree for some samples, overlayed by the expression of HLADR
pdf("./05_SPADEVizR-figures/treeViewer-PBsamples-HLADR.pdf", width = 15, height = 15)
samples <- c("PBD28_BB078", "PBD28_BB231", "PBD28_BC641", "PBD28_BD619", "PBD28_BD620")
treeViewer(spade_results, samples = samples, marker = "HLADR")
dev.off()

# Treeviewer for NK clusters with DAC for BP vs PP condition
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/treeViewer_k", k_value, "_BPvsPP_NK.pdf"), width = 15, height = 15)
treeViewer(NK_results, highlight = resultsDAC_BPvsPP_NK)
dev.off()



## 4.5 CC Identification ====================================================
variable <- c("PPD00_BB078" = 50, "PPD00_BB231" = 50, "PPD00_BC641" = 50, "PPD00_BD619" = 50, "PPD00_BD620" = 50, "PBD08_BB078" = 32541, "PBD08_BB231" = 16769, "PBD08_BC641" = 16987, "PBD08_BD619" = 11592, "PBD08_BD620" = 7419, "PBD28_BB078" = 14621, "PBD28_BB231" = 7030, "PBD28_BC641" = 1048, "PBD28_BD619" = 3369, "PBD28_BD620" = 3881)
resultsCC <- identifyCC(spade_results, variable = variable, th.correlation = 0.8, th.pvalue = 0.05)

# diplays a char representations of the CC results
plot(resultsCC)

# exports the results in a text file
export(resultsCC, "./SPADEVizR-export/resultsCC.txt")




## 4.6 Identify Abundance Clusters =========================================
# Condition BP
resultsAC <- identifyAC(spade_results,
  samples = condition_BP,
  mu = 1,
  th.pvalue = 0.01
)
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/AbundanceClusters_k", k_value, "_BP.pdf"), height = 10, width = 17)
abundantClustersViewer(resultsAC)
dev.off()

# Condition PP
resultsAC <- identifyAC(spade_results,
  samples = condition_PP,
  mu = 1,
  th.pvalue = 0.01
)
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/AbundanceClusters_k", k_value, "_PP.pdf"), height = 10, width = 17)
abundantClustersViewer(resultsAC)
dev.off()

# Condition PB
resultsAC <- identifyAC(spade_results,
  samples = condition_PB,
  mu = 1,
  th.pvalue = 0.01
)
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/AbundanceClusters_k", k_value, "_PB.pdf"), height = 10, width = 17)
abundantClustersViewer(resultsAC)
dev.off()




## 4.6 MDS Representations =================================================
# MDS representation at the sample level with all clusters

# MDS BP vs PP
clusters.matrix <- resultsDAC_BPvsPP_NK@results
selected.clusters <- clusters.matrix[clusters.matrix$significant == TRUE, ]$cluster
print(selected.clusters)
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/MDS-Samples_k", k_value, "_BPvsPP_NK-clusters.pdf"), height = 10, width = 17)
MDSViewer(NK_results, space = "samples", clusters = selected.clusters, samples = c(condition_BP, condition_PP))
dev.off()

# MDS PP vs PB
clusters.matrix <- resultsDAC_PPvsPB_NK@results
selected.clusters <- clusters.matrix[clusters.matrix$significant == TRUE, ]$cluster
print(selected.clusters)
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/MDS-Samples_k", k_value, "_PPvsPB_NK-clusters.pdf"), height = 10, width = 17)
MDSViewer(NK_results, space = "samples", clusters = selected.clusters, samples = c(condition_PP, condition_PB))
dev.off()

# MDS BP vs PB
clusters.matrix <- resultsDAC_BPvsPB_NK@results
selected.clusters <- clusters.matrix[clusters.matrix$significant == TRUE, ]$cluster
print(selected.clusters)
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/MDS-Samples_k", k_value, "_BPvsPB_NK-clusters.pdf"), height = 10, width = 17)
MDSViewer(NK_results, space = "samples", clusters = selected.clusters, samples = c(condition_BP, condition_PB))
dev.off()

# MDS representation at the cluster level
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/MDS-Clusters_k", k_value, "NK-clusters.pdf"), height = 10, width = 17)
MDSViewer(NK_results)
dev.off()




## 4.7 Distogram ===========================================================
cluster = c("27")
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/Distogram_cluster_", cluster, ".pdf"))
distogramViewer(spade_results, samples = condition_PP, clusters = cluster)
dev.off()

cluster = M.III
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/Distogram_cluster_", cluster, ".pdf"))
distogramViewer(cluster_results, samples = condition_BP, clusters = M.III)
dev.off()


## 4.8 Kinetics visualization for NK metaclusters =============================

pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/kinetics_k", k_value, "_metaclusters_III.pdf"), height = 10, width = 10)
kineticsViewer(cluster_results, clusters = c("M.III"))
dev.off()

pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/kinetics_k", k_value, "_metaclusters_IV.pdf"), height = 10, width = 10)
kineticsViewer(cluster_results, clusters = c("M.IV"))
dev.off()

pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/kinetics_k", k_value, "_metaclusters_V.pdf"), height = 10, width = 10)
kineticsViewer(cluster_results, clusters = c("M.V"))
dev.off()

pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/kinetics_k", k_value, "_metaclusters_VI.pdf"), height = 10, width = 10)
kineticsViewer(cluster_results, clusters = c("M.VI"))
dev.off()



## 4.9 PhenoViewer -------------------------------------------------------------

pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/PhenoViewer_k", k_value, "_NK-metaclusters_III.pdf"), height = 8, width = 10)
phenoViewer(spade_results, clusters = M.III)
dev.off()

pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/PhenoViewer_k", k_value, "_NK-metaclusters_IV.pdf"), height = 8, width = 10)
phenoViewer(spade_results, clusters = M.IV)
dev.off()

pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/PhenoViewer_k", k_value, "_NK-metaclusters_V.pdf"), height = 8, width = 10)
phenoViewer(spade_results, clusters = M.V)
dev.off()

pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/PhenoViewer_k", k_value, "_NK-metaclusters_VI.pdf"), height = 8, width = 10)
phenoViewer(spade_results, clusters = M.VI)
dev.off()




# BROUILLON --------------------

### PHENOVIEWER REPRESENTATIONS
# select the significant clusters and print them
clusters.matrix <- resultsCC@spade_results
selected.clusters <- clusters.matrix[clusters.matrix$significant == TRUE, ]$cluster
print(clusters)

# displays a parallel coordinates representation for the CC
phenoViewer(spade_results, cluster = selected.clusters[1])




# GENERATE MASTER REPORT -------------------------------------------------------
# generates a report with the main SPADEVizR functionalities
createReport(spade_results, PDFfile = "./04_SPADEVizR/SPADEVizR-report-DAC_BPvsPP_metaclusters.pdf", select.plots = c("tree", "count", "heatmap", "kinetics_pheno", "distogram"), verbose = TRUE)

# Generate a report with the main SPADEVizR functionalities with DAC and CC reports for BP vs PP condition for NK clusters
createReport(NK_results, PDFfile = "./04_SPADEVizR/SPADEVizR-report-DAC_BPvsPP_NK_clusters.pdf", select.plots = c("count", "heatmap", "kinetics_pheno", "distogram", resultsDAC_BPvsPP_NK), verbose = TRUE)

# Generate a report with the main SPADEVizR functionalities with DAC and CC reports for BP vs PP condition for metaclusters
createReport(cluster_results, PDFfile = "./04_SPADEVizR/SPADEVizR-report-DAC_BPvsPP_metaclusters.pdf", select.plots = c("count", "kinetics_pheno", "distogram", resultsDAC_BPvsPP), verbose = TRUE)
