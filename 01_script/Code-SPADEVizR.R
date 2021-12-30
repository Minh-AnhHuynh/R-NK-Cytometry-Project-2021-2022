# 1. Setup -------------------------------------------------------------------
# Make sure to install librarian before launching the script.
# Make sure to have opened with an R project file first.
librarian::shelf(tchitchek-lab/SPADEVizR, data.table, ggplot2, stringr, dplyr, flowCore)
set.seed(123) # for reproducibility
seed_number = 123

# 2. Data loading ------------------------------------------------------------
clustering_markers = read.delim('./02_data/NK_panel.txt')
# Filtering markers
excluded_markers = clustering_markers[c(1:7,18,24,29,39,43,44,45),1]
excluded_markers2 = c("File Number", "density", "cells-(Ir191)Di", "cells-(Ir193)Di")
excluded_markers = c(excluded_markers, excluded_markers2)

## Select k_value
k_value = 100 

## Load SPADE results
output_dir = paste0("./03_spade_analysis/spade_k", k_value)
spade_results = importResultsFromSPADE(output_dir,
                                       exclude.markers = excluded_markers)
# To check if spade_results is correct :
# head(spade_results@cluster.phenotypes)


# 3. Quality Control ---------------------------------------------------------

# Generate a QC report to detect cell cluster phenotype with uniform marker distribution expressions 
accuracyQC <- qcUniformClusters(spade_results,
                                density.PDFfile = paste0("./04_SPADEVizR/QualityControls/UniformClusters_density_k",k_value,"_seed", seed_number,".pdf"),
                                heatmap.PDFfile = paste0("./04_SPADEVizR/QualityControls/UniformClusters_density_k",k_value,"_seed", seed_number,".pdf"),
                                uniform.test   = "unimodality")

# Display the percentage of accuracy clusters
print(accuracyQC$perc)

# Save the accuracy matrix in a text file
write.table(accuracyQC$accuracy.matrix, paste0("./04_SPADEVizR/QualityControls/UniformClusters-accuracy-matrix_k_", k_value,".txt"),quote=FALSE,sep="\t",col.names=NA)

# Generate a QC report to detect cell cluster phenotypes with low number of cell clusters
smallclusterQC <- qcSmallClusters(spade_results, PDFfile = paste0("./04_SPADEVizR/QualityControls/QCreport-SmallClusters_heatmap_k",k_value,".pdf"), th = 500)

# Display the percentage of cell cluster phenotypes with low number of cell clusters 
print(smallclusterQC$perc)

# Save the small cluster matrix in a text file
write.table(smallclusterQC$small.clusters, paste0("./04_SPADEVizR/QualityControls/UniformClusters-smallclusters-matrix_k", k_value,".txt"),quote=FALSE,sep="\t",col.names=NA)
# Notes : Cluster 57 and 73 are small clusters

# 4. Analysis and Figures ----------------------------------------------------

## 4.1 Heatmap generation ====================================================

pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/heatmap_spade_k",k_value,".pdf"), height=10,width=17)
heatmapViewer(spade_results)
dev.off()

## 4.2 Identification of Differentially Abundant Clusters ==========================================================================

### 4.2.1 Data frame creation : row = cell ; col = markers ##########################################################################
# This step is done in order to get the assignment data.frame without selecting samples manually

#Load data using flowSet
files = list.files("./02_data/NK_FCS/", full.names = TRUE)
FlowSet = read.flowSet(files)

# Select markers of interest
clustering_markers = read.delim('./02_data/NK_panel.txt')
markers2norm = clustering_markers[-c(1:7,18,24,29,39,43:45),1]

# Arcsinh transformation of markers to be able to use Flowset
ArcTrans = arcsinhTransform()
TransList = transformList(from = markers2norm, tfun = ArcTrans)
FlowSet = TransList %on% FlowSet
# Generation of files
data = data.frame()
for (i in 1:length(FlowSet)) {
  print(paste0("fichier numero ", i))
  data_sub <- exprs(FlowSet[[i]])
  data_sub <- data.frame(data_sub)
  data_sub$sample <- basename(files)[i]
  data <- rbind(data, data_sub)
  
}

### 4.2.2 Assigning sample biological conditions (metadata information) ####################
assignments = data.frame(bc = str_sub(data$sample,1,8),
                         tp = str_sub(data$sample,1,8),
                         #hour = str_sub(data$sample,6,8),
                         ind = str_sub(data$sample,10,14),
                         sample_names = str_sub(data$sample,1,23))

# Filtering samples to obtain only 42 values corresponding to 42 samples
assignments = assignments %>% group_by(sample_names) %>% slice(1)

# Put sample names in row names
sample = assignments$sample_names
assignments[4]=NULL
row.names(assignments) = sample

# Put assignment data.frame tidied into spade_results
spade_results = assignContext(spade_results, assignments = assignments)


### 4.2.3 Select monkey samples according to their conditions respectively : ##########################################################################

# Before Prime
condition_BP = c("BPD19H00_BB078_CD3-CD8+", "BPD19H00_BB231_CD3-CD8+", "BPD19H00_BC641_CD3-CD8+", "BPD19H00_BD620_CD3-CD8+")

# Post Prime
condition_PP = c("PPD00H00_BB078_CD3-CD8+", "PPD00H00_BB231_CD3-CD8+", "PPD00H00_BD620_CD3-CD8+", "PPD00H03_BB078_CD3-CD8+", "PPD00H03_BB231_CD3-CD8+", "PPD00H03_BC641_CD3-CD8+", "PPD00H03_BD620_CD3-CD8+", "PPD00H06_BB078_CD3-CD8+", "PPD00H06_BB231_CD3-CD8+", "PPD00H06_BC641_CD3-CD8+", "PPD01H00_BB078_CD3-CD8+", "PPD01H00_BB231_CD3-CD8+", "PPD01H00_BC641_CD3-CD8+", "PPD01H00_BD620_CD3-CD8+", "PPD03H00_BB078_CD3-CD8+", "PPD03H00_BB231_CD3-CD8+", "PPD03H00_BC641_CD3-CD8+", "PPD14H00_BB078_CD3-CD8+", "PPD14H00_BB231_CD3-CD8+", "PPD14H00_BC641_CD3-CD8+", "PPD14H00_BD620_CD3-CD8+")

# Post Boost
condition_PB = c("PBD00H00_BB078_CD3-CD8+", "PBD00H00_BB231_CD3-CD8+", "PBD00H00_BD620_CD3-CD8+", "PBD00H03_BB078_CD3-CD8+", "PBD00H03_BB231_CD3-CD8+", "PBD00H03_BC641_CD3-CD8+", "PBD00H06_BB078_CD3-CD8+", "PBD00H06_BB231_CD3-CD8+", "PBD00H06_BC641_CD3-CD8+", "PBD00H06_BD620_CD3-CD8+", "PBD01H00_BB078_CD3-CD8+", "PBD01H00_BB231_CD3-CD8+", "PBD01H00_BC641_CD3-CD8+", "PBD01H00_BD620_CD3-CD8+", "PBD03H00_BB078_CD3-CD8+", "PBD03H00_BB231_CD3-CD8+", "PBD03H00_BD620_CD3-CD8+")


#condition_BP = spade_results@assignments[grepl("^BP", rownames(spade_results@assignments)),1]

#rownames(spade_results@assignments) = row.names(assignments)

#condition_BP = rownames(spade_results@assignments[spade_results@assignments$bc == "BPD19H00",])

#condition_BP = rownames(spade_results@assignments[spade_results@assignments$tp == "BPD19H00", ])

# selects macaque samples post-prime vaccination
#condition_PP = assignments$sample_names[startsWith(assignments$sample_names, "PP")]

# selects macaque samples post-boost vaccination
#condition_PB = assignments$sample_names[startsWith(assignments$sample_names, "PB")]



## 4.3 Volcano plot  =====================================================
# Volcano plot Condition BP vs PP 
resultsDAC_BPvsPP <- identifyDAC(
  spade_results,
  condition1 = condition_BP,
  condition2 = condition_PP,
  th.pvalue = 0.05,
  th.fc = 2,
  method.paired = FALSE
)
export(resultsDAC_BPvsPP, filename = paste0("./04_SPADEVizR/ResultsDAC/resultsDAC_k", k_value,"_BPvsPP.txt"))
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/VolcanoDAC_k",k_value,"_BPvsPP.pdf"), height=10,width=17)
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
export(resultsDAC_PPvsPB, filename = paste0("./04_SPADEVizR/ResultsDAC/resultsDAC_k", k_value,"_PPvsPB.txt"))
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/VolcanoDAC_k",k_value,"_PPvsPB.pdf"), height=10,width=17)
volcanoViewer(resultsDAC_PPvsPB)
dev.off()

# Volcano plot Condition BP vs PB
resultsDAC_BPvsPB<- identifyDAC(
  spade_results,
  condition1 = condition_BP,
  condition2 = condition_PB,
  th.pvalue = 0.05,
  th.fc = 2,
  method.paired = FALSE
)
export(resultsDAC_BPvsPB, filename = paste0("./04_SPADEVizR/ResultsDAC/resultsDAC_k", k_value,"_BPvsPB.txt"))
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/VolcanoDAC_k",k_value,"_BPvsPB.pdf"), height=10,width=17)
volcanoViewer(resultsDAC_BPvsPB)
dev.off()


## 4.4 Tree Viewer ============================================================

# displays in a pdf the SPADE tree for sample PBD28_BB078
#pdf("./05_figures/treeViewer-PBD28_BB078.pdf", width=15, height=15)
#treeViewer(spade_results,samples="PBD28_BB078")
#dev.off()

# displays in a pdf the aggregated SPADE tree for all samples
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/treeViewer_k", k_value,".pdf"), width=15, height=15)
treeViewer(spade_results, highlight = resultsDAC)
dev.off()

# displays in a pdf the aggregated SPADE tree for some samples
pdf("./04_SPADEVizR/SPADEVizR-figures/treeViewer-condition-BPvsPP.pdf", width=15, height=15)
treeViewer(spade_results, samples = condition_BP, highlight = resultsDAC_BPvsPP)
dev.off()

# displays in a pdf the aggregated SPADE tree for some samples, overlayed by the expression of HLADR 
pdf("./05_SPADEVizR-figures/treeViewer-PBsamples-HLADR.pdf", width=15, height=15)
samples <- c("PBD28_BB078", "PBD28_BB231", "PBD28_BC641", "PBD28_BD619", "PBD28_BD620")
treeViewer(spade_results, samples = samples, marker = "HLADR")
dev.off()


## 4.5 CC IDENTIFICATION ========================================================
variable <- c("PPD00_BB078" = 50, "PPD00_BB231" = 50, "PPD00_BC641" = 50, "PPD00_BD619" = 50, "PPD00_BD620" = 50, "PBD08_BB078" = 32541, "PBD08_BB231" = 16769, "PBD08_BC641" = 16987, "PBD08_BD619" = 11592, "PBD08_BD620" = 7419, "PBD28_BB078" = 14621, "PBD28_BB231" = 7030, "PBD28_BC641" = 1048, "PBD28_BD619" = 3369, "PBD28_BD620" = 3881)
resultsCC <- identifyCC(spade_results, variable = variable, th.correlation = 0.8, th.pvalue = 0.05)

# diplays a char representations of the CC results
plot(resultsCC)

# exports the results in a text file
export(resultsCC,"./SPADEVizR-export/resultsCC.txt")


## 4.6 Identify Abundance Clusters =============================================

resultsAC <- identifyAC(spade_results,
           samples = condition_BP,
           mu = 1,
           th.pvalue = 0.05)
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/AbundanceClusters_k",k_value,"_BP.pdf"), height=10,width=17)
abundantClustersViewer(resultsAC)
dev.off()

## 4.6 MDS REPRESENTATIONS ============================================
# MDS representation at the sample level with all clusters

# Transformation of assignments from a tibble into a data.frame to make MDSViewer work
spade_results@assignments <- as.data.frame(spade_results@assignments)

clusters.matrix   <- resultsDAC_BPvsPP@results
selected.clusters <- clusters.matrix[clusters.matrix$significant==TRUE,]$cluster
print(selected.clusters)
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/MDS-Samples_k",k_value,"_BPvsPP.pdf"), height=10,width=17)
MDSViewer(spade_results, space = "samples", clusters = selected.clusters)
dev.off()

# MDS PP vs PB
spade_results@assignments <- as.data.frame(spade_results@assignments)

clusters.matrix   <- resultsDAC_PPvsPB@results
selected.clusters <- clusters.matrix[clusters.matrix$significant==TRUE,]$cluster
print(selected.clusters)
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/MDS-Samples_k",k_value,"_PPvsPB.pdf"), height=10,width=17)
MDSViewer(spade_results, space = "samples", clusters = selected.clusters)
dev.off()

# MDS BP vs PB
clusters.matrix   <- resultsDAC_BPvsPB@results
selected.clusters <- clusters.matrix[clusters.matrix$significant==TRUE,]$cluster
print(selected.clusters)
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/MDS-Samples_k",k_value,"_BPvsPB.pdf"), height=10,width=17)
MDSViewer(spade_results, space = "samples", clusters = selected.clusters)
dev.off()

selected.clusters <- rownames(clusters.matrix)
print(selected.clusters)
MDSViewer(spade_results, space = "samples", clusters = selected.clusters)


# MDS representation at the cluster level
pdf(paste0("./04_SPADEVizR/SPADEVizR-figures/MDS-Clusters_k",k_value,"_BPvsPP.pdf"), height=10,width=17)
MDSViewer(spade_results)
dev.off()

## 4.7 Distogram ==============================================================
# Cluster 21 is the only AC cluster in identifyAC() for BP condition
sample = condition_BP
cluster = c("21")
distogramViewer(spade_results, samples = condition_BP, clusters = cluster)



### PHENOVIEWER REPRESENTATIONS
# select the significant clusters and print them
clusters.matrix   <- resultsCC@spade_results
selected.clusters <- clusters.matrix[clusters.matrix$significant==TRUE,]$cluster
print(clusters)

# displays a parallel coordinates representation for the CC
phenoViewer(spade_results,cluster=selected.clusters[1])
######

### GENERATE MASTER REPORT
# generates a report with the main SPADEVizR functionalities
createReport(spade_results, PDFfile = "./SPADEVizR-masterReport/SPADEVizR-report.pdf", select.plots = c("tree","count", "heatmap", "kinetics_pheno", "distogram"), verbose = TRUE)

# generates a report with the main SPADEVizR functionalities with DAC and CC reports
createReport(spade_results, PDFfile = "./SPADEVizR-masterReport/SPADEVizR-report-withDACandCC.pdf", select.plots = c("tree","count", "heatmap", "kinetics_pheno", "distogram", resultsDAC, resultsCC), verbose = TRUE)
######								 

### CLUSTER MANIPULATION
# merges the abundances and the phenotypes of clusters 1 and 2 into a new cluster in a Results object
newresults <- mergeClusters(spade_results, clusters=c("1","2"),name="combined")
print(newresults@cluster.names)

# deletes the clusters 1 and 2 from a Results object
newresults <- removeClusters(spade_results, clusters=c("1","2"))
print(newresults@cluster.names)
######

### CLUSTER ANNOTATIONS
# defines an annotation dataframe
annotations <- data.frame()
annotations["resting_memory","CD21"] <- "c(2,3,4,5)"
annotations["resting_memory","CD27"] <- "c(2,3,4,5)"
annotations["activated_memory","CD21"] <- "1"
annotations["activated_memory","CD27"] <- "c(2,3,4,5)"
annotations["naive_B","CD21"] <- "c(2,3,4,5)"
annotations["naive_B","CD27"] <- "1"
annotations["tissuelike_memory","CD21"] <- "1"
annotations["tissuelike_memory","CD27"] <- "1"
print(annotations)

# annotates the cell clusters in a Results object#
# cell clusters are renamed according to the population names
spade_results <- annotateClusters(spade_results, annotations=annotations)
spade_results@cluster.names
###