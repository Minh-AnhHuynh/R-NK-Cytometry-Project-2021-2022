librarian::shelf(SPADEVizR, data.table, ggplot2, stringr, dplyr, flowCore)

### DATA LOADING

clustering_markers = read.delim('./02_data/NK_panel.txt')
excluded_markers = clustering_markers[c(1:7,18,24,29,39,43,44,45),1]


# loads the SPADE results contained in the folder "SPADE" on the Desktop

################
k_value = 100 ##
################

output_dir = paste0("./03_spade_analysis/spade_k", k_value)

spade_results = importResultsFromSPADE(output_dir,
                                       exclude.markers = excluded_markers)

# génération de la heatmap 
heatmapViewer(spade_results)


pdf(paste0("./05_figures/heatmap_spade_k_",k_value,".pdf"), height=10,width=17)
plot(spade_results)
dev.off()


### QC CONTROL
set.seed(123) #for reproducibility

# generates a QC report to detect cell cluster phenotypes with uniform marker distribution expressions 
accuracyQC <- qcUniformClusters(spade_results,
                                density.PDFfile = paste0("./04_quality_controls/UniformClusters_density_k_",k_value,"seed",".pdf"),
                                heatmap.PDFfile = paste0("./04_quality_controls/UniformClusters_heatmap_k_",k_value,"seed",".pdf"),
                                uniform.test   = "unimodality")


# displays the percentage of accuracy clusters
print(accuracyQC$perc)

# saves the accuracy matrix in a text file
write.table(accuracyQC$accuracy.matrix, paste0("./04_quality_controls/UniformClusters-accuracy-matrix_k_", k_value,".txt"),quote=FALSE,sep="\t",col.names=NA)

# generates a QC report to detect cell cluster phenotypes with low number of cell clusters 								  
smallclusterQC <- qcSmallClusters(spade_results, PDFfile = paste0("./04_quality_controls/QCreport-SmallClusters_heatmap_k_",k_value,".pdf"), th = 500)

# displays the percentage of cell cluster phenotypes with low number of cell clusters 
print(smallclusterQC$perc)

# saves the small cluster matrix in a text file
write.table(smallclusterQC$small.clusters, paste0("./04_quality_controls/UniformClusters-smallclusters-matrix_k_", k_value,".txt"),quote=FALSE,sep="\t",col.names=NA)
###








### DAC IDENTIFICATION 

# creation d'une dataframe avec ligne = cellule et colonne = marqueur
#data = data.frame()

for (i in 1:length(FlowSet)) {
  print(paste0("fichier numero ", i))
  data_sub <- exprs(FlowSet[[i]])
  data_sub <- data.frame(data_sub)
  data_sub$sample <- basename(files)[i]
  data <- rbind(data, data_sub)
  
}

####Création des différentes conditions pour chaque échantillon
assignments = data.frame(bc = str_sub(data$sample,1,8),
                         tp = str_sub(data$sample,1,8),
                         #hour = str_sub(data$sample,6,8),
                         ind = str_sub(data$sample,10,14),
                         sample_names = data$sample)

#Tri des samples pour obtenir 42 echantillons
assignments = assignments %>% group_by(sample_names) %>% slice(1) 
sample = assignments$sample_names
assignments[4]=NULL
row.names(assignments) = sample

# selects macaque samples before prime vaccination
condition_BP = assignments$sample_names[startsWith(assignments$sample_names, "BP")]

# selects macaque samples post-prime vaccination
condition_PP = assignments$sample_names[startsWith(assignments$sample_names, "PP")]

# selects macaque samples post-boost vaccination
condition_PB = assignments$sample_names[startsWith(assignments$sample_names, "PB")]

spade_results = assignContext(spade_results, assignments = assignments)

# identification des clusters diff?rentiellements exprim?s entre deux conditions
resultsDAC <- identifyDAC(spade_results, 
                          condition1 = condition_BP, condition2 = condition_PP, 
                          th.pvalue = 0.05, 
                          th.fc = 2, 
                          method.paired = FALSE)



# diplays a volcan plot representations of the DAC results
plot(resultsDAC)

# exports the results in a text file
export(resultsDAC,"./SPADEVizR-export/resultsDAC.txt")
######









### CC IDENTIFICATION 
variable <- c("PPD00_BB078" = 50, "PPD00_BB231" = 50, "PPD00_BC641" = 50, "PPD00_BD619" = 50, "PPD00_BD620" = 50, "PBD08_BB078" = 32541, "PBD08_BB231" = 16769, "PBD08_BC641" = 16987, "PBD08_BD619" = 11592, "PBD08_BD620" = 7419, "PBD28_BB078" = 14621, "PBD28_BB231" = 7030, "PBD28_BC641" = 1048, "PBD28_BD619" = 3369, "PBD28_BD620" = 3881)
resultsCC <- identifyCC(spade_results, variable = variable, th.correlation = 0.8, th.pvalue = 0.05)

# diplays a char representations of the CC results
plot(resultsCC)

# exports the results in a text file
export(resultsCC,"./SPADEVizR-export/resultsCC.txt")
######

heatmapViewer(spade_results)




### MDS REPRESENTATIONS
# displays a MDS representation at the sample level with all clusters
MDSViewer(spade_results, space = "samples")

# select the significant clusters and print them
clusters.matrix   <- resultsDAC@results
selected.clusters <- clusters.matrix[clusters.matrix$significant==TRUE,]$cluster
print(selected.clusters)

# displays a MDS representation at the sample level with DAC
MDSViewer(spade_results, space = "samples", clusters=selected.clusters)
######







### PHENOVIEWER REPRESENTATIONS
# select the significant clusters and print them
clusters.matrix   <- resultsCC@results
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