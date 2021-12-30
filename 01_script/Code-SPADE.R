# 1. Setup -------------------------------------------------------------------
# Make sure to install librarian before launching the script.
# install.packages("librarian")
# Make sure to have opened with an R project file first.
librarian::shelf(igraph, spade)
set.seed(123) # for reproducibility


# 2. Load data and select clustering markers ---------------------------------
clustering_markers = read.delim('./02_data/NK_panel.txt')
clustering_markers = clustering_markers[-c(1:7, 18, 24, 29, 39, 43:45), 1]
files = list.files("./02_data/NK_FCS/",
                   pattern = "*.fcs",
                   full.names = TRUE)

# Select parameters
DS = 0.05 #the downsampling parameter
K = 50 # the number of cell clusters to identify
output_dir = paste0("./03_spade_analysis/spade_k", K, "/")

# 3. Launch SPADE.driver -----------------------------------------------------
SPADE.driver(
  files,
  out_dir = output_dir,
  cluster_cols = clustering_markers,
  k = K,
  downsampling_target_pctile = DS
)
# 4. Save results in PDF format ----------------------------------------------
LAYOUT_TABLE = read.table(paste0(output_dir, "layout.table"))
mst_graph = igraph:::read.graph(paste(output_dir,
                                      "mst.gml",
                                      sep = .Platform$file.sep),
                                      format = "gml")
SPADE.plot.trees(
  mst_graph,
  output_dir,
  file_pattern = "*fcs*Rsave",
  layout = as.matrix(LAYOUT_TABLE),
  out_dir = paste0(output_dir, "/PDF/"),
  size_scale_factor = 1
)
