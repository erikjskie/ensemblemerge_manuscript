library(SingleCellExperiment)
library(SummarizedExperiment)
library(Seurat)
library(magrittr)
library(ensemblemerge)
library(scmetrics)
library(dplyr) #dataframe management

batch = "batch"
cell_type = "cell_type"
method = "Uncorrected"
n_iterations = 1000	n_neighbors = 15L #neighbors for clustering analysis
resolution_min = 1e-3 #resolution min change from one step to the next
resolution_max = 1 #resolution max change from one step to the next
resolution_length = 40L #number of steps in resolving resolution
metrics_max = c('ARI') #metric used to determine performance
downsample = 0.25

scores = data.frame(score_method=character(), score=double(),  Datasets=character(), seed = integer(), method=character()) #set dataframe with correct dimensions for analysis output

for(i in 1:18){
	system(sprintf("wget -x -c -nH https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_%i.rds", i)) #download dataset from aws
	data = readRDS(sprintf("skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_%i.rds", i)) # read dataset into memory as SummarizedExperiment object
	data = as(data, "SingleCellExperiment") #convert to SingleCellExperiment
	keep_feature <- rowSums(counts(data) > 0) > 0 #Identify genes of relevance
	data <- data[keep_feature, ] #trim extra genes left over from merging with other datasets
	data <- as.Seurat(data, data = NULL) #convert to Seurat object
	batches = data@meta.data[[batch]] %>% unique() %>% length() # get number of batches in dataset for downsampling
	Idents(data) = batch #set identity to batch
	data <- subset(data, downsample = floor((ncol(data)/batches)*downsample)) # reduce number of cells to specified downsampling
	stopifnot(batch %in% colnames(data@meta.data)) # stop if batch isn't present
	stopifnot(cell_type %in% colnames(data@meta.data)) # stop if cell type isn't present
	stopifnot(!is.na(data@meta.data$cell_type)) # stop if cell_type has missing entries
	stopifnot(!is.na(data@meta.data$batch)) # stop if batch has missing entries

	#running uncorrected results
	params = setParams(method = method, return = "Seurat") #set params object for ensemblemerge
	Uncorrected = Merge(params, data) #run uncorrected integration of batches with ensemblemerge
	for(seed in 1:n_iterations){
		uncorrected <- Seurat::RunUMAP(Uncorrected, reduction = "pca", dims = 1:20, seed.use = seed)
		uncorrected <- FindNeighbors(uncorrected, reduction = 'umap', dims = 1:2, verbose = FALSE)
		uncorrected <- FindClusters(uncorrected, algorithm = 3, resolution = 0.01, verbose = FALSE)
		uncorrected$cluster <- Idents(uncorrected)
		#convert
		uncorrected_graph = uncorrected[[Graphs(uncorrected)[grep("_nn", Graphs(uncorrected))]]] #assign graph to variable
		uncorrected = as.SingleCellExperiment(uncorrected) # convert to SingleCellExperiment object
		metadata(uncorrected)[["graph"]] = uncorrected_graph # assign graph back to SingleCellExperiment Object
		#SCmetrics testing pipeline

		params <- new('ClusterLouvain', 
									n_neighbors = n_neighbors, 
									resolution_min = resolution_min, 
									resolution_max = resolution_max, 
									resolution_length = resolution_length)

		r <- 'UMAP'

		names = method

		tests <- suppressWarnings(cluster_cells(uncorrected, params, reduction = r, label=  cell_type, cluster = 'cluster', metric = metrics_max))
		metrics = data.frame(NMI = metrics(tests, new('MetricsNMI'), label = cell_type, cluster = 'cluster'),
					ARI_Cell = metrics(tests, new('MetricsARI'), label = cell_type, cluster = 'cluster'),
					ARI_Batch = 1 - metrics(tests, new('MetricsARI'), label = batch, cluster = 'cluster'), # 1-metric for batch metrics
					ASW_Cell = (metrics(tests, new('MetricsASW'), label = cell_type, reduction = r)+1)/2, #rescaling asw 
					ASW_Batch = 1-abs(metrics(tests, new('MetricsASW'), label = batch, reduction = r)), # 1-metric for batch and abs() for rescaling
					GC = metrics(tests, new('MetricsGraphConnectivity'), label = batch, graph = "graph")) %>%
					t() %>% as.data.frame() %>% #transpose
					rownames_to_column(var = "score_method") #set rownames to data in table
		colnames(metrics)[2] = "score"
		sprintf('%s Score:', names[j]) %>% message()
		sprintf('NMI: %.3f', metrics[which(metrics$score_method == 'NMI'),"score"]) %>% message()
		sprintf('ARI(cell): %.3f', metrics[which(metrics$score_method == 'ARI_Cell'),"score"]) %>% message()
		sprintf('ASW(cell): %.3f', metrics[which(metrics$score_method == 'ASW_Cell'),"score"]) %>% message()
		sprintf('ARI(batch): %.3f', metrics[which(metrics$score_method == 'ARI_Batch'),"score"]) %>% message()
		sprintf('ASW(batch): %.3f', metrics[which(metrics$score_method == 'ASW_Batch'),"score"]) %>% message()
		sprintf('GC: %.3f', metrics[which(metrics$score_method == 'GC'),"score"]) %>% message()

		metrics$Datasets = sprintf("Dataset_%i", j)
		metrics$seed = seed
		metrics$method = method
		scores = rbind(scores, metrics)
	}
}

write.csv(scores, file = "fig9.csv")
