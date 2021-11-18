library(SingleCellExperiment)
library(SummarizedExperiment)
library(Seurat)
library(magrittr)
library(ensemblemerge)
library(scmetrics)

set.seed(1)
individual_methods = TRUE
j = 0 #track dataset #
scores_full = data.frame(score=double(), score_method=character(), method=character(), Datasets=character())
for(i in c(1:11, 13:19)){
	j = j + 1
	system(sprintf("wget -x -c -nH https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_%i.rds", i))
	data = readRDS(sprintf("skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=20211008a/dataset_%i.rds", i))
	data = as(data, "SingleCellExperiment")
	keep_feature <- rowSums(counts(data) > 0) > 0
	data <- data[keep_feature, ]
	data <- as.Seurat(data, data = NULL)
	if (ncol(data) > 20000){
		batches = data@meta.data$batch %>% unique() %>% length()
		Idents(data) = "batch"
		data <- subset(data, downsample = 10000/batches)
	}
	print(sprintf("i is %i",j))
	stopifnot("batch" %in% colnames(data@meta.data))
	stopifnot("cell_type" %in% colnames(data@meta.data))
	stopifnot(!is.na(data@meta.data$cell_type))
	stopifnot(!is.na(data@meta.data$batch))
		
	#running Scanorama results
	params = setParams(method = "fastMNN", return = "Seurat")
	fastMNN = Merge(params, data)
	fastMNN <- RunUMAP(fastMNN, reduction = "mnn", dims = 1:20, seed.use = 1)
	fastMNN <- FindNeighbors(fastMNN, reduction = 'umap', dims = 1:2, verbose = FALSE)
	fastMNN <- FindClusters(fastMNN, algorithm = 3, resolution = 0.01, verbose = FALSE)
	fastMNN$cluster <- Idents(fastMNN)
	#convert
	fastMNN_graph = fastMNN[[Graphs(fastMNN)[grep("_nn", Graphs(fastMNN))]]]
	fastMNN = as.SingleCellExperiment(fastMNN)
	metadata(fastMNN)[["graph"]] = fastMNN_graph
	
	#SCmetrics testing pipeline
	config <- list(
	  label_field = "cell_type",
	  batch_field = "batch",
	  n_neighbors = 15L,
	  resolution_min = 1e-3,
	  resolution_max = 1,
	  resolution_length = 40L,
	  metrics_max = c('ARI')
	) %>%
	  expand.grid(stringsAsFactors = FALSE)

	ARGS<- config[1, ]

	params <- new('ClusterLouvain', n_neighbors = ARGS[['n_neighbors']], resolution_min = ARGS[['resolution_min']], resolution_max = ARGS[['resolution_max']], resolution_length = ARGS[['resolution_length']])

	r <- 'UMAP'

	names = "fastMNN"

	nmi_list = list()
	ari_cell_list = list()
	ari_batch_list = list()
	asw_cell_list = list()
	asw_batch_list = list()
	graph_con_list = list()
	
	tests = fastMNN

	tests <- suppressWarnings(cluster_cells(tests, params, reduction = r, label=  ARGS[['label_field']], cluster = 'cluster', metric = ARGS[['metrics_max']]))
	nmi <- metrics(tests, new('MetricsNMI'), label = ARGS[['label_field']], cluster = 'cluster')
	nmi_list = append(nmi_list, nmi)
	ari_cell <- metrics(tests, new('MetricsARI'), label = ARGS[['label_field']], cluster = 'cluster')
	ari_cell_list = append(ari_cell_list, ari_cell)
	ari_batch <- metrics(tests, new('MetricsARI'), label = ARGS[['batch_field']], cluster = 'cluster')
	ari_batch <- 1 - ari_batch
	ari_batch_list = append(ari_batch_list, ari_batch)
	asw_cell <- metrics(tests, new('MetricsASW'), label = ARGS[['label_field']], reduction = r)
	asw_cell <- (asw_cell + 1) / 2
	asw_cell_list = append(asw_cell_list, asw_cell)
	asw_batch <- metrics(tests, new('MetricsASW'), label = ARGS[['batch_field']], reduction = r)
	asw_batch <- 1 - abs(asw_batch)
	asw_batch_list = append(asw_batch_list, asw_batch)
	graph_con <- metrics(tests, new('MetricsGraphConnectivity'), label = ARGS[['batch_field']], graph = "graph")
	graph_con_list = append(graph_con_list, graph_con)
	sprintf('%s Score:', names[i]) %>% message()
	sprintf('NMI: %.3f', nmi) %>% message()
	sprintf('ARI(cell): %.3f', ari_cell) %>% message()
	sprintf('ASW(cell): %.3f', asw_cell) %>% message()
	sprintf('ARI(batch): %.3f', ari_batch) %>% message()
	sprintf('ASW(batch): %.3f', asw_batch) %>% message()
	sprintf('GC: %.3f', graph_con) %>% message()

	names(ari_cell_list) = names
	names(nmi_list) = names
	names(ari_cell_list) = names
	names(ari_batch_list) = names
	names(asw_cell_list) = names
	names(asw_batch_list) = names
	names(graph_con_list) = names

	nmi_list = as.data.frame(nmi_list)
	ari_cell_list = as.data.frame(ari_cell_list)
	ari_batch_list = as.data.frame(ari_batch_list)
	asw_cell_list = as.data.frame(asw_cell_list)
	asw_batch_list = as.data.frame(asw_batch_list)
	graph_con_list = as.data.frame(graph_con_list)

	nmi_list = t(nmi_list)
	ari_cell_list = t(ari_cell_list)
	ari_batch_list = t(ari_batch_list)
	asw_cell_list = t(asw_cell_list)
	asw_batch_list = t(asw_batch_list)
	graph_con_list = t(graph_con_list)

	colnames(ari_cell_list) = c("score")
	colnames(nmi_list) = c("score")
	colnames(ari_batch_list) = c("score")
	colnames(asw_cell_list) = c("score")
	colnames(asw_batch_list) = c("score")
	colnames(graph_con_list) = c("score")

	nmi_list = as.data.frame(nmi_list)
	ari_cell_list = as.data.frame(ari_cell_list)
	ari_batch_list = as.data.frame(ari_batch_list)
	asw_cell_list = as.data.frame(asw_cell_list)
	asw_batch_list = as.data.frame(asw_batch_list)
	graph_con_list = as.data.frame(graph_con_list)

	nmi_list$score_method = "NMI"
	ari_cell_list$score_method = "ARI_Cell"
	ari_batch_list$score_method = "ARI_Batch"
	asw_cell_list$score_method = "ASW_Cell"
	asw_batch_list$score_method = "ASW_Batch"
	graph_con_list$score_method = "GC"

	nmi_list$method = rownames(nmi_list)
	ari_cell_list$method = rownames(ari_cell_list)
	ari_batch_list$method = rownames(ari_batch_list)
	asw_cell_list$method = rownames(asw_cell_list)
	asw_batch_list$method = rownames(asw_batch_list)
	graph_con_list$method = rownames(graph_con_list)

	nmi_list = nmi_list[order(-nmi_list$score),]
	ari_cell_list = ari_cell_list[order(-ari_cell_list$score),]
	ari_batch_list = ari_batch_list[order(-ari_batch_list$score),]
	asw_cell_list = asw_cell_list[order(-asw_cell_list$score),]
	asw_batch_list = asw_batch_list[order(-asw_batch_list$score),]
	graph_con_list = graph_con_list[order(-graph_con_list$score),]

	scores = list(nmi_list, ari_cell_list, ari_batch_list, asw_cell_list, asw_batch_list, graph_con_list)
	scores = Reduce(rbind, scores)
	scores$Datasets = sprintf("Dataset_%i", j)
	write.csv(scores, file = sprintf("Dataset_%i_Metric_Scores.csv", j))
	rownames(scores_full) = NULL
	scores_full = rbind(scores_full, scores)
}
write.csv(scores_full, file = "fastMNN.csv")