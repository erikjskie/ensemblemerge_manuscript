library(SingleCellExperiment)
library(SummarizedExperiment)
library(Seurat)
library(magrittr)
library(ensemblemerge)
library(scmetrics)

set.seed(1)
individual_methods = FALSE
scores_list = list()

for(h in 1:5){
	for(i in 1:18){
		system(sprintf("wget -x -c -nH https://s3.msi.umn.edu/skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=2021219a/dataset_%i.rds", i))
		data = readRDS(sprintf("skiex003/datasets/dataset=JC_benchmark_scRNAseq_version=2021219a/dataset_%i.rds", i))
		data = as(data, "SingleCellExperiment")
		keep_feature <- rowSums(counts(data) > 0) > 0
		data <- data[keep_feature, ]
		data <- as.Seurat(data, data = NULL)
		if(i == 3){
			data@meta.data$CellType = data@meta.data$Group ##correct for synthetic data
		}
		
		if(individual_methods){
			#running uncorrected results
			params = setParams(method = "Uncorrected", return = "Seurat")
			uncorrected = Merge(params, data)
			uncorrected <- RunUMAP(uncorrected, reduction = "pca", dims = 1:20, seed.use = 1)
			uncorrected <- FindNeighbors(uncorrected, reduction = 'umap', dims = 1:2, verbose = FALSE)
			uncorrected <- FindClusters(uncorrected, algorithm = 3, resolution = 0.01, verbose = FALSE)
			uncorrected$cluster <- Idents(uncorrected)

			#running Harmony results
			params = setParams(method = "Harmony", return = "Seurat")
			Harmony = Merge(params, data)
			Harmony <- RunUMAP(Harmony, reduction = "harmony", dims = 1:20, seed.use = 1)
			Harmony <- FindNeighbors(Harmony, reduction = 'umap', dims = 1:2, verbose = FALSE)
			Harmony <- FindClusters(Harmony, algorithm = 3, resolution = 0.01, verbose = FALSE)
			Harmony$cluster <- Idents(Harmony)

			#running Seurat results
			params = setParams(method = "Seurat", return = "Seurat")
			Seurat = Merge(params, data)
			Seurat <- RunUMAP(Seurat, reduction = "pca", dims = 1:20, seed.use = 1)
			Seurat <- FindNeighbors(Seurat, reduction = 'umap', dims = 1:2, verbose = FALSE)
			Seurat <- FindClusters(Seurat, algorithm = 3, resolution = 0.01, verbose = FALSE)
			Seurat$cluster <- Idents(Seurat)

			#running BBKNN results
			params = setParams(method = "BBKNN", return = "Seurat")
			BBKNN = Merge(params, data)
			BBKNN <- RunUMAP(BBKNN, nn.name = "bbknn", seed.use = 1)
			BBKNN <- FindNeighbors(BBKNN, reduction = 'umap', dims = 1:2, verbose = FALSE)
			BBKNN <- FindClusters(BBKNN, algorithm = 3, resolution = 0.01, verbose = FALSE)
			BBKNN$cluster <- Idents(BBKNN)
			
			#running Liger results
			params = setParams(method = "Liger", return = "Seurat")
			Liger = Merge(params, data)
			Liger <- RunUMAP(Liger, reduction = "iNMF", dims = 1:10, seed.use = 1)
			Liger <- FindNeighbors(Liger, reduction = 'umap', dims = 1:2, verbose = FALSE)
			Liger <- FindClusters(Liger, algorithm = 3, resolution = 0.01, verbose = FALSE)
			Liger$cluster <- Idents(Liger)
			
			#running Scanorama results
			params = setParams(method = "Scanorama", return = "Seurat")
			Scanorama = Merge(params, data)
			Scanorama <- RunUMAP(Scanorama, reduction = "scanorama", dims = 1:10, seed.use = 1)
			Scanorama <- FindNeighbors(Scanorama, reduction = 'umap', dims = 1:2, verbose = FALSE)
			Scanorama <- FindClusters(Scanorama, algorithm = 3, resolution = 0.01, verbose = FALSE)
			Scanorama$cluster <- Idents(Scanorama)
			
			#convert to SingleCellExpirement for metrics
			uncorrected = as.SingleCellExperiment(uncorrected)
			Seurat = as.SingleCellExperiment(Seurat)
			Harmony = as.SingleCellExperiment(Harmony)
			BBKNN = as.SingleCellExperiment(BBKNN)
		}
		

		#running ensemblemerge results
		test = EnsembleMerge(data, methods = c("Seurat", "Harmony", "Liger", "BBKNN", "Scanorama"), return = "Seurat", file = sprintf("Dataset_%i.csv", i))
		test = RunUMAP(test, graph = "EnsembleMerge")
		test <- FindNeighbors(test, reduction = 'umap', dims = 1:2, verbose = FALSE)
		test <- FindClusters(test, algorithm = 3, resolution = 0.01, verbose = FALSE)
		
		#convert to SingleCellExpirement for metrics
		test_graph = test[["RNA_nn"]]
		test = as.SingleCellExperiment(test)
		metadata(test)[["RNA_nn"]] = test_graph
		
		#SCmetrics testing pipeline
		config <- list(
		  label_field = "CellType",
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

		names = if (individual_methods) c("Uncorrected", "EnsembleMerge", "Seurat", "Harmony", "BBKNN", "Liger", "Scanorama") else c("EnsembleMerge")
		i = 1

		nmi_list = list()
		ari_cell_list = list()
		ari_batch_list = list()
		asw_cell_list = list()
		asw_batch_list = list()
		graph_con_list = list()
		
		methods = if (individual_methods) list(uncorrected, test, Seurat, Harmony, BBKNN, Liger, Scanorama) else list(test)

		for(tests in methods){
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
		  graph_con <- metrics(tests, new('MetricsGraphConnectivity'), label = ARGS[['batch_field']], graph = "RNA_nn")
		  graph_con_list = append(graph_con_list, graph_con)
		  sprintf('%s Score:', names[i]) %>% message()
		  sprintf('NMI: %.3f', nmi) %>% message()
		  sprintf('ARI(cell): %.3f', ari_cell) %>% message()
		  sprintf('ASW(cell): %.3f', asw_cell) %>% message()
		  sprintf('ARI(batch): %.3f', ari_batch) %>% message()
		  sprintf('ASW(batch): %.3f', asw_batch) %>% message()
		  sprintf('GC: %.3f', graph_con) %>% message()
		  i = i + 1
		}

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
		scores$seed = h
		scores_list = append(scores_list, scores)
	}
}

scores = Reduce(rbind, scores_list)
write.csv(scores, file = "fig9.csv")