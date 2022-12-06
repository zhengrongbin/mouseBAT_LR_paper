args=commandArgs(T)
library(Seurat)
library(CellChat)

#cond = args[1]


#RT.obj = readRDS(paste0(cond, '_seurat.rds'))
#RT.obj@active.assay = 'RNA'
#inputDataf = args[1]## expression matrix, csv
#metaFile = args[2] ## meta matrix, csv
#cond = args[3]

cond = 'allcell'

inputData = NULL
meta = NULL
for (tem in c('TN', 'RT', 'cold2', 'cold7')){
print(tem)
inputDataf = paste0('scanpy_exp/', tem, '_sc_exp_updated.csv.gz')
tmp = read.csv(gzfile(inputDataf), row.names = 1, sep = '\t')
metaFile =  paste0('scanpy_exp/', tem, '_sc_meta_updated.csv')
tmpm = read.csv(metaFile, row.names = 1, sep = '\t')

inputData = rbind(inputData, tmp)
meta = rbind(meta, tmpm)
}

rownames_adj = NULL
for (i in 1:nrow(inputData)){
	x = rownames(inputData)[i]
	x = strsplit(x, '\\-')[[1]][1]
	rownames_adj = c(rownames_adj, paste0(x, '-', i))
}

rownames(inputData) = rownames_adj
rownames(meta) = rownames_adj

# combine labels as cell_type~cond
meta$cell_type_new = apply(meta[,c('cell_type', 'cond')], 1, function(x) paste0(x, collapse = '~'))


## matrix
inputData = t(as.matrix(inputData))
meta = as.matrix(meta)


## create cellchat object
#cellchat <- createCellChat(object = RT.obj@assays$RNA@data, meta = RT.obj@meta.data, group.by = "cluster_id")
cellchat <- createCellChat(object = inputData, meta = meta, group.by = 'cell_type_new')

## cellchat
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- CellChatDB#subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multicore", workers = 8) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

##Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

## Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
## Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, paste0('cellchat_res/', cond, '_updated_celltype_cond_updated_cellchat.newrun.rds'))
