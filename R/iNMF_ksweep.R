#' Run iNMF via liger
#' @export

#' @title Multiplies two numbers, adds a constant, and returns the result
#' @name iNMF_ksweep
#' @param seurat.object A seurat object containing raw RNA counts
#' @param assay_slot Which seurat assay_slot to perform NMF on, defaults to RNA counts. For advanced users, see documentation for discussion of performing NMF on SCTransformed count data.
#' @param Sample.thresh Minimum number of samples a given cell type must be present across for NMF
#' @param Type.thresh Minimum number of cells per cell type in n Sample.thresh samples to be included for NMF
#' @param Batch Sequencing batch
#' @param scale.factor Scale.factor for count normalization
#' @param log.norm Perform log normalization on counts
#' @param k.min,k.max Max and min k's for k.sweep
#' @param n.reps Number of reps for consensus nmf
#' @param nn.pt Resolution for nearest neighbors clustering for consensus nmf
#' @param max.cores N cores to use for parallelization.
#' @param output.reps Save individual NMF reps. Default FALSE.
#' @param return.scale.data Save scaled data. Defualt FALSE
#'
#' @return a list containing k sweep consensus NMF results
#'
#' @import Seurat
#' @import rliger
#' @import Matrix
#' @import matrixStats
#' @import parallel
#' @import parallelly
#' @import cluster
#' @import stats
#' @import rlist
#' @import liger
#' @import RANN
#' @export


iNMF_ksweep <- function(seurat.object, assay_slot = 'RNA', Type.thresh = 100, Sample.thresh = 10, Batch = TRUE,
                            scale.factor = 10000, log.norm = TRUE, k.min = 2, k.max = 40, n.reps = 20, nn.pt = 0.3,
                            max.cores = parallelly::availableCores(), output.reps = FALSE, return.scale.data = F){
  if (!"Sample"%in%colnames(seurat.object@meta.data) |
      !"CellType"%in%colnames(seurat.object@meta.data)){
    stop("metadata slot must contain colnames 'Sample' and 'CellType'")
  }
  if (Batch){
    if (!"Batch"%in%colnames(seurat.object@meta.data)){
      stop("metadata slot must contain colname 'Batch'")
    }
  }
  if (round(n.reps*nn.pt)<=1){
    stop("n.reps or nn.pt too low to find consensus results")
  }
  if (!assay_slot%in%seurat.object@assays){
    stop('Supplied assay slot not in sequencing object')
  DefaultAssay(seurat.object)=assay_slot
  }
  # Process Seurat object and determine which cell types to analyze
  # Only run analysis on cell types that represented across at least 10 samples with at least 100 cells per sample (or adjust inputs above)
  cell.types <- names(which(colSums(table(seurat.object@meta.data[,c("Sample","CellType")])>=Type.thresh)>=Sample.thresh))
  if (length(cell.types)==0){
    stop("Not enough cells to meet input Sample/Type thresholds")
  } else if (length(cell.types)==1){
    warning("Only one cell type met input Sample/Type thresholds")
  }
  print(paste0("Running consensus iNMF on the following ", length(cell.types), " cell types: ", paste(cell.types, collapse = ", ")))

  # SET UP LIGER OBJECTS AND RUN K SWEEP
  print(paste0('Creating Liger objects and running k sweep...'))
  names(cell.types) = cell.types
  res <- lapply(cell.types, function(x){
    print(paste0("Cell type ", x))
    seu_sub <- subset(seurat.object, CellType==x)
    if (Batch){
      if (length(table(seurat.object$Batch))>1){
        Batch.list <- SplitObject(seu_sub, split.by="Batch")
        Liger.setup=list()
        for (j in 1:length(Batch.list)){
          Liger.setup[[j]]=GetAssayData(Batch.list[[j]], slot = 'counts') #updating so single function interoperable with multiple data slots, using improved interaction functions for SeuratObject
        }
        names(Liger.setup)=names(Batch.list)
        Liger <- createLiger(Liger.setup)
      }
    } else {
      seu_sub <- AddMetaData(seu_sub, col.name = 'Batch', metadata = rep('Batch1', nrow(seu_sub@meta.data)))
      Liger <- createLiger(list(Batch1 = GetAssayData(seu_sub, slot = 'counts'))) #even for SCT normalization, residuals converted back to 'corrected' UMI counts, pre log-normalization per SeuratV5 Documentation
    }
    # normalize so all cells have same total counts
    Liger <- rliger::normalize(Liger) #still need to select variable genes for object regardless of normalization method
    Liger <- selectGenes(Liger)
    if (log.norm){
      # log normalize (this combined with the normalize step above is the same as LogNormalize in Seurat)
      for (k in 1:length(Liger@norm.data)){
        Liger@norm.data[[k]]=as.sparse(as.matrix(log1p(Liger@norm.data[[k]]*scale.factor)))}
    } else{
      for (k in 1:length(Liger@norm.data)){
        Liger@norm.data[[k]]=as.sparse(as.matrix(Liger@norm.data[[k]]))} #should always be counts slot regardless of normalization/pre-processing method
    }
    # scale without centering
    Liger <- rliger::scaleNotCenter(Liger) # does need to be re-scaled (if integrated, pull from integrated assay's data slot not integrated assay's 'scale.data' slot)

    # adjust online_iNMF parameters based on dataset size
    minibatch_size <- min(table(seu_sub@meta.data$Batch))
    if (minibatch_size > 5000) {
      minibatch_size = 5000} else if (minibatch_size > 1000) {
        minibatch_size = floor(minibatch_size/1000)*1000} else if (minibatch_size > 500) {
          minibatch_size = floor(minibatch_size/500)*500} else if (minibatch_size > 100) {
            minibatch_size = floor(minibatch_size/100)*100} else minibatch_size = minibatch_size
    h5_chunk_size <- min(1000, minibatch_size)

    Liger_list = list()
    for (k in c(k.min:k.max)){ #because for loop is outside forking step where individual reps are run,
      #each previous k in Liger_list (with reps.k) is stored and then exported to the forked workers!!!
      #exponentially increasing memory load. that's why she crashes!!!!!!
      print(paste0('Running K = ', k))
      reps.k = rep(k, n.reps)
      names(reps.k) = paste0("rep", 1:length(reps.k))
      cl <- makeCluster(max.cores)
      clusterExport(cl=cl, c('Liger', 'reps.k', 'minibatch_size', 'h5_chunk_size'),
                    unclass(lsf.str(envir = asNamespace("deciphR"), all = T)),
                    envir = as.environment(asNamespace("deciphR")))
      Liger_list[[paste0("R", k)]] = clusterMap(cl, LIGER.run, K = reps.k, Liger = c(Liger),
                                              minibatch_size = minibatch_size, h5_chunk_size = h5_chunk_size,
                                              SIMPLIFY = F)
      stopCluster(cl)
      for (i in 1:length(Liger_list[[paste0("R", k)]])){
        colnames(Liger_list[[paste0("R", k)]][[i]][["W"]]) =  paste0("Rep", i, "_", colnames(Liger_list[[paste0("R", k)]][[i]][["W"]]))
        colnames(Liger_list[[paste0("R", k)]][[i]][["H"]]) = paste0("Rep", i, "_", colnames(Liger_list[[paste0("R", k)]][[i]][["H"]]))
        Liger_list[[paste0("R", k)]][[i]][["V"]] = lapply(Liger_list[[paste0("R", k)]][[i]][["V"]], function(y) {
          colnames(y) = paste0("Rep", i, "_", colnames(y))
          return(y)})
      }
      print(paste0('Done with ', x, ' K = ', k))
    }

    # FIND CONSENSUS RESULTS
    print(paste0('Finding consensus results'))
    Consensus.results = mclapply(Liger_list, function(x){
      W_list = lapply(x, '[[', 'W')
      W_2norm = rapply(W_list, f = function(y) {apply(y, MARGIN = 2, FUN = function(y){norm(y, type ="2")})},
                       how = "replace")
      W_list = rapply(W_list, f = function(y) {t(apply(y, MARGIN = 2, FUN = function(y) {y/norm(y, type = "2")}))},
                      how = "replace")
      W.mat = list.rbind(W_list)
      kmeans_clusts = dim(W.mat)[1]/n.reps

      # find mean distance to nearest neighbors
      nn.dist = rowMeans(RANN::nn2(as.matrix(W.mat), k = n.reps)$nn.dists[,c(1:round(n.reps*nn.pt))+1])
      names(nn.dist) = rownames(W.mat)

      # filter outliers
      min = diff(quantile(nn.dist, probs = c(0.25, 0.75)))*0.5 + quantile(nn.dist, probs = c(0.75))
      if (sum(nn.dist>min)>0) {
        W.mat.filt = W.mat[-which(nn.dist>min),]
      } else {
        W.mat.filt = W.mat
      }

      # perform kmeans clustering
      km.res = kmeans(W.mat.filt, centers = kmeans_clusts, nstart = 100, iter.max = 100)

      # find consensus W (shared gene loadings)
      W_consensus = matrix(nrow = kmeans_clusts, ncol = ncol(W.mat.filt))
      for (i in seq(kmeans_clusts)){
        row.ind = which(km.res$cluster==i)
        if (length(row.ind) > 1){
          W_consensus[i,] = colMedians(W.mat.filt[row.ind,])
        } else W_consensus[i,] = W.mat.filt[row.ind,]
      }
      rownames(W_consensus) = paste0("R", kmeans_clusts, "_Program", seq(kmeans_clusts))
      colnames(W_consensus) = colnames(W.mat.filt)

      # find consensus V (batch-specific gene loadings)
      V_list = lapply(x, `[[`, "V")
      for (i in names(V_list)){
        V_list[[i]] = lapply(V_list[[i]], function(x){t(t(x)*(1/W_2norm[[i]]))})
      }
      V.mat = t(list.cbind(lapply(V_list, list.rbind)))
      if (sum(nn.dist>min)>0) {
        V.mat.filt = V.mat[-which(nn.dist>min),]
      } else {
        V.mat.filt = V.mat
      }
      V_consensus = matrix(nrow = kmeans_clusts, ncol = ncol(V.mat.filt))
      for (i in seq(kmeans_clusts)){
        row.ind = which(km.res$cluster==i)
        if (length(row.ind) > 1){
          V_consensus[i,] = colMedians(V.mat.filt[row.ind,])
        } else V_consensus[i,] = V.mat.filt[row.ind,]
      }
      rownames(V_consensus) = paste0("R", kmeans_clusts, "_Program", seq(kmeans_clusts))
      colnames(V_consensus) = colnames(V.mat.filt)
      Batch.info = 1:length(V_list$rep1)
      names(Batch.info) = sort(names(V_list$rep1))
      V_consensus = lapply(Batch.info, function(x){
        V = V_consensus[,((x-1)*length(colnames(W_consensus))+1):(x*length(colnames(W_consensus)))]
      })

      # Solve for H (activity program expression scores) using consensus W and V initializations
      H_consensus_list = lapply(Batch.info, function(x){
        H = solveNNLS(rbind(t(W_consensus) + t(V_consensus[[x]]),
                            sqrt(5) * t(V_consensus[[x]])), rbind(t(Liger@scale.data[[x]]), matrix(0,
                                                                                                   dim(W_consensus)[2], dim(Liger@raw.data[[x]])[2])))
      })
      H_consensus_list = lapply(H_consensus_list, t)
      H_consensus = list.rbind(H_consensus_list)
      rownames(H_consensus) = rownames((lapply(x, `[[`, "H"))[[1]])
      colnames(H_consensus) = paste0("R", kmeans_clusts, "_Program", seq(kmeans_clusts))
      consensus_res = list(H = H_consensus, W = W_consensus, V = V_consensus)
      return(consensus_res)
    })
    print(paste0('Done with ', x))

    res = list()
    if (return.scale.data){
      res[["scale.data"]] = Liger@scale.data
    }
    if (output.reps) {
      res[["Liger.list"]] = Liger_list
    }
    res[["consensus.results"]] = Consensus.results
    return(res)
  })
  return(res)
}

LIGER.run <- function(K, Liger, minibatch_size, h5_chunk_size, seed = NULL){
  # sample seeds and report for reproducibility
  if (is.null(seed)) {seed = sample(1:1000000, 1)}
  Liger <- online_iNMF(Liger, k = K, max.epochs=10, seed = seed, miniBatch_size = minibatch_size, h5_chunk_size = h5_chunk_size, verbose = F)
  cat(".")
  W_res = t(Liger@W)
  colnames(W_res) = paste0("R", K, "_Program", 1:K)
  H_res = list.rbind(Liger@H)
  colnames(H_res) = paste0("R", K, "_Program", 1:K)
  V_res = Liger@V
  V_res = lapply(V_res, function(y){
    rownames(y) = paste0("R", K, "_Program", 1:K)
    y = t(y)
    return(y)
  })
  params = c(minibatch_size = minibatch_size, h5_chunk_size = h5_chunk_size)
  return(list(W = W_res, H = H_res, V = V_res))
}
