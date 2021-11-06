sil_stats <- function(emb_list, labels, all_sils = FALSE, npc = ncol(emb_list[[1]])){
  labfactor <- as.factor(labels)
  n <- min(npc, min(unlist(lapply(emb_list[[nam]], ncol))))
  emb_list <- lapply(emb_list, FUN = function(x) x[,1:n])
  distlist <- lapply(emb_list, dist)
  cell_embs_sils <- lapply(distlist, 
                           cluster::silhouette, 
                           x = as.numeric(labfactor))
  
  if(all_sils == TRUE){
    return(cell_embs_sils)
  }
  
  cell_sil_summary <- lapply(cell_embs_sils, summary)
  
  cell_sil_clustasw <- lapply(cell_sil_summary, function(x) x$clus.avg.widths)
  clustasw_df <- Reduce(cbind, cell_sil_clustasw)
  colnames(clustasw_df) <- names(cell_sil_clustasw)
  rownames(clustasw_df) <- levels(labfactor)
  
  cell_sil_summary[['clustasw_df']] <- clustasw_df
  return(cell_sil_summary)
}
  
  
plot_sil_stats <- function(clustasw_df, returngg = TRUE, dataset = NULL){
  if(is.null(dataset)){
    plottitle = 'ASW by method, celltype'
  }
  else{
    plottitle = paste0('ASW by method, celltype: ', dataset)
  }
  melted <- melt(clustasw_df)
  colnames(melted) <- c('celltype','method','value')
  ggobj <- ggplot(melted, aes(x = celltype, fill = method, 
                                y = value)) +
    geom_bar(stat='identity', position = 'dodge') +
    labs(title = plottitle) +
    coord_flip()
  
  show(ggobj)
  
  if(returngg){
    return(ggobj)
  }
}

plot_asw_df <- function(asw_df, title = "ASW by dataset, preprocessing, and dimreduction method", returngg = TRUE){
  ggobj <- ggplot(asw_df, aes(x = inp_mat, y = ASW, fill = method)) + 
    geom_bar(stat='identity', position = 'dodge') + 
    facet_wrap(~ dataset) + ggtitle(title)
  if(returngg){
    return(ggobj)
  }
}

# For use with silstats(..., all_sils = FALSE)
flatten_sils <- function(sil_list, datname){
  sil_list <- sil_list[which(names(sil_list) != 'clustasw_df')]
  asws <- matrix(unlist(lapply(sil_list, '[','avg.width')), ncol = 1)
  rownames(asws) <- paste0(datname,'_', names(sil_list))
  dat <- cbind(asws, 
               data.frame(Reduce(rbind,strsplit(names(sil_list), split = '_'))),
               rep(datname, length(sil_list)))
  colnames(dat) <- c('ASW','inp_mat', 'method', 'dataset')
  return(dat)
}

# For use with silstats(..., all_sils = TRUE)
flatten_full_sils <- function(sil_list, datname){
  splitnames <- strsplit(names(sil_list), split = '_')
  names(splitnames) <- names(sil_list)
  for(sil in names(sil_list)){
    sil_list[[sil]] <- data.frame(sil_list[[sil]][,])
    sil_list[[sil]]$inp_mat <- splitnames[[sil]][1]
    sil_list[[sil]]$method <- splitnames[[sil]][2]
  }
  full_sils <- Reduce(rbind, sil_list)
  full_sils$dataset <- datname
  return(full_sils)
}


plot_sw_distribution <- function(sw_df, title = 'Distribution of individual sil. widths', returngg = TRUE){
  ggobj <- ggplot(sw_df, aes(x = inp_mat, y = sil_width, fill = method)) + 
    geom_boxplot() + 
    facet_wrap(~ dataset) + ggtitle(title)  
  if(returngg){
    return(ggobj)
  }
}

supscale <- function(mat, row = FALSE){
  if(row){
    rw <- rowSums(mat)
    supscaled <- sweep(mat, MARGIN = 1, STATS = rw, FUN = '/')
  }
  else{
  cs <- colSums(mat)
  supscaled <- sweep(mat, MARGIN = 2, STATS = cs, FUN = '/')
  }
  return(supscaled)
}

projemb <- function(gene_embeds, d, newmat){
  feat_embeds <- gene_embeds %*% diag(d)
  procmat <- t(supscale(newmat))
  return(procmat %*% feat_embeds)
}

projcorral <- function(corral_obj, newmat){
  projemb(gene_embeds = corral_obj$u, 
          d = corral_obj$d, 
          newmat = newmat)
}



moonplot_corral <- function(corral_obj, color_vec, text_vec, nfeat = 20, xpc = 1, ellip = TRUE, pts = FALSE, plot_title = 'Moonplot', coords = c('svd','PC','SC','PCSC', 'svdSC'), brewpalette = 'Paired', ptsize = .5){
  n <- nfeat
  
  coords <- match.arg(coords, c('svd','PC','SC', 'PCSC', 'svdSC'))
  
  umat <- 'u'
  vmat <- 'v'
  
  if(coords == 'PCSC'){
    # this mode doesn't work right now
    umat <- paste0('SC', umat)
    vmat <- paste0('PC', vmat)
  }
  else if(coords == 'svdSC'){
    vmat <- paste0('PC', vmat)
  }
  else if(coords != 'svd'){
    umat <- paste0(coords, umat)
    vmat <- paste0(coords, vmat)
  }
  
  rad <- sort(abs(corral_obj[[vmat]][,xpc:(xpc + 1)]), decreasing = T)[10]
  border = rad*.15
  
  moonplot <- ggplot(data = data.frame(corral_obj[[vmat]][,xpc:(xpc + 1)], color_vec = color_vec), 
                     aes(x = X1, y = X2, color = color_vec)) + 
    geom_hline(yintercept = 0, color = 'gray') + 
    geom_vline(xintercept = 0, color = 'gray')
  
  if(ellip){
    moonplot <- moonplot + stat_ellipse()
  }
  if(pts){
    moonplot <- moonplot + geom_point(size = ptsize)
  }
  
  gene_dists <- sqrt(corral_obj[[umat]][,xpc]^2 + corral_obj[[umat]][,xpc + 1]^2)
  gene_dists_ordinds <- order(gene_dists, decreasing = TRUE)
  
  inflgenes <- corral_obj[[umat]][gene_dists_ordinds[1:n],]
  rownames(inflgenes) <- text_vec[gene_dists_ordinds][1:n]
  
  gene_df <- data.frame(rad*inflgenes[,xpc:(xpc + 1)]/gene_dists[gene_dists_ordinds[1:n]])
  gene_df$Name <- rownames(gene_df)
  gene_df$dists <- gene_dists[gene_dists_ordinds[1:n]]
  gene_df$sizes <- gene_df$dists * (1 / (max(gene_df$dists) - min(gene_df$dists)))
  gene_df$sizes <- gene_df$sizes - min(gene_df$sizes) + 1
  gene_df$drot <- 90 - atan(gene_df$X1/gene_df$X2) * (180/pi)
  gene_df$drot[which(gene_df$X2 < 0)] <- 180 + gene_df$drot[which(gene_df$X2 < 0)]
  
  moonplot <- moonplot + coord_fixed(ratio = 1) + 
    xlim(-1*(rad + border), rad + border) + ylim(-1*(rad + border), rad + border) +
    geom_function(fun = function(x) sqrt(rad^2 - x^2), inherit.aes = FALSE, xlim = c(-rad, rad), color = 'gray') +
    geom_function(fun = function(x) -1*sqrt(rad^2 - x^2), inherit.aes = FALSE, xlim = c(-rad, rad), color = 'gray') +
    #geom_circle(aes(x0 = 0, y0 = 0, r = rad), inherit.aes = FALSE) + 
    theme_void() +
    guides(size = FALSE) + labs(color = 'cell type') + ggtitle(plot_title) + 
    scale_color_brewer(palette=brewpalette) + scale_size(range = c(2,4))
  
  if(nrow(gene_df[which(gene_df$X1 > 0),] > 0)){
    moonplot <- moonplot + geom_text(data = gene_df[which(gene_df$X1 > 0),], 
                                     aes(x = X1, y = X2, label=as.character(Name), size = sizes, angle = drot),
                                     hjust = 0, vjust = 0, inherit.aes = FALSE,
                                     position=position_jitter(width=.01,height=.01) )
  }
  
  if(nrow(gene_df[which(gene_df$X1 < 0),] > 0)){
    moonplot <- moonplot + geom_text(data = gene_df[which(gene_df$X1 < 0),], 
                                     aes(x = X1, y = X2, label=as.character(Name), size = sizes, angle = I(drot + 180)),
                                     hjust = 1, vjust = 0, inherit.aes = FALSE,
                                     position=position_jitter(width=.01,height=.01))
  }
  
  return(moonplot)
}


biplot_corral <- function(corral_obj, color_vec, text_vec, feat_name = '(genes)', nfeat = 20, xpc = 1, plot_title = 'Biplot', text_size = 2, xjitter = .005, yjitter = .005, coords = c('svd','PC','SC','PCSC', 'svdSC'), scale_weights = FALSE,...){
  n <- nfeat
  coords <- match.arg(coords, c('svd','PC','SC', 'PCSC', 'svdSC'))
  
  umat <- 'u'
  vmat <- 'v'
  
  if(coords == 'PCSC'){
    # this mode doesn't work right now
    umat <- paste0('SC', umat)
    vmat <- paste0('PC', vmat)
  }
  else if(coords == 'svdSC'){
    vmat <- paste0('SC', vmat)
  }
  else if(coords != 'svd'){
    umat <- paste0(coords, umat)
    vmat <- paste0(coords, vmat)
  }
  
  gene_dists <- sqrt(corral_obj[[umat]][,xpc]^2 + corral_obj[[umat]][,xpc + 1]^2)
  gene_dists_ordinds <- order(gene_dists, decreasing = TRUE)
  
  scaling_factor <- 1
  if(scale_weights){
    cell_dists <- sqrt(corral_obj[[vmat]][,xpc]^2 + corral_obj[[vmat]][,xpc + 1]^2)
    scaling_factor <- quantile(as.numeric(cell_dists), .9) / gene_dists[gene_dists_ordinds[n]]
  }
  
  inflgenes <- corral_obj[[umat]][gene_dists_ordinds[1:n],] * scaling_factor
  rownames(inflgenes) <- text_vec[gene_dists_ordinds][1:n]
  biplot_labs_filt <- c(color_vec, rep(feat_name,n))
  biplot_dat_filt <- rbind(corral_obj[[vmat]], inflgenes)
  bipfilt_gg <- plot_embedding(biplot_dat_filt, 
                               color_vec = biplot_labs_filt, returngg = T,
                               dimname = 'corral', xpc = xpc, showplot = F,
                               plot_title = plot_title, ...)
  
  bipfilt_gg$data$Name <- rownames(bipfilt_gg$data)
  
  bipfilt_gg <- bipfilt_gg + 
    geom_text(aes(label=ifelse(color_vec == feat_name,as.character(Name),'')),
              hjust=0,vjust=0, 
              size = text_size, 
              position=position_jitter(width=xjitter,height=yjitter))
  
  return(bipfilt_gg)
}

ridge_plot <- function(mat, gene, color_vec, legend = TRUE, plot_title = gene){
  ridge_df <- data.frame(ct = color_vec,
                         log10_exp = log10(1 + mat[gene,]))
  ggobj <- ggplot(ridge_df, aes(x = log10_exp, y = ct, fill = ct)) +
    geom_density_ridges() +
    theme_ridges() + ggtitle(plot_title)
  
  if(!legend){
    ggobj <- ggobj + theme(legend.position = "none")
  }
  
  return(ggobj)
}


add_HVGsub <- function(scelist, modelfunction = modelGeneVar, fdr.thresh = .05){
  # adding in HVG subsetted SCEs to the sce list
  for(i in 1:length(scelist)){
    dat <- scelist[[i]]
    fns <- names(scelist)
    genevarmod <- modelfunction(dat)
    hvg_inds <- getTopHVGs(genevarmod, fdr.threshold = fdr.thresh)
    scelist[[paste0(fns[i],'_HVGsub')]] <- dat[hvg_inds,]
  } 
  return(scelist)
}

test_corral_preproc <- function(inp, rtype = c('ft_p','ft_x','varstab','ftalt_p','neyman','kullback','deflate','anscombe', 'std', 'anscombe_scale','ft_p_scale','ft_x_scale','anscombe_std', 'ft_x_std', 'adj_pearson','ft_p_std','varstab_noCA','varstab_scale'), deflate_alpha = .9){
  rtype <- match.arg(rtype, c('ft_p','ft_x','varstab','ftalt_p','neyman','kullback','deflate','anscombe', 'std', 'anscombe_scale','ft_p_scale','ft_x_scale','anscombe_std', 'ft_x_std', 'adj_pearson','ft_p_std','varstab_noCA','varstab_scale'))
  if(!is(inp, "dgCMatrix")){
    sp_mat <- Matrix::Matrix(inp, sparse = TRUE)
  } else {sp_mat <- inp}
  N <- sum(sp_mat)
  xmat <- sp_mat
  pmat <- sp_mat / N
  ws <- get_weights(sp_mat)
  row.w <- ws$row.w
  col.w <- ws$col.w
  row.sum <- rowSums(xmat)
  col.sum <- colSums(xmat)
  expectedp <- row.w %*% t(col.w)
  expectedx <- row.sum %*% t(col.sum)
  if(rtype == 'std'){
    return(corral_preproc(inp = sp_mat))
  }
  if(rtype == 'ft_p'){
    return((pmat^.5 - expectedp^.5)^2)
  }
  if(rtype == 'ft_x'){
    return(xmat^.5 + (xmat + 1)^.5)
  }
  if(rtype == 'varstab'){
    return(corral_preproc(xmat^.5))
  }
  if(rtype == 'varstab_noCA'){
    return((xmat^.5))
  }
  if(rtype == 'varstab_scale'){
    return(t(scale(t(xmat^.5))))
  }
  if(rtype == 'ftalt_p'){
    return(pmat^.5 + (pmat + 1/N)^.5 - (4*expectedp + 1/N)^.5)
  }
  if(rtype == 'neyman'){
    res <- pmat - expectedp
    res <- res / pmat^.5
    return(res)
  }
  if(rtype == 'kullback'){
    return(expectedp * log(expectedp / pmat))
  }
  if(rtype == 'deflate'){
    corral_pp <- corral_preproc(inp)
    return(sign(corral_pp) * (abs(corral_pp)^deflate_alpha))
  }
  if(rtype == 'anscombe'){
    return(2 * (xmat + 3/8)^.5)
  }
  if(rtype == 'anscombe_scale'){
    return(t(scale(2 * (t(xmat) + 3/8)^.5)))
  }
  if(rtype == 'anscombe_std'){
    return(corral_preproc(2 * (xmat + 3/8)^.5))
  }
  if(rtype == 'ft_p_scale'){
    return(t(scale(t(pmat^.5 - expectedp^.5)^2)))
  }
  if(rtype == 'ft_p_std'){
    return(corral_preproc((pmat^.5 - expectedp^.5)^2))
  }
  if(rtype == 'ft_x_scale'){
    return(t(scale(t((xmat^.5 + (xmat + 1)^.5)))))
  }
  if(rtype == 'ft_x_std'){
    return(corral_preproc(xmat^.5 + (xmat + 1)^.5))
  }
  if(rtype == 'adj_pearson'){
    ppmat <- corral_preproc(inp = xmat)
    ppmat <- ppmat * sqrt(1 - row.w)
    ppmat <- sweep(ppmat, 2, sqrt(1 - col.w), "*")
    return(ppmat)
  }
}


plot_embedding_continuous <- function(emb, feat_vec, plot_title, xpc = 1, ypc = xpc + 1, ptsize = .5, xlab = paste0('dim',xpc), ylab = paste0('dim',ypc), colhi = 'red', collo = 'blue'){
  grad_df <- data.frame(cbind(emb[,xpc], emb[,ypc], t(feat_vec)))
  names(grad_df) <- c('xpc','ypc','feat_level')
  ggplot(data = grad_df, aes(x = xpc, y = ypc, color = feat_level)) +
    theme_classic() + geom_point(size = ptsize) +
    geom_hline(yintercept=0, color = 'gray') + geom_vline(xintercept=0, color = 'gray') + 
    scale_colour_gradient(low=collo, high=colhi) + 
    labs(x = xlab, 
         y = ylab, 
         title = plot_title) + 
    theme(axis.text.x = element_text(size = rel(1.4), colour = 'black'), axis.text.y = element_text(size = rel(1.4), colour = 'black'), axis.title = element_text(size = rel(1.4), colour = 'black'))
}

compsvd_mod <- function(mat, method = c('irl','svd'), ncomp = 30, ...){
  method <- match.arg(method, c('irl','svd'))
  if(method == 'irl'){
    result <- irlba::irlba(mat, nv = ncomp, ...)
  }
  else if(method == 'svd'){
    result <- svd(mat, nv = ncomp, nu = ncomp, ...)
  }
  else {
    print('Provided method was not understood; used irlba.')
    result <- irlba::irlba(mat, nv = ncomp, ...)
  }
  result[['eigsum']] <- sum(mat^2)
  if(!is.null(rownames(mat))){rownames(result$u) <- rownames(mat)}
  if(!is.null(colnames(mat))){rownames(result$v) <- colnames(mat)}
  
  # add percent variance explained
  pct_var_exp <- t(data.frame('percent.Var.explained' = result$d^2 / result$eigsum))
  colnames(pct_var_exp) <- paste0(rep('PC',ncomp),seq(1,ncomp,1))
  result$pct_var_exp <- rbind(pct_var_exp,t(data.frame('cumulative.Var.explained' = cumsum(pct_var_exp[1,]))))
  
  return(result)
}


seurat2sce <- function(seurat_obj, slotname = 'counts'){
  count_mtx <- GetAssayData(object = seurat_obj, slot = slotname)
  coldat <- seurat_obj@meta.data
  sce <- SingleCellExperiment(counts = count_mtx, colData = DataFrame(coldat))
  rowData(sce)$symbol <- rownames(sce)
  return(sce)
}


plot_var_exp <- function(corral_list, whichvar = c('cumulative.Var.explained','percent.Var.explained')){
  npc <- length(corral_list[[1]]$d)
  vlist <- lapply(corral_list, FUN = '[[', 'pct_var_exp')
  vlist <- lapply(names(corral_list), FUN = function(x) data.frame(cbind(vlist[[x]][whichvar,], 1:npc ,rep(x, npc))))
  df <- Reduce(rbind, vlist)
  names(df) <- c('cumvar','ind', 'mth')
  df$cumvar <- as.numeric(df$cumvar)
  df$ind <- as.numeric(df$ind)
  ggplot(df, aes(x = ind, y = cumvar, color = mth)) + geom_line() + theme_classic()
}

build_corral_list_only <- function(sce_name, scelist = tscelist, mthkeep = c('std', 'varstab', 'anscombe_std', 'ft_x_std', 'ft_p_std','ftalt_p', 'adj_pearson', 'deflate'), ncores = 6, npc = 50){
  print(sce_name)
  sce <- scelist[[sce_name]]
  inpmat <- counts(sce)
  
  corral_ppcompmats <- mclapply(as.list(mthkeep), FUN = function(x) test_corral_preproc(inpmat, x), mc.cores = ncores)
  names(corral_ppcompmats) <- mthkeep
  
  corral_ppcomp <- list()
  
  corral_ppcomp <- mclapply(corral_ppcompmats, compsvd_mod, ncomp = npc, mc.cores = ncores)
  
  return(corral_ppcomp)
}

glmpca_routine <- function(sce, nfeat = 1500, npc = 50, devmod = 'poisson', seed = 1234, mb = 'none'){
  set.seed(seed)
  runtime <- system.time({
    sce_filt <- filterDev(sce, nkeep = min(nrow(sce),nfeat), dev = devmod)
    res <- glmpca::glmpca(counts(sce_filt), L = npc, fam = 'poi', minibatch = mb)
  })
  res$runtime <- runtime
  res
}

CA_routine <- function(sce, nfeat = 1500, npc = 50, devmod = 'poisson'){
  runtime <- system.time({
    sce_filt <- filterDev(sce, nkeep = min(nrow(sce),nfeat), dev = devmod)
    res <- corral(counts(sce_filt), ncomp = npc)
  })
  res$runtime <- runtime
  res
}


filter_low_sf <- function(sce, minreads = 100){
  sce[,which(colSums(counts(sce)) > minreads)]
}


clusters_acrossPCs <- function(emb_list, ct_vec, pcs = c(2, seq(5, ncol(emb_list[[1]]), 5)), dataset_name, blusparam = NNGraphParam()){
  final_list <- list()
  ARI_df <- data.frame()
  try({
    for(pcnum in pcs){
      emb_list_sub <- lapply(emb_list, FUN = function(X, pc_ind = pcnum) X[,1:pc_ind])
      clust_list <- lapply(emb_list_sub, clusterRows, BLUSPARAM = blusparam)
      final_list[[paste0(pcnum, 'PC')]] <- clust_list
      ARI_list <- lapply(clust_list, FUN = pairwiseRand, ref = as.factor(ct_vec), mode = 'index')
      pc_df <- cbind(data.frame(ARI = unlist(ARI_list), 
                                method = names(ARI_list),
                                npc = rep(pcnum, length(ARI_list))))
      ARI_df <- rbind(ARI_df, pc_df)
    }
  })
  #names(ARI_df_corral) <- c('ARI','method')
  ARI_df$dataset <- rep(dataset_name, nrow(ARI_df))
  
  final_list$ARI_df <- ARI_df
  
  return(final_list)
}


checkdir <- function(dirname) ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)


scran_feat_selection <- function(sce, fdr_thresh = .1, nfeat = 3000, mth = 'fdr'){
  genevarmod <- modelGeneVar(sce)
  
  if(mth == 'fdr'){
    hvg_inds <- getTopHVGs(genevarmod, fdr.threshold = fdr_thresh)
  }
  else{
    hvg_inds <- getTopHVGs(genevarmod, n = nfeat)
  }
  
  sce[hvg_inds,]
}


