source('../R/load_pkgs.R')
source('../R/utils.R')
source('../R/glmpca_feat_selection.R')

library(scRNAseq)

mthkeep <- c('std', 'varstab', 'anscombe_std', 'ft_x_std', 'ft_p_std',
             'ftalt_p', 'adj_pearson', 'deflate', 'deflate_95', 'deflate_98')

save_dir <- '../data/0_make_all_transformations_data_7.2.21topfeats/'

checkdir(save_dir)

ncomps <- 30

mcc <- 4

# setting up scelist

try(load('../data/1_scelist_extended.rda')) # code to make this file in reproduce_figs/1_3-figures.Rmd

tscelist <- scelist_ext[c("zhengmix4eq","zhengmix4uneq", "zhengmix8eq",
                          'pbmcsca_10x','pbmcsca_dropseq','pbmcsca_indrops')]

new_sces_full <- list(baron_panc = BaronPancreasData('human', location = FALSE),
                      lawlor_panc = LawlorPancreasData(),
                      muraro_panc = MuraroPancreasData(ensembl = FALSE, location = FALSE),
                      campbell_brain = CampbellBrainData(location = FALSE),
                      chen_brain = ChenBrainData(location = FALSE),
                      darmanis_brain = DarmanisBrainData(location = FALSE),
                      xenopus = AztekinTailData()
)

new_sces_full <- lapply(new_sces_full, filter_low_sf)

new_sces_full <- lapply(new_sces_full, logNormCounts)

new_sces <- lapply(new_sces_full, scran_feat_selection)

scelist <- append(tscelist, new_sces)

saveRDS(scelist, file = paste0(save_dir, 'scelist.rds'))

# Function to get all decompositions, including time comparisons
# (time comparisons in glmpca_routine and CA_routine output)

all_decomps <- function(sce){
  res <- build_corral_list_only('sce', list('sce' = sce), ncores = 1, npc = ncomps)
  
  # seeds <- c(1234, 12345, 987,9876,98765, 2020, 2021, 123, 31415)
  # for(seed in seeds){
  #   res[[paste0('glmpca_seed',seed)]]<- glmpca_routine(sce, npc = ncomps, seed = seed)
  #   res[[paste0('glmpca_mem_seed',seed)]]<- glmpca_routine(sce, npc = ncomps, seed = seed, mb = 'memoized')
  #   res[[paste0('glmpca_stoch_seed',seed)]]<- glmpca_routine(sce, npc = ncomps, seed = seed, mb = 'stochastic')
  # }
  
  # to just run with glmpca once
  res$glmpca <- glmpca_routine(sce, npc = ncomps, seed = 987)
  
  res$CA <- CA_routine(sce, npc = ncomps)
  
  # adding in the deflate param comparisons
  res$deflate_95 <- compsvd_mod(test_corral_preproc(counts(sce), rtype = 'deflate', deflate_alpha = .95), ncomp = ncomps)
  res$deflate_98 <- compsvd_mod(test_corral_preproc(counts(sce), rtype = 'deflate', deflate_alpha = .98), ncomp = ncomps)
  
  # adding in pca on logcounts
  res$pca <- irlba::prcomp_irlba(x = t(logcounts(sce)), n = ncomps, center = TRUE, scale. = TRUE)
  
  res
}

decomp_list <- mclapply(scelist, all_decomps, mc.cores = mcc)

saveRDS(decomp_list, file = paste0(save_dir, 'decomp_list.rds'))

# embedding list

decomps2embeds <- function(decomp_list, multi_glmpca = FALSE){
  res <- lapply(decomp_list[mthkeep], FUN = '[[', 'v')
  if(multi_glmpca){
    for(gnam in grep){
      res[[gnam]] <- decomp_list[[gnam]][['factors']]
    }
  }
  else{
    res$glmpca <- decomp_list$glmpca$factors 
  }
  res$pca <- decomp_list$pca$x
  res
}

emb_list <- mclapply(decomp_list, decomps2embeds, mc.cores = mcc)

saveRDS(emb_list, file = paste0(save_dir, 'emb_list.rds'))

