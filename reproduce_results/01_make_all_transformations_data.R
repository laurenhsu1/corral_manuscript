source('../R/load_pkgs.R')
source('../R/utils.R')
source('../R/glmpca_feat_selection.R')

library(scRNAseq)

mthkeep <- c('std', 'varstab', 'anscombe_std', 'ft_x_std', 'ft_p_std',
             'ftalt_p', 'adj_pearson', 'deflate', 'deflate_95', 'deflate_98',
             'sqrt_ftres')

checkdir('../data')

save_dir <- '../data/0_make_all_transformations_data_50pcs/'

checkdir(save_dir)

ncomps <- 50

mcc <- 4

# setting up scelist

tscelist <- list(zhengmix4eq = sce_full_Zhengmix4eq(),
                 zhengmix4uneq = sce_full_Zhengmix4uneq(),
                 zhengmix8eq = sce_full_Zhengmix8eq())

tscelist[['zhengmix8eq']] <- tscelist[['zhengmix8eq']][which(rowSums(counts(tscelist[['zhengmix8eq']])) > 0),]
tscelist <- lapply(tscelist, logNormCounts)

new_sces_full <- list(baron_panc = BaronPancreasData('human', location = FALSE),
                      lawlor_panc = LawlorPancreasData(),
                      muraro_panc = MuraroPancreasData(ensembl = FALSE, location = FALSE),
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
  
  # # to just run with glmpca once
  # res$glmpca <- glmpca_routine(sce, npc = ncomps, seed = 987)
  
  res$CA <- CA_routine(sce, npc = ncomps)
  
  # adding in the deflate param comparisons
  res$deflate_95 <- compsvd_mod(test_corral_preproc(counts(sce), rtype = 'deflate', deflate_alpha = .95), ncomp = ncomps)
  res$deflate_98 <- compsvd_mod(test_corral_preproc(counts(sce), rtype = 'deflate', deflate_alpha = .98), ncomp = ncomps)
  
  # adding in pca on logcounts
  res$pca <- irlba::prcomp_irlba(x = t(logcounts(sce)), n = ncomps, center = TRUE, scale. = TRUE)
  
  # adding in sqrt + FT-residuals
  res$sqrt_ftres <- compsvd_mod(test_corral_preproc(counts(sce), rtype = 'sqrt_ftres'), ncomp = ncomps)
  
  res
}

decomp_list_CA <- mclapply(scelist, all_decomps, mc.cores = mcc)
saveRDS(decomp_list_CA, file = paste0(save_dir, 'decomp_list_CA.rds'))

# add glmpca

decomp_list <- list()

for(dat in names(scelist)){
  print(dat)
  sce <- scelist[[dat]]
  seeds <- list(1234, 12345, 987,9876,98765, 2020, 2021, 123, 31415, 4321)
  gorig_res <- mclapply(seeds, FUN = function(s) try(glmpca_routine(sce, npc = ncomps, seed = s)), mc.cores = mcc)
  names(gorig_res) <- paste0('glmpca_seed',seeds)
  decomp_list[[dat]] <- append(decomp_list_CA[[dat]], gorig_res)
}

saveRDS(decomp_list, file = paste0(save_dir, 'decomp_list.rds'))

# embedding list

decomps2embeds <- function(decomp_list, multi_glmpca = FALSE){
  res <- lapply(decomp_list[mthkeep], FUN = '[[', 'v')
  if(multi_glmpca){
    for(gnam in grep('glmpca',names(decomp_list), value = TRUE)){
      res[[gnam]] <- decomp_list[[gnam]][['factors']]
    }
  }
  else{
    res$glmpca <- decomp_list$glmpca$factors 
  }
  res$pca <- decomp_list$pca$x
  res
}

emb_list <- mclapply(decomp_list, decomps2embeds, multi_glmpca = TRUE, mc.cores = mcc)

saveRDS(emb_list, file = paste0(save_dir, 'emb_list.rds'))
