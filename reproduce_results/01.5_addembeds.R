save_dir <- '../data/0_make_all_transformations_data_7.4.21topfeats/'

source('../R/load_pkgs.R')
source('../R/utils.R')
source('../R/glmpca_feat_selection.R')

library(scRNAseq)

mthkeep <- c('std', 'varstab', 'anscombe_std', 'ft_x_std', 'ft_p_std',
             'ftalt_p', 'adj_pearson', 'deflate', 'deflate_95', 'deflate_98')


datkeep <- c('zhengmix4eq', "zhengmix4uneq","zhengmix8eq", "baron_panc", 
             "muraro_panc", "lawlor_panc","chen_brain","darmanis_brain", 'xenopus')


checkdir(save_dir)

ncomps <- 50

mcc <- 4

scelist <- readRDS(file = paste0(save_dir, 'scelist.rds'))
old_decomp_list <- readRDS(file = paste0(save_dir, 'decomp_list.rds'))

seeds <- list(1234, 12345, 987, 9876, 98765, 2020, 2021, 123, 31415)
old_mth_keep <- c(paste0('glmpca_seed',seeds), 'CA', 'pca', 'deflate_95')

decomp_list <- lapply(old_decomp_list, FUN = function(x) x[old_mth_keep])

for(dat in names(scelist)){
  sce <- scelist[[dat]]
  print(dat)
  res <- try(build_corral_list_only(sce_name = 'sce', 
                                    scelist = list(sce = sce), 
                                    ncores = mcc, npc = ncomps))

  
    
  res <- append(res, decomp_list[[dat]])
  decomp_list[[dat]] <- res
  saveRDS(decomp_list, file = paste0(save_dir, 'decomp_list.rds'))
}


# embedding list

decomps2embeds <- function(decomp_list){
  res <- lapply(decomp_list[mthkeep], FUN = '[[', 'v')
  for(gnam in grep('glmpca', names(decomp_list), value = T)){
    res[[gnam]] <- decomp_list[[gnam]][['factors']]
  }
#  res$glmpca <- decomp_list$glmpca$factors
  res$pca <- decomp_list$pca$x
  res
}

emb_list <- lapply(decomp_list, decomps2embeds)

saveRDS(emb_list, file = paste0(save_dir, 'emb_list.rds'))

