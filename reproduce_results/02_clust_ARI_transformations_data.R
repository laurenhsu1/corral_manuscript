save_dir <- '../data/0_make_all_transformations_data_50pcs/'
source('../R/load_pkgs.R')
source('../R/utils.R')
source('../R/glmpca_feat_selection.R')

# number of cores to use
mcc <- 5

# reading in data from 0_make_all_transformations_data

scelist <- readRDS(file = paste0(save_dir, 'scelist.rds'))
decomp_list <- readRDS(file = paste0(save_dir, 'decomp_list.rds'))
emb_list <- readRDS(file = paste0(save_dir, 'emb_list.rds'))


# making ctlist

old_ctlist <- lapply(scelist[1:3], FUN = function(sce) sce$phenoid)

xeno.sce <- scelist$xenopus
xeno.sce$ct <- xeno.sce$cluster
xeno.sce$ct[grep('Erythrocyte', xeno.sce$ct)] <- 'Erythrocyte'
xeno.sce$ct[grep('Interneuron', xeno.sce$ct)] <- 'Interneuron'
xeno.sce$ct[grep('Lymphoid', xeno.sce$ct)] <- 'Lymphoid'
xeno.sce$ct[grep('Motor', xeno.sce$ct)] <- 'Motor neuron'
xeno.sce$ct[grep('Myeloid', xeno.sce$ct)] <- 'Myeloid'

new_ctlist <- list(baron_panc = scelist$baron_panc$label,
                   lawlor_panc = scelist$lawlor_panc$`cell type`,
                   muraro_panc = scelist$muraro_panc$label,
                   chen_brain = scelist$chen_brain$SVM_clusterID,
                   darmanis_brain = scelist$darmanis_brain$cell.type,
                   xenopus = xeno.sce$ct)

ctlist <- append(old_ctlist, new_ctlist)

convert_na <- function(x){
  x[which(is.na(x))] <- 'NA'
  x
}

ctlist <- lapply(ctlist, convert_na)

saveRDS(ctlist, file = paste0(save_dir, 'ctlist.rds'))

dat_names <- intersect(names(ctlist), names(emb_list))
names(dat_names) <- dat_names


# clustering + ARI results


clust_list <- mclapply(as.list(dat_names),
                       FUN = function(nam) try(clusters_acrossPCs(emb_list[[nam]], ctlist[[nam]], dataset_name = nam)),
                       mc.cores = mcc)

saveRDS(clust_list, file = paste0(save_dir, 'nngraph_clust.rds'))

# louvclust_list <- mclapply(as.list(dat_names),
#                        FUN = function(nam) try(clusters_acrossPCs(emb_list[[nam]], ctlist[[nam]], dataset_name = nam, blusparam = NNGraphParam(cluster.fun="louvain"))),
#                        mc.cores = mcc)
# 
# saveRDS(louvclust_list, file = paste0(save_dir, 'louvain_clust.rds'))
