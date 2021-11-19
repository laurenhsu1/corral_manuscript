mthkeep <- c('std', 'varstab', 'anscombe_std', 'ft_x_std', 'ft_p_std',
             'ftalt_p', 'adj_pearson', 'deflate')

mthnames <- list('std' = 'CA (standard)',
                 'varstab' = 'CA on sqrt counts',
                 'anscombe_std' = 'CA on Anscombe counts',
                 'ft_x_std' = 'CA on Freeman Tukey counts',
                 'ft_p_std' = 'CA on Freeman Tukey abundances',
                 'ftalt_p' = 'Freeman Tukey residuals',
                 'adj_pearson' = 'Adjusted Pearson residuals',
                 'deflate' = 'Power deflation (a=.9)', 
                 'glmpca' = 'glmpca',
                 'pca' = 'PCA on logcounts',
                 'deflate_95' = 'Power deflation (a=.95)',
                 'deflate_98' = 'Power deflation (a=.98)',
                 'sqrt_ftres' = 'CA with FT residuals on sqrt counts')

seeds <- list(1234, 12345, 987,9876,98765, 2020, 2021, 123, 31415)
glmpca_nameslist <- as.list(paste0('glmpca.',seeds))
names(glmpca_nameslist) <- paste0('glmpca_seed',seeds)

mthnames <- append(mthnames, glmpca_nameslist)

mthcolors <- c('std' = '#DEA0FD',
               'varstab' = '#FEAF16',
               'anscombe_std' = '#D55E00',
               'ft_x_std' = '#90AD1C',
               'ft_p_std' = '#AAF400',
               'ftalt_p' = '#1D6914',
               'adj_pearson' = '#782AB6',
               'deflate' = '#64ACA7', 
               'glmpca' = '#3283FE',
               'pca' = 'gray', #000000
               'deflate_95' = '#7ED7D1',
               'deflate_98' = '#BEEBE7',
               'CA' = '#DEA0FD',
               'sqrt_ftres' = '#5D3FD3')

namedmthcolors <- mthcolors
names(namedmthcolors) <- mthnames[names(mthcolors)]


# originally from the transformations_fig files (4_...)
umap_plotpanel <- function(x, cts, mthname, ell = FALSE){
  if(length(unique(cts)) > 8){
    if(ell){
      ggobj <- plot_embedding(embedding = x, color_vec = cts, returngg = T, plot_title = mthname, dimname = 'umap', showplot = F, color_pal_vec = pals::alphabet2(),
                              ellipse_vec = cts)
    }
    else{
      ggobj <- plot_embedding(embedding = x, color_vec = cts, returngg = T, plot_title = mthname, dimname = 'umap', showplot = F, color_pal_vec = pals::alphabet2())
    }
  }
  else{
    if(ell){
      ggobj <- plot_embedding(embedding = x, color_vec = cts, returngg = T, plot_title = mthname, dimname = 'umap', showplot = F,
                              ellipse_vec = cts)
    }
    else{
      ggobj <- plot_embedding(embedding = x, color_vec = cts, returngg = T, plot_title = mthname, dimname = 'umap', showplot = F)
    }
  }
  ggobj <- ggobj +
    theme(legend.key.height = unit(.1, 'cm'),
          legend.key.width = unit(.1, 'cm'),
          legend.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))
  return(ggobj)
}

umap_plotlist <- function(umap_list, subtitle_text = NULL, ctl, ...){
  outlist <- list()
  for(dat in names(umap_list)){
    print(dat)
    ggobjs <- lapply(as.list(names(umap_list[[dat]])), FUN = function(x) umap_plotpanel(umap_list[[dat]][[x]], cts = ctl[[dat]], mthname = mthnames[[x]], ...))
    tiled <- cowplot::plot_grid(plotlist = ggobjs, ncol = 3)
    
    title_gg <- cowplot::ggdraw() + labs(title = dat, subtitle = subtitle_text) + theme_classic()
    
    final <- cowplot::plot_grid(title_gg, tiled, ncol = 1, rel_heights = c(0.05, 1))
    outlist[[dat]] <- final
  }
  return(outlist)
}


plot_datfacetARI <- function(ARI_df, plot_title = '', xvar = 'npc', yvar = 'ARI', colvar = 'method', facetby = 'dataset', dat_cols = mthcolors, vline = TRUE){
  max_npc <- max(ARI_df[xvar])
  if(vline) {vlines <- seq(10, max_npc, 10)}
  else {vlines <- c()}
  ARI_df <- ARI_df[order(ARI_df[,xvar]),]
  ggobj <- ggplot(ARI_df, aes_string(x = xvar, y = yvar, color = colvar)) + 
    geom_line() + theme_classic() + ggtitle(plot_title) + 
    facet_wrap(as.formula(paste("~", facetby))) + 
    scale_color_manual(values = dat_cols[unique(ARI_df$method)], 
                       labels = mthnames[unique(ARI_df$method)]) +
    geom_vline(xintercept = vlines, 
               size = .1, color = 'lightgray') + 
    geom_hline(yintercept = c(.25, .5, .75), 
               size = .1, color = 'lightgray')
  return(ggobj)
}

savegg2pdf <- function(ggobj, nam, wdim = 12, hdim = 10, dir = plot_dir){
  cowplot::ggsave2(filename = paste0(dir, nam, '.pdf'),
                   plot = ggobj, 
                   units = 'in',
                   width = wdim,
                   height = hdim)
}
