# corral_manuscript

Code to reproduce results and figures in the manuscript describing `corral`.

## How to reproduce figures

### Figures 1, 3, and 4: 

Knit `reproduce_figs/1_4-figures.Rmd`

### Figure 2:

1. Generate the embeddings (run `reproduce_results/01_make_all_transformations_data.R`)
2. Cluster embeddings (run `reproduce_results/02_clust_ARI_transformations_data.R`)
3. Run chunks in `reproduce_figs/2-CAvariations.Rmd`

### Figure 5:

Knit `reproduce_figs/5-time_fig_zm8.Rmd`


### Supplementary figures:

- **S1**: Run chunks in `reproduce_figs/suppl_ftresiduals.Rmd`
- **S2**: Re-run steps from **Figure 2** above, setting the `ncomps` variable to 30 in the two `reproduce_results` scripts
- **S3**: Included in `reproduce_figs/2-CAvariations.Rmd`

