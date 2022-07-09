# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)
#library(fishualize)

# Get variables

## Debug
#CONC = as.list(list.files("/hps/nobackup/birney/users/ian/pilot/hmm_concordance_recode/",full.names = T, recursive = T))
#KW = "/hps/nobackup/birney/users/ian/pilot/kruskal_wallis/out.rds"

## True
CONC = snakemake@input[["conc"]]
KW = snakemake@input[["kw"]]
OUT_PNG = snakemake@output[["png"]]
OUT_PDF = snakemake@output[["pdf"]]

# Read and process

kw = readRDS(KW) %>% 
  dplyr::mutate(INTERVAL = as.numeric(INTERVAL)) %>% 
  dplyr::mutate(N_STATES = as.numeric(N_STATES))

conc = purrr::map_dfr(CONC, readRDS) %>% 
  dplyr::group_by(INTERVAL, N_STATES) %>% 
  # get mean across groups A and B
  dplyr::summarise(MEAN_CONC = mean(CONCORDANCE))

df = dplyr::left_join(kw, conc, by = c("INTERVAL", "N_STATES"))

# Plot

## Get palette
#devtools::source_gist("b68d094672bdde2f0144b755b4721fe6")
#pal = grDevices::colorRampPalette(electro_angler)
pal = colorspace::sequential_hcl(length(unique(df$N_STATES)), palette = "ag_Sunset")

out_plot = df %>% 
  dplyr::mutate(N_STATES = factor(N_STATES, levels = sort(unique(N_STATES)))) %>% 
  dplyr::mutate(INTERVAL = factor(INTERVAL, levels = sort(unique(INTERVAL)))) %>% 
  # make new column combining `INTERVAL` and `N_STATES`
  tidyr::unite(INTERVAL, N_STATES,
               col = "INT_STATE",
               sep = ";",
               remove = F) %>% 
  ggplot(aes(MEAN_CONC, SUM_KW_STAT)) +
    geom_point(aes(size = INTERVAL, colour = N_STATES),
               alpha = 0.8) +
    ggrepel::geom_text_repel(aes(label = INT_STATE),
              size = 2,
              ) +
    theme_bw() +
    guides(size = guide_legend(title = "Interval\n(seconds)"),
           colour = guide_legend(title = "N states")) +
    scale_colour_manual(values = pal) +
    xlab("mean concordance between cross-validated HMM states") +
    ylab("summed Kruskal-Wallis statistic comparing frequency\nof time spent in each HMM state across medaka lines")

#out_plot

ggsave(OUT_PNG,
       out_plot,
       device = "png",
       width = 10,
       height = 8.5,
       units = "in",
       dpi = 400)

ggsave(OUT_PDF,
       out_plot,
       device = "pdf",
       width = 10,
       height = 8.5,
       units = "in",
       dpi = 400)
