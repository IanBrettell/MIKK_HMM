# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = as.list(list.files("/hps/nobackup/birney/users/ian/MIKK_HMM/coverage/hdrr/bwamem2",
                full.names = T))
### Remove KC
to_remove = IN %>% 
  unlist() %>% 
  basename() %>% 
  stringr::str_remove(".txt") %>% 
  grepl("K", .) %>% 
  which()
IN = IN[-to_remove]

## True
IN = snakemake@input
OUT_PNG = snakemake@output[["png"]]
OUT_PDF = snakemake@output[["pdf"]]

# Read in data

names(IN) = IN %>% 
  unlist() %>% 
  basename() %>% 
  stringr::str_remove(".txt") 

df = purrr::map_dfr(IN, readr::read_tsv, .id = "SAMPLE")

# Create palette

pal = fishualize::fish(n = 24, option = "Lepomis_megalotis")

# Generate figure

fig = df %>% 
  dplyr::rename(CHROM = `#rname`) %>% 
  dplyr::filter(CHROM != "MT") %>% 
  dplyr::mutate(CHROM = as.numeric(CHROM),
                CHROM = factor(CHROM, levels = 1:24)) %>% 
  ggplot() +
  geom_histogram(aes(meandepth, fill = CHROM), bins = 40) +
  facet_wrap(~CHROM,nrow = 6, ncol = 4) +
  cowplot::theme_cowplot() +
  scale_fill_manual(values = pal) +
  guides(fill = "none") +
  xlab("mean depth")

# Save

ggsave(OUT_PNG,
       fig,
       device = "png",
       width = 6,
       height = 8,
       units = "in",
       dpi = 400)

ggsave(OUT_PDF,
       fig,
       device = "pdf",
       width = 6,
       height = 8,
       units = "in",
       dpi = 400)
