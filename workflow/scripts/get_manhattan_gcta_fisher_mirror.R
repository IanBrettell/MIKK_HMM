# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug

IN_DGE = "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco/fisher_method/hdrr/dist_angle/0.05/15/5000/0.8/dge/invnorm/open_field/time-quadrant/nyholt.csv"
IN_SGE = "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco/fisher_method/hdrr/dist_angle/0.05/15/5000/0.8/sge/invnorm/open_field/time-quadrant/nyholt.csv"
BIN_LENGTH = "5000" %>% 
  as.numeric()
COV = "0.8" %>% 
  as.numeric()
COVARS = "time-quadrant" %>% 
  stringr::str_replace("-", ", ")
ASSAY = "novel_object"
TRANS = "invnorm"
PHENO = "state_freq"
ADJ_METHOD = "nyholt"
## True

IN_DGE = snakemake@input[["dge"]]
IN_SGE = snakemake@input[["sge"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
COV = snakemake@params[["cov"]] %>% 
  as.numeric()
TRANS = snakemake@params[["transformation"]]
ASSAY = snakemake@params[["assay"]]
COVARS =snakemake@params[["covars"]] %>% 
  stringr::str_replace("-", ", ")
PHENO = snakemake@params[["pheno"]]
ADJ_METHOD = snakemake@params[["adj_method"]]
OUT = snakemake@output[["man"]]
TSV = snakemake@output[["sig"]]


########################
# Plotting parameters
########################

# Palettes 
dge_pal = c("red", "#F5D547", "#FCA17D")
sge_pal = c("red", "#DA627D", "#9A348E")

# Set names
names(dge_pal) = c("target", "even chr", "odd chr")
names(sge_pal) = c("target", "even chr", "odd chr")

if (TRANS == "notrans"){
  TRANS = "None"
} else if (TRANS == "invnorm"){
  TRANS = "inverse-normalised"
}

if (COVARS == "All"){
  COVARS = "date, time, quadrant, tank side"
}

ASSAY = stringr::str_replace(ASSAY, "_", " ")

########################
# HdrR chromosome data
########################
# Get chromosome lengths
med_chr_lens = readr::read_csv(here::here("config/hdrr_chrom_lengths.csv"),
                               col_names = c("chr", "end"))
# Add start
med_chr_lens$start = 1
# Reorder
med_chr_lens = med_chr_lens %>% 
  dplyr::select(chr, start, end) %>% 
  # remove MT
  dplyr::filter(chr != "MT") %>% 
  # convert to integer
  dplyr::mutate(chr = as.integer(chr)) %>% 
  # Add cumulative bases
  dplyr::mutate(CUMSUM = cumsum(end),
                TOT = CUMSUM - end) %>% 
  # Add midpoint for each chr
  dplyr::mutate(MID_TOT = TOT + (end / 2))

########################
# Read in files
########################

# Read in and process data

## DGE

df_dge = readr::read_csv(IN_DGE,
                         col_types = c("icid")) %>% 
  # join chromosome lengths
  dplyr::left_join(med_chr_lens, by = c("Chr" = "chr")) %>% 
  # add x-coord
  dplyr::mutate(X_COORD = bp + TOT) %>% 
  # change column names
  dplyr::rename(CHROM = Chr)

## SGE

df_sge = readr::read_csv(IN_DGE,
                         col_types = c("icid")) %>% 
  # join chromosome lengths
  dplyr::left_join(med_chr_lens, by = c("Chr" = "chr")) %>% 
  # add x-coord
  dplyr::mutate(X_COORD = bp + TOT) %>% 
  # change column names
  dplyr::rename(CHROM = Chr)

# Get significance levels

## Bonferroni

BONF_SIG = 0.05 / nrow(df)

# Get maximum limit for y-axis

Y_MAX = min(BONF_SIG, df_dge$FISHER_P, df_sge$FISHER_P)

# Set title

if (PHENO == "mean_speed"){
  TITLE = paste(DGE_SGE,
                "\nAssay: ",
                ASSAY,
                "\nPhenotype: mean speed",
                sep = ""
                )  
} else if (PHENO == "state_freq"){
  TITLE = paste(DGE_SGE,
                "\nFisher method adjustment method: ",
                ADJ_METHOD,
                "\nAssay: ",
                ASSAY,
                "\nPhenotype: state frequency",
                sep = ""
  )
}


SUBTITLE = paste("Transformation: ",
                 TRANS,
                 "\nCovariates: ",
                 COVARS,
                 "\nEmission covariances: ",
                 COV,
                 "\nBin length: ",
                 BIN_LENGTH,
                 sep = "")

########################
# Manhattan plot function
########################

plot_man = function(df, title = NULL, subtitle = NULL, gwas_pal, size = 0.5, alpha = 0.5, med_chr_lens, perm_sig = NULL, bonf_sig = NULL){
  # Create palette
  pal = rep_len(gwas_pal, length.out = nrow(med_chr_lens))
  names(pal) = med_chr_lens$chr
  
  df = df %>% 
    # create `COLOUR` vector
    dplyr::mutate(COLOUR = dplyr::case_when(FISHER_P < bonf_sig ~ gwas_pal[["target"]],
                                            gtools::even(CHROM) ~ gwas_pal[["even chr"]],
                                            gtools::odd(CHROM) ~ gwas_pal[["odd chr"]])) %>% 
    dplyr::mutate(CHROM = factor(CHROM, levels = med_chr_lens$chr)) 
#  %>% 
#    dplyr::slice_sample(n = 1e5)
  
  out_plot = df %>% 
    ggplot(aes(x = X_COORD,
               y = -log10(FISHER_P),
               label = bp)) + 
    geom_point(colour = df$COLOUR,
               size = size,
               alpha = alpha) +
    # permutations significance level
    #geom_hline(yintercept = -log10(perm_sig), colour = "#60D394", linetype = "dashed") +
    #geom_text(aes(MID_TOT[1], -log10(perm_sig), label = "permutations", vjust = 1), size = 3, colour = "#60D394") + 
    # bonferroni significance level
    geom_hline(yintercept = -log10(bonf_sig), colour = "#F06449", linetype = "dashed") +
    geom_text(aes(MID_TOT[1], -log10(bonf_sig), label = "bonferroni", vjust = 1), size = 3, colour = "#F06449") +
    scale_x_continuous(breaks = med_chr_lens$MID_TOT, 
                       labels = med_chr_lens$chr) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
    ) +
    guides(colour = "none") +
    labs(title = title,
         subtitle = subtitle) +
    xlab("Chromosome") +
    ylab("-log10(p-value)")

  
  return(out_plot)
  
}

########################
# Plot and save
########################

# Plot
dge_plot = plot_man(df_dge %>% 
                      dplyr::slice_sample(n = 1e5),
                    title = TITLE,
                    subtitle = SUBTITLE,
                    gwas_pal = dge_pal,
                    med_chr_lens = med_chr_lens,
                    bonf_sig = BONF_SIG)

sge_plot = plot_man(df_sge,
                    title = TITLE,
                    subtitle = SUBTITLE,
                    gwas_pal = sge_pal,
                    med_chr_lens = med_chr_lens,
                    bonf_sig = BONF_SIG) +
    scale_y_reverse()

# Save
## Make sure the directory exists
dir.create(dirname(OUT), recursive = T, showWarnings = F)
## Write
ggsave(OUT,
       out_plot,
       device = "png",
       width = 9.6,
       height = 6,
       units = "in",
       dpi = 400)


########################
# Pull out significant hits and save
########################

sig_hits = df %>% 
  dplyr::filter(p < PERM_SIG) %>% 
  dplyr::select(-c(start, end, CUMSUM, TOT, MID_TOT, X_COORD)) 

if (nrow(sig_hits) == 0){
  file.create(TSV)
} else{
  readr::write_tsv(sig_hits, TSV)
}

