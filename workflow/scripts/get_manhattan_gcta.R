# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug

IN = "/hps/nobackup/birney/users/ian/somites/gcta/mlma_loco/true/hdrr/None/5000/0.8/unsegmented_psm_area/Microscope.loco.mlma"
MIN_P = "/hps/nobackup/birney/users/ian/somites/gcta/mlma_loco/min_p/hdrr/None/5000/0.8/unsegmented_psm_area/Microscope.csv"
BIN_LENGTH = "5000" %>% 
  as.numeric()
COV = "0.8" %>% 
  as.numeric()
PHENO = "unsegmented_psm_area"
COVARS = "Microscope"

## True

IN = snakemake@input[["res"]]
MIN_P = snakemake@input[["min_p"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
COV = snakemake@params[["cov"]] %>% 
  as.numeric()
PHENO = snakemake@params[["phenotype"]]
COVARS =snakemake@params[["covars"]]
OUT = snakemake@output[["man"]]

########################
# Plotting parameters
########################

gwas_pal = c("#2B2D42", "#F7B267", "#F25C54")
names(gwas_pal) = c("target", "even chr", "odd chr")

# Intercept
intercept_pal = c("#EF476F", "#8D99AE", "#2b2d42")
names(intercept_pal) = c("target", "even chr", "odd chr")

# PSM
unsegmented_psm_area_pal = c("#E59500", "#9C6FC3", "#401F3E")
names(unsegmented_psm_area_pal) = c("target", "even chr", "odd chr")

########################
# HdrR chromosome data
########################
# Get chromosome lengths
med_chr_lens = read.table(here::here("data",
                                     "Oryzias_latipes.ASM223467v1.dna.toplevel.fa_chr_counts.txt"),
                          col.names = c("chr", "end"))
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

df = readr::read_tsv(IN) %>% 
  # Add POS
  dplyr::mutate(BIN_START = (bp * BIN_LENGTH) + 1,
                BIN_END = (bp + 1) * BIN_LENGTH) %>% 
  # join chromosome lengths
  dplyr::left_join(med_chr_lens, by = c("Chr" = "chr")) %>% 
  # add x-coord
  dplyr::mutate(X_COORD = BIN_START + TOT) %>% 
  # change column names
  dplyr::rename(CHROM = Chr)

# Get significance levels

## Permutations

PERM_SIG = readr::read_csv(MIN_P) %>% 
  dplyr::pull(MIN_P) %>% 
  min(.)

## Bonferroni

BONF_SIG = 0.05 / nrow(df)

# Parse covariates

covars = stringr::str_split(COVARS, "-") %>% 
  unlist()

# Set title

TITLE = paste("Phenotype: ",
              PHENO,
              "\nCovariates: ",
              paste(covars, collapse = ", "))

SUBTITLE = paste("Emission covariances: ",
                 COV,
                 "\nBin length: ",
                 BIN_LENGTH)

# Set palette

pal = eval(as.name(paste(PHENO, "_pal", sep = "")))

########################
# Manhattan plot function
########################

plot_man = function(df, title = NULL, subtitle = NULL, gwas_pal, size = 0.5, alpha = 0.5, med_chr_lens, perm_sig = NULL, bonf_sig = NULL){
  # Create palette
  pal = rep_len(gwas_pal, length.out = nrow(med_chr_lens))
  names(pal) = med_chr_lens$chr
  
  df = df %>% 
    # create `COLOUR` vector
    dplyr::mutate(COLOUR = dplyr::case_when(!is.null(perm_sig) & p < perm_sig ~ gwas_pal[1],
                                            gtools::even(CHROM) ~ gwas_pal[2],
                                            gtools::odd(CHROM) ~ gwas_pal[3])) %>% 
    dplyr::mutate(CHROM = factor(CHROM, levels = med_chr_lens$chr)) 
  
  out_plot = df %>% 
    ggplot(aes(x = X_COORD,
               y = -log10(p),
               label = BIN_START,
               label2 = BIN_END)) + 
    geom_point(colour = df$COLOUR,
               size = size,
               alpha = alpha) +
    #scale_color_manual(values = gwas_pal) +
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
    ylab("-log10(p-value)") + 
    # permutations significance level
    geom_hline(yintercept = -log10(perm_sig), colour = "#60D394", linetype = "dashed") +
    geom_text(aes(MID_TOT[1], -log10(perm_sig), label = "permutations", vjust = 1), size = 3, colour = "#60D394") + 
    # bonferroni significance level
    geom_hline(yintercept = -log10(bonf_sig), colour = "#F06449", linetype = "dashed") +
    geom_text(aes(MID_TOT[1], -log10(bonf_sig), label = "bonferroni", vjust = 1), size = 3, colour = "#F06449")
  
  return(out_plot)
  
}

########################
# Plot and save
########################

# Plot
out_plot = plot_man(df,
                    title = TITLE,
                    subtitle = SUBTITLE,
                    gwas_pal = pal,
                    med_chr_lens = med_chr_lens,
                    perm_sig = PERM_SIG,
                    bonf_sig = BONF_SIG)

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
