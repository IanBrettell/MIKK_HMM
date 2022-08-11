# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
### State freq
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco_state_freq_consol/true/hdrr/dist_angle/0.05/15/5000/0.8/dge/invnorm/novel_object/1/time-quadrant.mlma"
MIN_P = "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco/min_p/hdrr/dist_angle/0.05/15/5000/0.8/dge/invnorm/novel_object/1/time-quadrant.csv"
BIN_LENGTH = "5000" %>% 
  as.numeric()
COV = "0.8" %>% 
  as.numeric()
COVARS = "time-quadrant" %>% 
  stringr::str_replace("-", ", ")
STATE = "1" %>% 
  as.numeric()
DGE_SGE = "dge" %>% 
  toupper()
ASSAY = "novel_object"
TRANS = "invnorm"
PHENO = "state_freq"
### Mean speed
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco_mean_speed_consol/true/hdrr/0.08/5000/0.8/sge/notrans/open_field/time-quadrant.mlma"
MIN_P = "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco/min_p/hdrr/0.08/5000/0.8/sge/notrans/open_field/time-quadrant.csv"
BIN_LENGTH = "5000" %>% 
  as.numeric()
COV = "0.8" %>% 
  as.numeric()
COVARS = "time-quadrant" %>% 
  stringr::str_replace("-", ", ")
DGE_SGE = "sge" %>% 
  toupper()
ASSAY = "open_field"
TRANS = "notrans"
PHENO = "mean_speed"

## True

IN = snakemake@input[["res"]]
MIN_P = snakemake@input[["min_p"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
COV = snakemake@params[["cov"]] %>% 
  as.numeric()
DGE_SGE = snakemake@params[["dge_sge"]] %>% 
  toupper()
TRANS = snakemake@params[["transformation"]]
ASSAY = snakemake@params[["assay"]]
STATE = snakemake@params[["state"]] %>% 
  as.numeric()
COVARS =snakemake@params[["covars"]] %>% 
  stringr::str_replace("-", ", ")
PHENO = snakemake@params[["pheno"]]
OUT = snakemake@output[["man"]]
TSV = snakemake@output[["sig"]]


########################
# Plotting parameters
########################

# Palettes for mean speed
if (PHENO == "mean_speed" & DGE_SGE == "DGE"){
  gwas_pal = c("red", "#A4036F", "#EE851B")
} else if (PHENO == "mean_speed" & DGE_SGE == "SGE"){
  gwas_pal = c("red", "#55DDE0", "#33658A")
  # Others for state frequency
} else if (PHENO == "state_freq" & DGE_SGE == "DGE"){
  # get palette
  col = viridis::viridis(n = 15)
  # get complementary colours
  gwas_pal = c("red", col[STATE], colortools::analogous(col[STATE])[2])
} else if (PHENO == "state_freq" & DGE_SGE == "SGE"){
  # get palette
  col = viridis::inferno(n = 15)
  # get complementary colours
  gwas_pal = c("red", col[STATE], colortools::analogous(col[STATE])[2])
  # lighten lower states because they're too dark
  if (STATE < 3){
    gwas_pal[3] = karyoploteR::lighter(gwas_pal[3], amount = 70)
  }
}

# Set names
names(gwas_pal) = c("target", "even chr", "odd chr")

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

print("chrom lengths set")

########################
# Read in files
########################

# Read in and process data

df = readr::read_tsv(IN,
                     col_types = c("iciccdddd")) %>% 
  # join chromosome lengths
  dplyr::left_join(med_chr_lens, by = c("Chr" = "chr")) %>% 
  # add x-coord
  dplyr::mutate(X_COORD = bp + TOT) %>% 
  # change column names
  dplyr::rename(CHROM = Chr)

print("mlma input read")

# Get significance levels

## Permutations

PERM_SIG = readr::read_csv(MIN_P) %>% 
  dplyr::pull(MIN_P) %>% 
  min(.)

## Bonferroni

BONF_SIG = 0.05 / nrow(df)

print("sig levels set")

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
                "\nAssay: ",
                ASSAY,
                "\nPhenotype: state frequency",
                "\nState: ",
                STATE,
                sep = ""
  )
}

print("title set")

SUBTITLE = paste("Transformation: ",
                 TRANS,
                 "\nCovariates: ",
                 COVARS,
                 "\nEmission covariances: ",
                 COV,
                 "\nBin length: ",
                 BIN_LENGTH,
                 sep = "")

print("subtitle set")

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
#  %>% 
#    dplyr::slice_sample(n = 1e5)
  
  out_plot = df %>% 
    ggplot(aes(x = X_COORD,
               y = -log10(p),
               label = bp)) + 
    geom_point(colour = df$COLOUR,
               size = size,
               alpha = alpha) +
    # permutations significance level
    geom_hline(yintercept = -log10(perm_sig), colour = "#60D394", linetype = "dashed") +
    geom_text(aes(MID_TOT[1], -log10(perm_sig), label = "permutations", vjust = 1), size = 3, colour = "#60D394") + 
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
out_plot = plot_man(df,
                    title = TITLE,
                    subtitle = SUBTITLE,
                    gwas_pal = gwas_pal,
                    med_chr_lens = med_chr_lens,
                    perm_sig = PERM_SIG,
                    bonf_sig = BONF_SIG)

print("manhattan plot generated")

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

print("plot saved")
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

