# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
BED = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds/F2/hdrr/5000/0.8.bed"
BIM = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds/F2/hdrr/5000/0.8.bim"
FAM = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds/F2/hdrr/5000/0.8.fam"
LINE_COLS = here::here("config/line_colours/line_colours_0.08.csv")
PHENOS = as.list(c("/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/10.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/11.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/12.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/13.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/14.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/15.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/1.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/2.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/3.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/4.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/5.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/6.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/7.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/8.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/novel_object/9.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/10.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/11.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/12.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/13.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/14.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/15.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/1.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/2.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/3.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/4.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/5.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/6.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/7.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/8.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/dge/invnorm/open_field/9.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/10.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/11.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/12.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/13.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/14.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/15.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/1.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/2.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/3.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/4.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/5.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/6.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/7.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/8.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/novel_object/9.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/10.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/11.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/12.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/13.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/14.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/15.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/1.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/2.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/3.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/4.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/5.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/6.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/7.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/8.phen",
                   "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.08/15/sge/invnorm/open_field/9.phen"))
SAMPLES = here::here("config/F2_samples_converted.csv")
SIGS = "/hps/nobackup/birney/users/ian/MIKK_HMM/sig_snps/hdrr/dist_angle/0.08/15/5000/0.8/invnorm/sigs.csv"
N_STATES = 15
OUT_DIR = here::here("book/figs/sig_snps_boxplots/dist_angle/0.08/15/invnorm")

##True
BED = snakemake@input[["bed"]]
BIM = snakemake@input[["bim"]]
FAM = snakemake@input[["fam"]]
LINE_COLS = snakemake@input[["line_cols"]]
PHENOS = snakemake@input[["phenos"]]
SAMPLES = snakemake@input[["samples_file"]]
SIGS = snakemake@input[["sigs"]]
N_STATES = snakemake@params[["n_states"]]
OUT_DIR = snakemake@params[["out_dir"]]
GREEDY_OUT = snakemake@output[["greedy"]]

####################

dir.create(OUT_DIR, recursive = T)

names(PHENOS) = purrr::map_chr(PHENOS, function(FILE){
  split_file = FILE %>% 
    stringr::str_split("/", simplify = T) %>% 
    as.vector()
  TOT_LEN = length(split_file)
  DGE_SGE = split_file[TOT_LEN - 3]
  ASSAY = split_file[TOT_LEN - 1]
  STATE = split_file[TOT_LEN] %>% 
    stringr::str_remove(".phen")
  OUT = paste(DGE_SGE, ASSAY, STATE, sep = "__")
  return(OUT)
})
# Read in files

bim = genio::read_bim(BIM)
fam = genio::read_fam(FAM)
bed = genio::read_bed(BED,
                      m_loci = nrow(bim),
                      n_ind = nrow(fam))


samples = readr::read_csv(SAMPLES) %>% 
  # add cross
  dplyr::mutate(CROSS = paste(pat_line, mat_line, sep = "x")) %>% 
  dplyr::select(SAMPLE = finclip_id, CROSS, pat_line, mat_line) %>% 
  dplyr::mutate(SAMPLE = as.character(SAMPLE)) %>% 
  # convert CROSS to factor so it's always in the same order
  dplyr::mutate(CROSS = as.factor(CROSS))

line_cols = readr::read_csv(LINE_COLS)
pal = line_cols$colour ; names(pal) = line_cols$line

phen_df = purrr::map_dfr(PHENOS, ~readr::read_tsv(.x, 
                                                    col_names = c("SAMPLE", "IID", "STATE_FREQ"),
                                                    col_types = c("ccd")),
                           .id = "GE_ASSAY_STATE") %>% 
  tidyr::separate(col = "GE_ASSAY_STATE",
                  into = c("DGE_SGE", "ASSAY", "STATE"),
                  sep = "__") 

sigs = readr::read_csv(SIGS, col_types = c("cciiicccccdddd"))

# Greedily filter significant SNPs

WINDOW_LEN = 1e5

# Function for greedy filtering
greedy_drop = function(df){
  # sort `df` by p-value
  df = df %>% 
    dplyr::arrange(p)
  # add first row (lowest p-value) to output
  out = df[1, ]
  # for all other rows...
  if (nrow(df) > 1){
    for (i in 2:nrow(df)){
      # if its absolute bp distance from any other kept SNP is greater than the window length
      if (min(abs(df[[i,"POS"]] - out$POS)) > WINDOW_LEN){
        # add that row to the output
        out = out %>% 
          tibble::add_row(df[i, ])
      }
    }    
  }
  return(out)
}

# Apply function to each combination
greedy_keep = split(sigs, ~ASSAY + DGE_SGE + STATE + CHROM) %>% 
  purrr::map(~ .x %>% 
               dplyr::arrange(p)) %>% 
  # drop empty data frames
  purrr::discard(~nrow(.x) == 0) %>% 
  purrr::map_dfr(~greedy_drop(.x)) %>% 
  # sort
  dplyr::arrange(CHROM, POS)

# Write greedy-filtered SNPs to file
readr::write_csv(greedy_keep, GREEDY_OUT)

## Create a list of DFs with phenotypes
#ALL_STATES = as.character(1:N_STATES); names(ALL_STATES) = ALL_STATES
#phenos_all = purrr::map(ALL_STATES, function(TARGET_STATE){
#  out = tibble::tibble(SAMPLE = fam$id) %>% 
#    # bind with SAMPLES
#    dplyr::left_join(samples,
#                     by = "SAMPLE") %>% 
#    # join phenos
#    dplyr::left_join(phen_list[[TARGET_STATE]] %>% 
#                       dplyr::select(SAMPLE, SF),
#                     by = "SAMPLE")
#  return(out)
#})

############################
# Loop over each significant SNP and save plot
############################

lapply(1:nrow(greedy_keep), function(INDEX){
  # Get locus for SNP
  #TARGET_SNP = "19:4921761"
  TARGET_SNP = greedy_keep$SNP[INDEX]

  # Get genetic effect
  #TARGET_GE = greedy_keep$DGE_SGE[INDEX]
  
  TARGET_ASSAY = greedy_keep$ASSAY[INDEX]
  # Get state
  #TARGET_STATE = 9
  TARGET_STATE = as.character(greedy_keep$STATE[INDEX])


  # Get index in bim for that SNP
  BIM_IND = which(bim$id == TARGET_SNP)
  
  # Pull out genotypes for that SNP
  GENOS = tibble::tibble(SAMPLE = fam$id,
                         GENOS = bed[BIM_IND, ]) %>% 
    # add cross
    dplyr::left_join(samples,
                     by = "SAMPLE")
  
  
  # Bind fam to genos
  
  df = GENOS %>% 
    # join phenos
    dplyr::left_join(phen_df %>%
                       dplyr::filter(ASSAY == TARGET_ASSAY & STATE == TARGET_STATE) %>% 
                       dplyr::select(SAMPLE, DGE_SGE, STATE_FREQ),
                     by = "SAMPLE") %>%
    # convert GENO to factor
    dplyr::mutate(GENOS = factor(GENOS, levels = 0:2)) %>% 
    # pivot longer
    #tidyr::pivot_longer(cols = c("DGE", "SGE"),names_to = "DGE_SGE",values_to = "STATE_FREQ") %>% 
    # remove NAs
    tidyr::drop_na() 
  
  
  # Plot
  
  fig = df %>% 
    ggplot(aes(GENOS, STATE_FREQ)) +
    gghalves::geom_half_boxplot(aes(fill = pat_line), side = "l") +
    gghalves::geom_half_boxplot(aes(fill = mat_line), side = "r") +
    ggbeeswarm::geom_beeswarm(groupOnX = T) +
    facet_grid(rows = vars(DGE_SGE),
               cols = vars(CROSS)) +
    cowplot::theme_cowplot() +
    ggtitle(paste("SNP: ", TARGET_SNP,
                  "\nAssay: ", stringr::str_replace(TARGET_ASSAY, "_", " "),
                  "\nState: ", TARGET_STATE,
                  sep = "")) +
    xlab("genotype") +
    ylab("inverse-normalised\nstate frequency") +
    guides(fill = "none") +
    #scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal)
  
  # Save
  ## Adjust plot width by number of unique crosses
  PLOT_WID  = 1.5 * length(unique(df$CROSS))
  OUT_PNG = file.path(OUT_DIR,
                      paste(TARGET_STATE, "-", TARGET_SNP, ".png", sep = ""))
  OUT_PDF = file.path(OUT_DIR,
                      paste(TARGET_STATE, "-", TARGET_SNP, ".pdf", sep = ""))
  ggsave(OUT_PNG,
         fig,
         device = "png",
         width = PLOT_WID,
         height = 6,
         dpi = 400)
  
  ggsave(OUT_PDF,
         fig,
         device = "pdf",
         width = PLOT_WID,
         height = 6,
         dpi = 400)
  
  message(paste("Index done:", INDEX))
  #return(fig)
})

