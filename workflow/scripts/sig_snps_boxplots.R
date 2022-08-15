# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

BED = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds/F2/hdrr/5000/0.8.bed"
BIM = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds/F2/hdrr/5000/0.8.bim"
FAM = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds/F2/hdrr/5000/0.8.fam"
LINE_COLS = here::here("config/line_colours/line_colours_0.05.csv")

## NO
PHENOS_DGE = as.list((list.files("/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.05/15/dge/invnorm/novel_object",
                         full.names = T)))
PHENOS_SGE = as.list((list.files("/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.05/15/sge/invnorm/novel_object",
                                 full.names = T)))
SAMPLES = here::here("config/F2_samples_converted.csv")
SIGS = "/hps/nobackup/birney/users/ian/MIKK_HMM/sig_snps/hdrr/dist_angle/dist_angle/15/5000/0.8/invnorm/novel_object_sigs.csv"
N_STATES = 15
OUT_DIR = here::here("book/figs/sig_snps_boxplots/dist_angle/0.05/15/invnorm/novel_object")

## OF
PHENOS_DGE = as.list((list.files("/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.05/15/dge/invnorm/open_field",
                                 full.names = T)))
PHENOS_SGE = as.list((list.files("/hps/nobackup/birney/users/ian/MIKK_HMM/phens_sf/true/dist_angle/0.05/15/sge/invnorm/open_field",
                                 full.names = T)))
SAMPLES = here::here("config/F2_samples_converted.csv")
SIGS = "/hps/nobackup/birney/users/ian/MIKK_HMM/sig_snps/hdrr/dist_angle/dist_angle/15/5000/0.8/invnorm/open_field_sigs.csv"
N_STATES = 15
OUT_DIR = here::here("book/figs/sig_snps_boxplots/dist_angle/0.05/15/invnorm/open_field")




####################

dir.create(OUT_DIR, recursive = T)

names(PHENOS_DGE) = PHENOS_DGE %>% 
  unlist() %>% 
  basename() %>% 
  stringr::str_remove(".phen")

names(PHENOS_SGE) = PHENOS_SGE %>% 
  unlist() %>% 
  basename() %>% 
  stringr::str_remove(".phen")


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

phen_list_dge = purrr::map(PHENOS_DGE, ~readr::read_tsv(.x, 
                                                col_names = c("SAMPLE", "IID", "SF"),
                                                col_types = c("ccd")))

phen_list_sge = purrr::map(PHENOS_SGE, ~readr::read_tsv(.x, 
                                                col_names = c("SAMPLE", "IID", "SF"),
                                                col_types = c("ccd")))

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
  purrr::map_dfr(~greedy_drop(.x))

# Create a list of DFs with phenotypes
ALL_STATES = as.character(1:N_STATES); names(ALL_STATES) = ALL_STATES
phenos_all = purrr::map(ALL_STATES, function(TARGET_STATE){
  out = tibble::tibble(SAMPLE = fam$id) %>% 
    # bind with SAMPLES
    dplyr::left_join(samples,
                     by = "SAMPLE") %>% 
    # join phenos
    dplyr::left_join(phen_list_dge[[TARGET_STATE]] %>% 
                       dplyr::select(SAMPLE, DGE = SF),
                     by = "SAMPLE") %>%
    dplyr::left_join(phen_list_sge[[TARGET_STATE]] %>% 
                       dplyr::select(SAMPLE, SGE = SF),
                     by = "SAMPLE") 
  return(out)
})

############################
# Loop over each significant SNP and save plot
############################

lapply(1:nrow(greedy_keep), function(INDEX){
  # Get locus for SNP
  TARGET_SNP = greedy_keep$SNP[INDEX]
  TARGET_SNP = "19:4921761"
  
  # Get state
  TARGET_STATE = as.character(greedy_keep$STATE[INDEX])
  TARGET_STATE = 9
  
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
    dplyr::left_join(phenos_all[[TARGET_STATE]] %>% 
                       dplyr::select(SAMPLE, DGE, SGE),
                     by = "SAMPLE") %>%
    # convert GENO to factor
    dplyr::mutate(GENOS = factor(GENOS, levels = 0:2)) %>% 
    # pivot longer
    tidyr::pivot_longer(cols = c("DGE", "SGE"),names_to = "DGE_SGE",values_to = "STATE_FREQ") %>% 
    # remove NAs
    tidyr::drop_na() 
  
  
  # Plot
  
  fig = df %>% 
    ggplot(aes(GENOS, STATE_FREQ)) +
    gghalves::geom_half_boxplot(aes(fill = pat_line), side = "l") +
    gghalves::geom_half_boxplot(aes(fill = mat_line), side = "r") +
    ggbeeswarm::geom_beeswarm() +
    facet_grid(rows = vars(DGE_SGE),
               cols = vars(CROSS)) +
    cowplot::theme_cowplot() +
    ggtitle(paste("SNP: ", TARGET_SNP,
                  "\nState: ", TARGET_STATE,
                  sep = "")) +
    xlab("genotype") +
    ylab("inverse-normalised\nstate frequency") +
    guides(colour = "none",
           fill = "none") +
    #scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal)
  
  # Save
  ## Adjust plot width by number of unique crosses
  PLOT_WID  = 2 * length(unique(df$CROSS))
  OUT_PATH = file.path(OUT_DIR,
                       paste(TARGET_STATE, "-", TARGET_SNP, ".png", sep = ""))
  ggsave(OUT_PATH,
         fig,
         device = "png",
         width = PLOT_WID,
         height = 6,
         dpi = 400)
  #return(fig)
})


  
