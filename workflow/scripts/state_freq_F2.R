# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out/0.05/dist_angle/15.csv"
DGE_SGE = "dge"

## True
IN = snakemake@input[["data"]]
OUT_CSV_NT = snakemake@output[["csv_notrans"]]
OUT_CSV_IV = snakemake@output[["csv_invnorm"]]
OUT_HIST_PNG = snakemake@output[["hist_png"]]
OUT_HIST_PDF = snakemake@output[["hist_pdf"]]
DGE_SGE = snakemake@params[["dge_sge"]]

# Get output directories for .csv files

OUT_DIR_NT = dirname(OUT_CSV_NT[[1]])
OUT_DIR_IV = dirname(OUT_CSV_IV[[1]])

# Read data

## NOTE: we're using the full dataset rather than split (between F0, F2 and KCC) 
## so that the recoded states are the same across those datasets (as they're
## sorted by mean speed, which may differ across datasets)

raw = readr::read_csv(IN) %>% 
  # convert `time` to character and add a 0 if only 3 characters
  dplyr::mutate(time = as.character(time),
                time = dplyr::if_else(nchar(time) == 3,
                                      paste("0", time, sep = ""),
                                      time))

# Add inverse-normalisation function
invnorm = function(x) {
  res = rank(x)
  res = qnorm(res/(length(res)+0.5))
  return(res)
}

# Clean

df = raw %>% 
  # Get individual
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv", 
               remove = F) %>% 
  # add `line` %>% 
  dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                        fish == "test" ~ test_fish)) %>% 
  # recode and order `assay` 
  dplyr::mutate(assay = stringr::str_replace(assay, pattern = "_", " "),
                assay = factor(assay, levels = c("open field", "novel object"))) %>% 
  # convert `date` to factor
  dplyr::mutate(date = factor(date))

# Recode states by mean distance

rank_df = df %>% 
  dplyr::group_by(state) %>% 
  dplyr::summarise(mean_dist = mean(distance)) %>% 
  # rank
  dplyr::arrange(mean_dist) %>% 
  dplyr::mutate(rank = 1:nrow(.))

recode_vec = rank_df %>% 
  dplyr::pull(rank)
names(recode_vec) = rank_df %>% 
  dplyr::pull(state)

# Recode `state`

df = df %>% 
  dplyr::mutate(state_recode = dplyr::recode(state, !!!recode_vec),
                state_recode = factor(state_recode, levels = recode_vec))


# Filter for F0
df = df %>% 
  dplyr::filter(dataset == "F2")

#######################
## DGE
#######################

if (DGE_SGE == "dge"){
  # Get proportion of time each fish spent in each state
  df_dge = df %>% 
    # remove iCab 
    dplyr::filter(fish == "test") %>% 
    ## count rows per fish per state. `.drop` added to include the states with 0 counts
    dplyr::count(indiv, assay, state_recode, .drop = F) %>% 
    # add total row count per fish
    dplyr::add_count(indiv, assay, wt = n, name = "nn") %>% 
    # get proportion of time fish spent in each state
    dplyr::mutate(state_freq = n / nn) %>% 
    # drop NAs (created by n == 0 and nn == 0, where the fish is not tracked at all)
    tidyr::drop_na()
  
  # Split into states
  df_list = df_dge %>% 
    split(., f = .$state_recode)
  # Write to files
  counter = 0
  lapply(df_list, function(DF){
    counter <<- counter + 1
    OUT_PATH = file.path(OUT_DIR_NT, paste(names(df_list)[counter],
                                           ".csv",
                                           sep = ""))
    readr::write_csv(DF, OUT_PATH)
  })
  
  # Split by assay
  
  dge_hist_pretrans = df_dge %>% 
    ggplot() + 
    geom_histogram(aes(state_freq, fill = state_recode),
                   bins = 40) +
    facet_grid(rows = vars(state_recode),
               cols = vars(assay)) +
    cowplot::theme_cowplot() +
    scale_fill_viridis_d() +
    guides(fill = "none")+
    xlab("state frequency")
  
  ### Inverse-normalise
  
  df_dge = df_dge %>% 
    dplyr::group_by(assay, state_recode) %>% 
    dplyr::mutate(state_freq = invnorm(state_freq)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(indiv, assay, state_recode)
  
  # Split by states
  df_list = df_dge %>% 
    split(., f = .$state_recode)
  # Write to files
  counter = 0
  lapply(df_list, function(DF){
    counter <<- counter + 1
    OUT_PATH = file.path(OUT_DIR_IV, paste(names(df_list)[counter],
                                           ".csv",
                                           sep = ""))
    readr::write_csv(DF, OUT_PATH)
  })
  
  # Plot
  
  dge_hist_posttrans = df_dge %>% 
    ggplot() + 
    geom_histogram(aes(state_freq, fill = state_recode),
                   bins = 40) +
    facet_grid(rows = vars(state_recode),
               cols = vars(assay)) +
    cowplot::theme_cowplot() +
    scale_fill_viridis_d() +
    guides(fill = "none") +
    xlab("state frequency (inverse-normalised per state/assay)")
  
  # Compile into single plot
  
  dge_hist = cowplot::plot_grid(dge_hist_pretrans,
                                dge_hist_posttrans,
                                align = "hv",axis = "tblr")
  
  ggsave(OUT_HIST_PNG,
         dge_hist,
         device = "png",
         width = 9,
         height = 11,
         units = "in",
         dpi = 400)
  
  ggsave(OUT_HIST_PDF,
         dge_hist,
         device = "pdf",
         width = 9,
         height = 11,
         units = "in",
         dpi = 400)
}

#######################
## SGE
#######################

if (DGE_SGE == "sge"){
  # Get proportion of time each fish spent in each state
  df_sge = df %>% 
    # remove iCab 
    dplyr::filter(fish == "ref") %>% 
    ## count rows per fish per state. `.drop` added to include the states with 0 counts
    dplyr::count(indiv, assay, state_recode, .drop = F) %>% 
    # add total row count per fish
    dplyr::add_count(indiv, assay, wt = n, name = "nn") %>% 
    # get proportion of time fish spent in each state
    dplyr::mutate(state_freq = n / nn) %>% 
    # drop NAs (created by n == 0 and nn == 0, where the fish is not tracked at all)
    tidyr::drop_na()
  
  # Split by states
  df_list = df_sge %>% 
    split(., f = .$state_recode)
  # Write to files
  counter = 0
  lapply(df_list, function(DF){
    counter <<- counter + 1
    OUT_PATH = file.path(OUT_DIR_NT, paste(names(df_list)[counter],
                                           ".csv",
                                           sep = ""))
    readr::write_csv(DF, OUT_PATH)
  })
  
  sge_hist_pretrans = df_sge %>% 
    ggplot() + 
    geom_histogram(aes(state_freq, fill = state_recode),
                   bins = 40) +
    facet_grid(rows = vars(state_recode),
               cols = vars(assay)) +
    cowplot::theme_cowplot() +
    scale_fill_viridis_d(option = "inferno") +
    guides(fill = "none") +
    xlab("state frequency")
  
  
  ### Inverse-normalise
  
  df_sge = df_sge %>% 
    dplyr::group_by(assay, state_recode) %>% 
    dplyr::mutate(state_freq = invnorm(state_freq)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(indiv, assay, state_recode)
  
  # Split by states
  df_list = df_sge %>% 
    split(., f = .$state_recode)
  # Write to files
  counter = 0
  lapply(df_list, function(DF){
    counter <<- counter + 1
    OUT_PATH = file.path(OUT_DIR_IV, paste(names(df_list)[counter],
                                           ".csv",
                                           sep = ""))
    readr::write_csv(DF, OUT_PATH)
  })
  
  sge_hist_posttrans = df_sge %>% 
    ggplot() + 
    geom_histogram(aes(state_freq, fill = state_recode),
                   bins = 40) +
    facet_grid(rows = vars(state_recode),
               cols = vars(assay)) +
    cowplot::theme_cowplot() +
    scale_fill_viridis_d(option = "inferno") +
    guides(fill = "none") +
    xlab("state frequency (inverse-normalised per state/assay)")
  
  # Compile into single plot
  
  sge_hist = cowplot::plot_grid(sge_hist_pretrans,
                                sge_hist_posttrans,
                                align = "hv",axis = "tblr")
  
  ggsave(OUT_HIST_PNG,
         sge_hist,
         device = "png",
         width = 9,
         height = 11,
         units = "in",
         dpi = 400)
  
  ggsave(OUT_HIST_PDF,
         sge_hist,
         device = "pdf",
         width = 9,
         height = 11,
         units = "in",
         dpi = 400)
}



