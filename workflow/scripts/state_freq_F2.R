# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out/0.05/dist_angle/15.csv"
LINE_COLS = here::here("config/line_colours/line_colours_0.05.csv")
N_STATES = 15
DGE_SGE = "dge"

## True
IN = snakemake@input[["data"]]
LINE_COLS = snakemake@input[["line_cols"]]
N_STATES = snakemake@params[["n_states"]]
OUT_CSV_NT = snakemake@output[["csv_notrans"]]
OUT_CSV_IV = snakemake@output[["csv_invnorm"]]
OUT_HIST = snakemake@output[["hist"]]


# Deauthorise google sheets so that it doesn't ask for prompt
googlesheets4::gs4_deauth()

# Read 

## NOTE: we're using the full dataset rather than split (between F0, F2 and KCC) 
## so that the recoded states are the same across those datasets (as they're
## sorted by mean speed, which may differ across datasets)

raw = readr::read_csv(IN)

line_cols = readr::read_csv(LINE_COLS)

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
  # rename outbred
  dplyr::mutate(line = dplyr::if_else(stringr::str_detect(line,
                                                          "outbred"),
                                      "outbred",
                                      line)) %>% 
  # order `line` by mean speed
  dplyr::mutate(line = factor(line, levels = line_cols$line)) %>% 
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

if (DGE_SGE = "dge"){
  # Get proportion of time each fish spent in each state
  df_dge = df %>% 
    # remove iCab 
    dplyr::filter(test_fish != "iCab") %>% 
    ## count rows per fish per state
    dplyr::count(indiv, assay, line, date, time, quadrant, tank_side, state_recode) %>% 
    # add total row count per fish
    dplyr::add_count(indiv, assay, line, date, time, quadrant, tank_side, wt = n, name = "nn") %>% 
    # get proportion of time fish spent in each state
    dplyr::mutate(state_freq = n / nn)
  
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
  
  # Add function
  invnorm = function(x) {
    res = rank(x)
    res = qnorm(res/(length(res)+0.5))
    return(res)
  }
  
  df_dge = df_dge %>% 
    dplyr::group_by(assay, state_recode) %>% 
    dplyr::mutate(state_freq_invnorm = invnorm(state_freq)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(indiv, assay, line, date, time, quadrant, tank_side, state_recode)
  
  dge_hist_posttrans = df_dge %>% 
    ggplot() + 
    geom_histogram(aes(state_freq_invnorm, fill = state_recode),
                   bins = 40) +
    facet_grid(rows = vars(state_recode),
               cols = vars(assay)) +
    cowplot::theme_cowplot() +
    scale_fill_viridis_d() +
    guides(fill = "none") +
    xlab("state frequency (inverse-normalised per state)")
  
  # Compile into single plot
  
  dge_hist = cowplot::plot_grid(dge_hist_pretrans,
                                dge_hist_posttrans,
                                align = "hv",axis = "tblr")
  
  ggsave(DGE_HIST,
         dge_hist,
         device = "png",
         width = 9,
         height = 11,
         units = "in",
         dpi = 400)
}

#######################
## SGE
#######################

if (DGE_SGE = "sge"){
  # Get proportion of time each fish spent in each state
  df_sge = df %>% 
    # take all iCab fishes
    dplyr::filter(line == "iCab") %>% 
    ## count rows per fish per state
    dplyr::count(indiv, assay, test_fish, date, time, quadrant, tank_side, state_recode) %>% 
    # add total row count per fish
    dplyr::add_count(indiv, assay, test_fish, date, time, quadrant, tank_side, wt = n, name = "nn") %>% 
    # get proportion of time fish spent in each state
    dplyr::mutate(state_freq = n / nn)
  
  # Split by assay
  
  sge_hist_pretrans = df_sge %>% 
    ggplot() + 
    geom_histogram(aes(state_freq, fill = state_recode),
                   bins = 40) +
    facet_grid(rows = vars(state_recode),
               cols = vars(assay)) +
    theme_bw() +
    scale_fill_viridis_d(option = "inferno") +
    guides(fill = "none") +
    xlab("state frequency")
  
  
  
  ### Inverse-normalise
  
  
  df_sge = df_sge %>% 
    dplyr::group_by(assay, state_recode) %>% 
    dplyr::mutate(state_freq_invnorm = invnorm(state_freq)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(indiv, assay, test_fish, date, time, quadrant, tank_side, state_recode)
  
  sge_hist_posttrans = df_sge %>% 
    ggplot() + 
    geom_histogram(aes(state_freq_invnorm, fill = state_recode),
                   bins = 40) +
    facet_grid(rows = vars(state_recode),
               cols = vars(assay)) +
    theme_bw() +
    scale_fill_viridis_d(option = "inferno") +
    guides(fill = "none") +
    xlab("state frequency (inverse-normalised per state)")
  
  # Compile into single plot
  
  sge_hist = cowplot::plot_grid(sge_hist_pretrans,
                                sge_hist_posttrans,
                                align = "hv",axis = "tblr")
  
  ggsave(SGE_HIST,
         sge_hist,
         device = "png",
         width = 9,
         height = 11,
         units = "in",
         dpi = 400)  
}



