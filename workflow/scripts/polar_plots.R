# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(wesanderson)
#library(ggbeeswarm)
library(rstatix)
#library(ggpubr)
library(cowplot)

# Get variables

## Debug
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out/0.05/dist_angle/16.csv"
N_STATES = 16
VARIABLES = "distance and angle of travel"
INTERVAL = 0.05
POLAR_ALL_DGE = here::here("book/figs/polar_plots/0.05/dist_angle/16/polar_all_dge.png")
POLAR_ALL_SGE = here::here("book/figs/polar_plots/0.05/dist_angle/16/polar_all_sge.png")

# Deauthorise google sheets so that it doesn't ask for prompt
googlesheets4::gs4_deauth()

## True
IN = snakemake@input[[1]]
N_STATES = snakemake@params[["n_states"]] %>% 
  as.numeric()
INTERVAL = snakemake@params[["interval"]] %>% 
  as.numeric()
VARIABLES = "distance and angle of travel"
POLAR_ALL_DGE = snakemake@output[["polar_all_dge"]]
POLAR_ALL_SGE = snakemake@output[["polar_all_sge"]]
POLAR_BOX_DGE = snakemake@output[["polar_box_dge"]]
POLAR_BOX_SGE = snakemake@output[["polar_box_sge"]]
POLAR_BOX_DGE_SGE = snakemake@output[["polar_box_dge_sge"]]
POLAR_ALL_DGE_SIG_OF = snakemake@output[["polar_all_dge_sig_of"]]
POLAR_ALL_DGE_SIG_NO = snakemake@output[["polar_all_dge_sig_no"]]
POLAR_ALL_SGE_SIG_OF = snakemake@output[["polar_all_sge_sig_of"]]
POLAR_ALL_SGE_SIG_NO = snakemake@output[["polar_all_sge_sig_no"]]
TILE_DGE = snakemake@output[["tile_dge"]]
TILE_SGE = snakemake@output[["tile_sge"]]
TILE_DGE_SGE = snakemake@output[["tile_dge_sge"]]

#######################
# Read in data
#######################

# Get number of rows (for plotting) based on number of states

if (N_STATES == 15 | N_STATES == 20){
  N_ROWS = 5
} else if (N_STATES == 12 | N_STATES == 16) {
  N_ROWS = 4
} else if (N_STATES == 17 | N_STATES == 18){
  N_ROWS = 6
} else if (N_STATES == 14 | N_STATES == 10){
  N_ROWS = 2
} else {
  N_ROWS = 1
}

# And dimensions of `polar_all...` figs
if (N_STATES == 14){
  POL_ALL_WID = 12
  POL_ALL_HEI = 4
} else {
  POL_ALL_WID = 7.5
  POL_ALL_HEI = 10
}

# Get figure height

HEIGHT = 2 * N_ROWS

# Read in file

df = readr::read_csv(IN) %>% 
  # recode angle to sit between 0 and 360
  dplyr::mutate(angle_recode = ifelse(angle < 0,
                                      180 + (180 + angle),
                                      angle))


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
  dplyr::mutate(state_recode = dplyr::recode(state, !!!recode_vec))

#######################
# Polar plots
#######################

polar_all_dge = df %>% 
  # select random sample of 1e5 rows
  dplyr::slice_sample(n = 1e5) %>% 
  # factorise `state_recode`
  #dplyr::mutate(state_recode = factor(state_recode, levels = recode_vec)) %>% 
  ggplot() +
  geom_point(aes(angle_recode, log10(distance), colour = state_recode),
             alpha = 0.3, size = 0.2) +
  coord_polar() +
  facet_wrap(~state_recode, nrow = N_ROWS) +
  #theme_dark(base_size = 8) +
  cowplot::theme_cowplot(font_size = 9) +
  scale_x_continuous(labels = c(0, 90, 180, 270),
                     breaks = c(0, 90, 180, 270)) +
  scale_color_viridis_c() +
  guides(colour = "none") +
  xlab("angle of travel") +
  ylab(expression(log[10]("distance travelled in pixels"))) 
  #ggtitle(TITLE)

ggsave(POLAR_ALL_DGE,
       polar_all_dge,
       device = "png",
       width = POL_ALL_WID,
       height = POL_ALL_HEI,
       units = "in",
       dpi = 400)

#########################

polar_all_sge = df %>% 
  # select random sample of 1e5 rows
  dplyr::slice_sample(n = 1e5) %>% 
  # factorise `state_recode`
  #dplyr::mutate(state_recode = factor(state_recode, levels = recode_vec)) %>% 
  ggplot() +
  geom_point(aes(angle_recode, log10(distance), colour = state_recode),
             alpha = 0.3, size = 0.2) +
  coord_polar() +
  facet_wrap(~state_recode, nrow = N_ROWS) +
  #theme_dark(base_size = 8) +
  cowplot::theme_cowplot(font_size = 9) +
  scale_x_continuous(labels = c(0, 90, 180, 270),
                     breaks = c(0, 90, 180, 270)) +
  scale_color_viridis_c(option = "inferno") +
  guides(colour = "none") +
  xlab("angle of travel") +
  ylab(expression(log[10]("distance travelled in pixels")))

ggsave(POLAR_ALL_SGE,
       polar_all_sge,
       device = "png",
       width = POL_ALL_WID,
       height = POL_ALL_HEI,
       units = "in",
       dpi = 400)

