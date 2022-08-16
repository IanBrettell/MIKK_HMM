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
### 0.05 interval
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out/0.05/dist_angle/15.csv"
N_STATES = 15
VARIABLES = "distance and angle of travel"
INTERVAL = 0.05
SHEET_ID = "118HNbI7ch_0rgRSaDB1b73mpRo5Z2g7uQvNhjZUD7WI"
POLAR_ALL_DGE_SIG_OF = here::here("book/figs/polar_plots_with_sig/0.05/dist_angle/15/polar_all_dge_of.png")
POLAR_ALL_DGE_SIG_NO = here::here("book/figs/polar_plots_with_sig/0.05/dist_angle/15/polar_all_dge_no.png")
POLAR_ALL_SGE_SIG_OF = here::here("book/figs/polar_plots_with_sig/0.05/dist_angle/15/polar_all_sge_of.png")
POLAR_ALL_SGE_SIG_NO = here::here("book/figs/polar_plots_with_sig/0.05/dist_angle/15/polar_all_sge_no.png")
### 0.08 interval
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out/0.08/dist_angle/15.csv"
N_STATES = 15
VARIABLES = "distance and angle of travel"
INTERVAL = 0.08
SHEET_ID = "1kD2ndgmTvGlQw8YZFMztO3YoYfOdSw0SI5_BtomDuaI"
POLAR_ALL_DGE_SIG_OF = here::here("book/figs/polar_plots_with_sig/0.08/dist_angle/15/polar_all_dge_of.png")
POLAR_ALL_DGE_SIG_NO = here::here("book/figs/polar_plots_with_sig/0.08/dist_angle/15/polar_all_dge_no.png")
POLAR_ALL_SGE_SIG_OF = here::here("book/figs/polar_plots_with_sig/0.08/dist_angle/15/polar_all_sge_of.png")
POLAR_ALL_SGE_SIG_NO = here::here("book/figs/polar_plots_with_sig/0.08/dist_angle/15/polar_all_sge_no.png")

# Create output directory
dir.create(dirname(POLAR_ALL_DGE_SIG_OF), recursive = T)
# Deauthorise google sheets so that it doesn't ask for prompt
googlesheets4::gs4_deauth()

## True
IN = snakemake@input[[1]]
N_STATES = snakemake@params[["n_states"]] %>% 
  as.numeric()
INTERVAL = snakemake@params[["interval"]] %>% 
  as.numeric()
VARIABLES = "distance and angle of travel"
SHEET_ID = snakemake@params[["sheet_id"]]
POLAR_ALL_DGE = snakemake@output[["polar_all_dge"]]
POLAR_ALL_SGE = snakemake@output[["polar_all_sge"]]

#######################
# Read in data
#######################

# Get number of rows (for plotting) based on number of states

if (N_STATES == 15){
  N_ROWS = 3
} else if (N_STATES == 20) {
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
  POL_ALL_WID = 8
  POL_ALL_HEI = 5.3
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
# Plot DGE
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
  cowplot::theme_cowplot(font_size = 8) +
  scale_x_continuous(labels = c(0, 90, 180, 270),
                     breaks = c(0, 90, 180, 270)) +
  scale_color_viridis_c() +
  guides(colour = "none") +
  xlab("angle of travel") +
  ylab(expression(log[10]("distance travelled in pixels"))) 
  #ggtitle(TITLE)

# With significant states - OPEN FIELD: NOTE - this doesn't work from within the script

## Get significant states
sig_dge_of = googlesheets4::read_sheet(
  ss = SHEET_ID,
  sheet = "DGE_OF") %>% 
  dplyr::filter(Variable == "line" & `p-value FDR-adj` < 0.05) %>% 
  dplyr::pull(State)

g <- ggplot_gtable(ggplot_build(polar_all_dge))
strips <- which(startsWith(g$layout$name,'strip'))

strip_vec = ifelse(1:N_STATES %in% sig_dge_of,
                   "#ED474A",
                   "#9EA7A9") 
strip_vec = c(strip_vec[11:N_STATES], strip_vec[6:10], strip_vec[1:5])

for (s in seq_along(strips)) {
  g$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- strip_vec[s]
}

grid::grid.newpage()

png(POLAR_ALL_DGE_SIG_OF,
    width = POL_ALL_WID, 
    height = POL_ALL_HEI, 
    units = "in", 
    res = 400)
grid::grid.draw(g)
dev.off()

# With significant states - NOVEL_OBJECT: NOTE - this doesn't work from within the script

## Get significant states
sig_dge_no = googlesheets4::read_sheet(
  ss = SHEET_ID,
  sheet = "DGE_NO") %>% 
  dplyr::filter(Variable == "line" & `p-value FDR-adj` < 0.05) %>% 
  dplyr::pull(State)

g <- ggplot_gtable(ggplot_build(polar_all_dge))
strips <- which(startsWith(g$layout$name,'strip'))

strip_vec = ifelse(1:N_STATES %in% sig_dge_no,
                   "#ED474A",
                   "#9EA7A9") 
strip_vec = c(strip_vec[11:N_STATES], strip_vec[6:10], strip_vec[1:5])

for (s in seq_along(strips)) {
  g$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- strip_vec[s]
}

grid::grid.newpage()

png(POLAR_ALL_DGE_SIG_NO,
    width = POL_ALL_WID, 
    height = POL_ALL_HEI, 
    units = "in", 
    res = 400)
grid::grid.draw(g)
dev.off()


#######################
# Plot SGE
#######################

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
  cowplot::theme_cowplot(font_size = 8) +
  scale_x_continuous(labels = c(0, 90, 180, 270),
                     breaks = c(0, 90, 180, 270)) +
  scale_color_viridis_c(option = "inferno") +
  guides(colour = "none") +
  xlab("angle of travel") +
  ylab(expression(log[10]("distance travelled in pixels")))

# With significant states - OPEN FIELD: NOTE - this doesn't work from within the script

## Get significant states
sig_sge_of = googlesheets4::read_sheet(
  ss = SHEET_ID,
  sheet = "SGE_OF") %>% 
  dplyr::filter(Variable == "test_fish" & `p-value FDR-adj` < 0.05) %>% 
  dplyr::pull(State)

g <- ggplot_gtable(ggplot_build(polar_all_sge))
strips <- which(startsWith(g$layout$name,'strip'))

strip_vec = ifelse(1:N_STATES %in% sig_sge_of,
                   "#ED474A",
                   "#9EA7A9") 
strip_vec = c(strip_vec[11:N_STATES], strip_vec[6:10], strip_vec[1:5])

for (s in seq_along(strips)) {
  g$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- strip_vec[s]
}

grid::grid.newpage()

png(POLAR_ALL_SGE_SIG_OF,
    width = POL_ALL_WID, 
    height = POL_ALL_HEI, 
    units = "in", 
    res = 400)
grid::grid.draw(g)
dev.off()

# With significant states - NOVEL_OBJECT: NOTE - this doesn't work from within the script

## Get significant states
sig_sge_no = googlesheets4::read_sheet(
  ss = SHEET_ID,
  sheet = "SGE_NO") %>% 
  dplyr::filter(Variable == "test_fish" & `p-value FDR-adj` < 0.05) %>% 
  dplyr::pull(State)

g <- ggplot_gtable(ggplot_build(polar_all_sge))
strips <- which(startsWith(g$layout$name,'strip'))

strip_vec = ifelse(1:N_STATES %in% sig_sge_no,
                   "#ED474A",
                   "#9EA7A9") 
strip_vec = c(strip_vec[11:N_STATES], strip_vec[6:10], strip_vec[1:5])

for (s in seq_along(strips)) {
  g$grobs[[strips[s]]]$grobs[[1]]$children[[1]]$gp$fill <- strip_vec[s]
}

grid::grid.newpage()

png(POLAR_ALL_SGE_SIG_NO,
    width = POL_ALL_WID, 
    height = POL_ALL_HEI, 
    units = "in", 
    res = 400)
grid::grid.draw(g)
dev.off()


