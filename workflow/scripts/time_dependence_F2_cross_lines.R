# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)
library(ggridges)
library(viridisLite)
library(googlesheets4)

# Get variables

## Debug
### 0.05 interval
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out_split/0.05/dist_angle/15/F0.csv"
LINE_COLS = here::here("config/line_colours/line_colours_0.05.csv")
SHEET_ID = "118HNbI7ch_0rgRSaDB1b73mpRo5Z2g7uQvNhjZUD7WI"
N_STATES = 15
SELECTED_LINES = list("8-2", "18-2", "21-2", "38-2", "40-1", "50-2") %>% 
  unlist()
POLAR_ALL_DGE_SIG_OF = here::here("book/figs/polar_plots_with_sig/0.05/dist_angle/15/polar_all_dge_of.png")
POLAR_ALL_DGE_SIG_NO = here::here("book/figs/polar_plots_with_sig/0.05/dist_angle/15/polar_all_dge_no.png")
POLAR_ALL_SGE_SIG_OF = here::here("book/figs/polar_plots_with_sig/0.05/dist_angle/15/polar_all_sge_of.png")
POLAR_ALL_SGE_SIG_NO = here::here("book/figs/polar_plots_with_sig/0.05/dist_angle/15/polar_all_sge_no.png")
OUT_DGE_OF = here::here("book/figs/time_dependence_F2_cross_lines/dist_angle/0.05_15_dge_of.png")
OUT_DGE_NO = here::here("book/figs/time_dependence_F2_cross_lines/dist_angle/0.05_15_dge_no.png")
OUT_SGE_OF = here::here("book/figs/time_dependence_F2_cross_lines/dist_angle/0.05_15_sge_of.png")
OUT_SGE_NO = here::here("book/figs/time_dependence_F2_cross_lines/dist_angle/0.05_15_sge_no.png")
### 0.08 interval
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out_split/0.08/dist_angle/15/F0.csv"
LINE_COLS = here::here("config/line_colours/line_colours_0.08.csv")
SHEET_ID = "1kD2ndgmTvGlQw8YZFMztO3YoYfOdSw0SI5_BtomDuaI"
N_STATES = 15
SELECTED_LINES = list("8-2", "18-2", "21-2", "38-2", "40-1", "50-2") %>% 
  unlist()
POLAR_ALL_DGE_SIG_OF = here::here("book/figs/polar_plots_with_sig/0.08/dist_angle/15/polar_all_dge_of.png")
POLAR_ALL_DGE_SIG_NO = here::here("book/figs/polar_plots_with_sig/0.08/dist_angle/15/polar_all_dge_no.png")
POLAR_ALL_SGE_SIG_OF = here::here("book/figs/polar_plots_with_sig/0.08/dist_angle/15/polar_all_sge_of.png")
POLAR_ALL_SGE_SIG_NO = here::here("book/figs/polar_plots_with_sig/0.08/dist_angle/15/polar_all_sge_no.png")
OUT_DGE_OF = here::here("book/figs/time_dependence_F2_cross_lines/dist_angle/select_0.08_15_dge_of.png")
OUT_DGE_NO = here::here("book/figs/time_dependence_F2_cross_lines/dist_angle/select_0.08_15_dge_no.png")
OUT_SGE_OF = here::here("book/figs/time_dependence_F2_cross_lines/dist_angle/select_0.08_15_sge_of.png")
OUT_SGE_NO = here::here("book/figs/time_dependence_F2_cross_lines/dist_angle/select_0.08_15_sge_no.png")


# create output directories
dir.create(dirname(OUT_DGE_OF), recursive = T)
# Deauthorise google sheets so that it doesn't ask for prompt
googlesheets4::gs4_deauth()

## True
IN = snakemake@input[["data"]]
LINE_COLS = snakemake@input[["line_cols"]]
SELECTED_LINES = snakemake@input[["selected_lines"]] %>% 
  unlist()
N_STATES = snakemake@params[["n_states"]] %>% 
  as.numeric()
SHEET_ID = snakemake@params[["sheet_id"]]
OUT_DGE_OF = snakemake@output[["out_dge_of"]]
OUT_DGE_NO = snakemake@output[["out_dge_no"]]
OUT_SGE_OF = snakemake@output[["out_sge_of"]]
OUT_SGE_NO = snakemake@output[["out_sge_no"]]

#######################
# Read in data
#######################

# Get number of rows (for plotting) based on number of states

if (N_STATES == 15){
  N_ROWS = 5
} else if (N_STATES == 12 | 16) {
  N_ROWS = 4
} else if (N_STATES == 17 | 18){
  N_ROWS = 6
}

# Read in line colours

line_cols = readr::read_csv(LINE_COLS)

# Read in file

raw = readr::read_csv(IN) %>% 
  # recode angle to sit between 0 and 360
  dplyr::mutate(angle_recode = ifelse(angle < 0,
                                      180 + (180 + angle),
                                      angle))

# Recode states by mean distance

rank_df = raw %>% 
  dplyr::group_by(state) %>% 
  dplyr::summarise(mean_dist = mean(distance)) %>% 
  # rank
  dplyr::arrange(mean_dist) %>% 
  dplyr::mutate(rank = 1:nrow(.))

recode_vec = rank_df %>% 
  dplyr::pull(rank)
names(recode_vec) = rank_df %>% 
  dplyr::pull(state)

# Further process (including filtering for lines selected in F2 cross)

df = raw %>% 
  # Recode `state`
  dplyr::mutate(state_recode = dplyr::recode(state, !!!recode_vec)) %>% 
  # Filter for selected lines
  dplyr::filter(test_fish %in% SELECTED_LINES) %>% 
  # recode `assay`
  dplyr::mutate(assay = stringr::str_replace(assay, "_", " "),
                assay = factor(assay, levels = c("open field", "novel object"))) %>% 
  # Add `line`
  dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                        fish == "test" ~ test_fish)) %>% 
  # recode outbred
  dplyr::mutate(dplyr::across(c("line", "test_fish"),
                              ~dplyr::if_else(stringr::str_detect(.x,
                                                                  "outbred"),
                                              "outbred",
                                              .x))) %>% 
  # factorise to order `line` and `test_fish`
  dplyr::mutate(line = factor(line, levels = line_cols$line),
                test_fish = factor(test_fish, levels = line_cols$line))



############################
# Get significant states
############################

SIGS_DGE_OF = googlesheets4::read_sheet(SHEET_ID, sheet = "DGE_OF") %>% 
  dplyr::filter(`p-value FDR-adj` < 0.05 & `Variable` == "line") %>% 
  dplyr::pull(State) %>% 
  as.integer()

SIGS_DGE_NO = googlesheets4::read_sheet(SHEET_ID, sheet = "DGE_NO") %>% 
  dplyr::filter(`p-value FDR-adj` < 0.05 & `Variable` == "line") %>% 
  dplyr::pull(State)%>% 
  as.integer()

SIGS_SGE_OF = googlesheets4::read_sheet(SHEET_ID, sheet = "SGE_OF") %>% 
  dplyr::filter(`p-value FDR-adj` < 0.05 & `Variable` == "test_fish") %>% 
  dplyr::pull(State)%>% 
  as.integer()

SIGS_SGE_NO = googlesheets4::read_sheet(SHEET_ID, sheet = "SGE_NO") %>% 
  dplyr::filter(`p-value FDR-adj` < 0.05 & `Variable` == "test_fish") %>% 
  dplyr::pull(State)%>% 
  as.integer()

#######################
# Medarkov matrices: DGE
#######################

SEC_INT = 2

dge_tile_df = df %>% 
  # remove iCab ref fishes (because DGE compares test fishes)
  dplyr::filter(!(line == "iCab" & fish == "ref")) %>% 
  #dplyr::slice_sample(n = 1e6) %>% 
  # add `indiv` column
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv",
               sep = "_") %>%
  # get top state per 2 seconds
  dplyr::mutate(seconds_bin = floor(seconds / SEC_INT)) %>% 
  dplyr::group_by(assay, indiv, line, seconds_bin) %>% 
  dplyr::count(state_recode) %>% 
  dplyr::slice_max(order_by = n, n = 1) %>% 
  dplyr::ungroup() %>% 
  # reverse order by `indiv` so that the earliest videos are at the top
  dplyr::arrange(indiv) %>% 
  # convert `seconds_bin` back to seconds
  dplyr::mutate(seconds = seconds_bin * SEC_INT)

# Open field

dge_tile_of = dge_tile_df %>% 
  dplyr::filter(assay == "open field") %>% 
  ggplot() +
  geom_tile(aes(seconds, indiv, fill = state_recode)) + 
  facet_grid(rows = vars(line), cols = vars(assay), scales = "free") +
  scale_fill_viridis_c() +
  scale_y_discrete(limits = rev) +
  guides(fill = "none") +
  ylab("individual fish") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) 

# Novel object

dge_tile_no = dge_tile_df %>% 
  dplyr::filter(assay == "novel object") %>% 
  ggplot() +
  geom_tile(aes(seconds, indiv, fill = state_recode)) + 
  facet_grid(rows = vars(line), cols = vars(assay), scales = "free") +
  scale_fill_viridis_c() +
  scale_y_discrete(limits = rev) +
  guides(fill = "none") +
  ylab("individual fish") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) 

# Combine 

dge_tiles_final = cowplot::plot_grid(dge_tile_of,
                                     dge_tile_no,
                                     align = "hv",
                                     labels = c("A", "B"))

# Save 

#ggsave(TILE_DGE,
#       dge_tiles_final,
#       device = "png",
#       width = 13,
#       height = 22,
#       units = "in",
#       dpi = 400)

#######################
# Medarkov matrices: SGE
#######################

sge_tile_df = df %>% 
  # take only iCab test fishes 
  dplyr::filter(fish == "ref") %>% 
  # add `indiv` column
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv",
               sep = "_") %>%
  # rename and reorder assay
  dplyr::mutate(assay = stringr::str_replace(assay, "_", " "),
                assay = factor(assay, levels = c("open field", "novel object"))) %>% 
  # get top state per 2 seconds
  dplyr::mutate(seconds_bin = floor(seconds / SEC_INT)) %>% 
  dplyr::group_by(assay, indiv, test_fish, seconds_bin) %>% 
  dplyr::count(state_recode) %>% 
  dplyr::slice_max(order_by = n, n = 1) %>% 
  dplyr::ungroup() %>% 
  # reverse order by `indiv` so that the earliest videos are at the top
  dplyr::arrange(indiv) %>% 
  # convert `seconds_bin` back to seconds
  dplyr::mutate(seconds = seconds_bin * SEC_INT)

# Open field

sge_tile_of = sge_tile_df %>% 
  dplyr::filter(assay == "open field") %>% 
  ggplot() +
  geom_tile(aes(seconds, indiv, fill = state_recode)) + 
  facet_grid(rows = vars(test_fish), cols = vars(assay), scales = "free") +
  scale_fill_viridis_c(option = "inferno") +
  scale_y_discrete(limits = rev) +
  guides(fill = "none") +
  ylab("individual fish") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) 

# Novel object

sge_tile_no = sge_tile_df %>% 
  dplyr::filter(assay == "novel object") %>% 
  ggplot() +
  geom_tile(aes(seconds, indiv, fill = state_recode)) + 
  facet_grid(rows = vars(test_fish), cols = vars(assay), scales = "free") +
  scale_fill_viridis_c(option = "inferno") +
  scale_y_discrete(limits = rev) +
  guides(fill = "none") +
  ylab("individual fish") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) 

# Combine 

sge_tiles_final = cowplot::plot_grid(sge_tile_of,
                                     sge_tile_no,
                                     align = "hv",
                                     labels = c("A", "B"))

# Save 

#ggsave(TILE_SGE,
#       sge_tiles_final,
#       device = "png",
#       width = 13,
#       height = 22,
#       units = "in",
#       dpi = 400)

##########################
# Time density - DGE
##########################

# Take viridis colours for significant states and add grey
pal_dge_of = viridisLite::viridis(n = N_STATES)
pal_dge_of = c(pal_dge_of[SIGS_DGE_OF], "#9da2ab")
names(pal_dge_of) = c(as.character(SIGS_DGE_OF), "other")

pal_dge_no = viridisLite::viridis(n = N_STATES)
pal_dge_no = c(pal_dge_no[SIGS_DGE_NO], "#9da2ab")
names(pal_dge_no) = c(as.character(SIGS_DGE_NO), "other")

ASSAY = "open field"
time_dens_dge_of = df %>% 
  dplyr::filter(assay == ASSAY) %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(fish == "test") %>% 
  # filter for target assay
  #dplyr::filter(assay == ASSAY) %>% 
  # recode state 
  dplyr::mutate(state_plot_recode = dplyr::case_when(state_recode %in% SIGS_DGE_OF ~ as.character(state_recode),
                                                     TRUE ~ "other"),
                state_plot_recode = factor(state_plot_recode, levels = c(as.character(SIGS_DGE_OF), "other"))) %>% 
  ggplot() +
  geom_density(aes(seconds, after_stat(count), fill = state_plot_recode),
               size = 0.25,
               position = "fill") +
  facet_grid(rows = vars(line),
             cols = vars(assay)) + 
  scale_fill_manual(values = pal_dge_of) +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        ) +
  guides(fill = "none") +
  ylab("frequency") +
  scale_x_continuous(breaks = c(0,200,400,600)) +
  scale_y_continuous(breaks = c(0,1))

ASSAY = "novel object"
time_dens_dge_no = df %>% 
  dplyr::filter(assay == ASSAY) %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(fish == "test") %>% 
  # filter for target assay
  #dplyr::filter(assay == ASSAY) %>% 
  # recode state 
  dplyr::mutate(state_plot_recode = dplyr::case_when(state_recode %in% SIGS_DGE_NO ~ as.character(state_recode),
                                                     TRUE ~ "other"),
                state_plot_recode = factor(state_plot_recode, levels = c(as.character(SIGS_DGE_NO), "other"))) %>% 
  ggplot() +
  geom_density(aes(seconds, after_stat(count),fill = state_plot_recode),
               size = 0.25,
               position = "fill") +
  facet_grid(rows = vars(line),
             cols = vars(assay)) + 
  scale_fill_manual(values = pal_dge_no) +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0)) +
  guides(fill = "none") +
  ylab("frequency") +
  scale_x_continuous(breaks = c(0,200,400,600)) +
  scale_y_continuous(breaks = c(0,1))

# Combine 

dge_dens_final = cowplot::plot_grid(time_dens_dge_of,
                                    time_dens_dge_no,
                                    align = "hv",
                                    labels = c("A", "B"))



##########################
# Time density: SGE
##########################


# Take viridis colours for significant states and add grey
pal_sge_of = viridisLite::viridis(n = N_STATES, option = "inferno")
pal_sge_of = c(pal_sge_of[SIGS_SGE_OF], "#9da2ab")
names(pal_sge_of) = c(as.character(SIGS_SGE_OF), "other")

pal_sge_no = viridisLite::viridis(n = N_STATES, option = "inferno")
pal_sge_no = c(pal_sge_no[SIGS_SGE_NO], "#9da2ab")
names(pal_sge_no) = c(as.character(SIGS_SGE_NO), "other")

ASSAY = "open field"
time_dens_sge_of = df %>% 
  dplyr::filter(assay == ASSAY) %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(fish == "ref") %>% 
  # filter for target assay
  #dplyr::filter(assay == ASSAY) %>% 
  # recode state 
  dplyr::mutate(state_plot_recode = dplyr::case_when(state_recode %in% SIGS_SGE_OF ~ as.character(state_recode),
                                                     TRUE ~ "other"),
                state_plot_recode = factor(state_plot_recode, levels = c(as.character(SIGS_SGE_OF), "other"))) %>% 
  ggplot() +
  geom_density(aes(seconds, after_stat(count), fill = state_plot_recode),
               size = 0.25,
               position = "fill") +
  facet_grid(rows = vars(test_fish),
             cols = vars(assay)) + 
  scale_fill_manual(values = pal_sge_of) +
  cowplot::theme_cowplot(font_size = 12,
                         rel_small = 10/14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        ) +
  ylab("frequency") +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,200,400,600)) +
  scale_y_continuous(breaks = c(0,1))

ASSAY = "novel object"
time_dens_sge_no = df %>% 
  dplyr::filter(assay == ASSAY) %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(fish == "ref") %>% 
  # filter for target assay
  #dplyr::filter(assay == ASSAY) %>% 
  # recode state 
  dplyr::mutate(state_plot_recode = dplyr::case_when(state_recode %in% SIGS_SGE_NO ~ as.character(state_recode),
                                                     TRUE ~ "other"),
                state_plot_recode = factor(state_plot_recode, levels = c(as.character(SIGS_SGE_NO), "other"))) %>% 
  ggplot() +
  geom_density(aes(seconds, after_stat(count), fill = state_plot_recode),
               size = 0.25,
               position = "fill") +
  facet_grid(rows = vars(test_fish),
             cols = vars(assay)) + 
  scale_fill_manual(values = pal_sge_no) +
  cowplot::theme_cowplot(font_size = 12,
                         rel_small = 10/14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        ) +
  ylab("frequency") +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,200,400,600)) +
  scale_y_continuous(breaks = c(0,1))

# Combine 

sge_dens_final = cowplot::plot_grid(time_dens_sge_of,
                                    time_dens_sge_no,
                                    align = "hv",
                                    labels = c("A", "B"))

# Save 

#ggsave(DENS_SGE,
#       sge_dens_final,
#       device = "png",
#       width = 13,
#       height = 22,
#       units = "in",
#       dpi = 400)


##########################
# Compose final figures
##########################

##########
# DGE
##########

final_dge = cowplot::plot_grid(dge_tile_of + 
                                 theme(strip.text.y = element_blank()),
                               time_dens_dge_of,
                               dge_tile_no + 
                                 theme(strip.text.y = element_blank()),
                               time_dens_dge_no,
                               nrow = 2, ncol = 2,
                               align = "hv",labels = c("B", "C", "D", "E"))

## Open field

final_dge_of = cowplot::plot_grid(dge_tile_of + 
                                    theme(strip.text.y = element_blank()),
                                  time_dens_dge_of,
                                  ncol = 2,
                                  align = "hv", axis = "tblr",
                                  labels = c("B", "C"),
                                  label_size = 14)

final_dge_of = cowplot::ggdraw() +
  cowplot::draw_image(POLAR_ALL_DGE_SIG_OF,
                      x = 0, y = 0.6,
                      width = 1, height = 0.4) +
  cowplot::draw_plot(final_dge_of,
                     x = 0, y = 0,
                     width = 1, height = 0.6) +
  cowplot::draw_plot_label('A',
                           size = 14)


ggsave(OUT_DGE_OF,
       final_dge_of,
       device = "png",
       width = 10,
       height = 15,
       units = "in",
       dpi = 400)

## Novel object

final_dge_no = cowplot::plot_grid(dge_tile_no + 
                                    theme(strip.text.y = element_blank()),
                                  time_dens_dge_no,
                                  ncol = 2,
                                  align = "hv", axis = "tblr",
                                  labels = c("B", "C"),
                                  label_size = 14)

final_dge_no = cowplot::ggdraw() +
  cowplot::draw_image(POLAR_ALL_DGE_SIG_NO,
                      x = 0, y = 0.6,
                      width = 1, height = 0.4) +
  cowplot::draw_plot(final_dge_no,
                     x = 0, y = 0,
                     width = 1, height = 0.6) +
  cowplot::draw_plot_label('A',
                           size = 14)


ggsave(OUT_DGE_NO,
       final_dge_no,
       device = "png",
       width = 10,
       height = 15,
       units = "in",
       dpi = 400)

##########
# SGE
##########

# Open field 

final_sge_of = cowplot::plot_grid(sge_tile_of + 
                                    theme(strip.text.y = element_blank()),
                                  time_dens_sge_of,
                                  ncol = 2,
                                  align = "hv", axis = "tblr",
                                  labels = c("B", "C"),
                                  label_size = 14)

final_sge_of = cowplot::ggdraw() +
  cowplot::draw_image(POLAR_ALL_SGE_SIG_OF,
                      x = 0, y = 0.6,
                      width = 1, height = 0.4) +
  cowplot::draw_plot(final_sge_of,
                     x = 0, y = 0,
                     width = 1, height = 0.6) +
  cowplot::draw_plot_label('A',
                           size = 14)


ggsave(OUT_SGE_OF,
       final_sge_of,
       device = "png",
       width = 10,
       height = 15,
       units = "in",
       dpi = 400)

# Novel object

final_sge_no = cowplot::plot_grid(sge_tile_no + 
                                    theme(strip.text.y = element_blank()),
                                  time_dens_sge_no,
                                  ncol = 2,
                                  align = "hv", axis = "tblr",
                                  labels = c("B", "C"),
                                  label_size = 14)

final_sge_no = cowplot::ggdraw() +
  cowplot::draw_image(POLAR_ALL_SGE_SIG_NO,
                      x = 0, y = 0.6,
                      width = 1, height = 0.4) +
  cowplot::draw_plot(final_sge_no,
                     x = 0, y = 0,
                     width = 1, height = 0.6) +
  cowplot::draw_plot_label('A',
                           size = 14)


ggsave(OUT_SGE_NO,
       final_sge_no,
       device = "png",
       width = 10,
       height = 15,
       units = "in",
       dpi = 400)
