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
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out_split/0.05/dist_angle/15/F0.csv"
LINE_COLS = here::here("config/line_colours/line_colours_0.05.csv")
SHEET_ID = "15hj3N59E4nCFvxH4lES16PPzQ00Iawyf0QWrjaPo3I0"
N_STATES = 15
POLAR_ALL_DGE_SIG_OF = here::here("book/figs/state_freq_F0/0.05/dist_angle/15/polar_dge_sig_of.png")
POLAR_ALL_DGE_SIG_NO = here::here("book/figs/state_freq_F0/0.05/dist_angle/15/polar_dge_sig_no.png")
POLAR_ALL_SGE_SIG_OF = here::here("book/figs/state_freq_F0/0.05/dist_angle/15/polar_sge_sig_of.png")
POLAR_ALL_SGE_SIG_NO = here::here("book/figs/state_freq_F0/0.05/dist_angle/15/polar_sge_sig_no.png")
OUT_DGE = here::here("book/figs/time_dependence/dist_angle/0.08_14_dge.png")
OUT_SGE = here::here("book/figs/time_dependence/dist_angle/0.08_14_sge.png")

# Deauthorise google sheets so that it doesn't ask for prompt
googlesheets4::gs4_deauth()

## True
IN = snakemake@input[["data"]]
LINE_COLS = snakemake@input[["line_cols"]]
OUT_DGE = snakemake@output[["dge"]]
OUT_SGE = snakemake@output[["sge"]]
N_STATES = snakemake@params[["n_states"]] %>% 
  as.numeric()
POLAR_ALL_DGE_SIG_OF = snakemake@input[["polar_all_dge_sig_of"]]
POLAR_ALL_DGE_SIG_NO = snakemake@input[["polar_all_dge_sig_no"]]
POLAR_ALL_SGE_SIG_OF = snakemake@input[["polar_all_sge_sig_of"]]
POLAR_ALL_SGE_SIG_NO = snakemake@input[["polar_all_sge_sig_no"]]

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

# Recode `state`

df = raw %>% 
  dplyr::mutate(state_recode = dplyr::recode(state, !!!recode_vec)) %>% 
  # recode `assay`
  dplyr::mutate(assay = stringr::str_replace(assay, "_", " "),
                assay = factor(assay, levels = c("open field", "novel object"))) %>% 
  # Add `line`
  dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                        fish == "test" ~ test_fish)) %>% 
  # recode outbred
  dplyr::mutate(line = dplyr::if_else(stringr::str_detect(line,
                                                          "outbred"),
                                      "outbred",
                                      line)) %>% 
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
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) 

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
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) 

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
               position = "fill") +
  facet_grid(rows = vars(line),
             cols = vars(assay)) + 
  scale_fill_manual(values = pal_dge_of) +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,200,400,600)) +
  scale_y_continuous(breaks = c(0,0.5,1))

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
  geom_density(aes(seconds, after_stat(count), fill = state_plot_recode),
               position = "fill") +
  facet_grid(rows = vars(line),
             cols = vars(assay)) + 
  scale_fill_manual(values = pal_dge_no) +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,200,400,600)) +
  scale_y_continuous(breaks = c(0,0.5,1))


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
               position = "fill") +
  facet_grid(rows = vars(test_fish),
             cols = vars(assay)) + 
  scale_fill_manual(values = pal_sge_of) +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,200,400,600)) +
  scale_y_continuous(breaks = c(0,0.5,1))

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
               position = "fill") +
  facet_grid(rows = vars(test_fish),
             cols = vars(assay)) + 
  scale_fill_manual(values = pal_sge_no) +
  cowplot::theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,200,400,600)) +
  scale_y_continuous(breaks = c(0,0.5,1))

##########################
# Spatial density: DGE
##########################

# Create density function

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


# OF

sdens_dge_df = df %>% 
  #dplyr::slice_sample(n = 1e5) %>% 
  # filter out iCab references when paired with a different line
  dplyr::filter(fish == "test") %>% 
  # group by line
  dplyr::group_by(assay, line) %>% 
  # get densities
  dplyr::mutate(density = get_density(x, y, n = 30)) %>% 
  dplyr::ungroup()

  # Plot
sdens_dge_of = sdens_dge_df %>% 
  ## take only states 1:4
  dplyr::filter(assay == "open field") %>% 
  ggplot() +
  geom_point(aes(x, y, colour = density),
             alpha = 0.1, size = 0.2) +
  facet_grid(cols = vars(assay), rows = vars(line)) +
  colorspace::scale_color_continuous_sequential(palette = "Mako", rev = F) +
  #scale_colour_viridis_c(option = "rocket") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_x_continuous(breaks = c(0,150,300)) +
  scale_y_continuous(breaks = c(0,150,300)) +
  xlab("x") +
  ylab("y")

# NO

sdens_dge_no = sdens_dge_df %>% 
  dplyr::filter(assay == "novel object") %>% 
  # Plot
  ggplot() +
  geom_point(aes(x, y, colour = density),
             alpha = 0.1, size = 0.2) +
  #coord_polar() +
  facet_grid(cols = vars(assay), rows = vars(line)) +
  colorspace::scale_color_continuous_sequential(palette = "Mako", rev = F) +
  #scale_colour_viridis_c(option = "rocket") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_x_continuous(breaks = c(0,150,300)) +
  scale_y_continuous(breaks = c(0,150,300)) +
  xlab("x") +
  ylab("y")

##########################
# Spatial density: SGE
##########################

# OF

sdens_sge_df = df %>% 
  #dplyr::slice_sample(n = 1e5) %>% 
  # filter out iCab references when paired with a different line
  dplyr::filter(fish == "ref") %>% 
  # group by line
  dplyr::group_by(assay, test_fish) %>% 
  # get densities
  dplyr::mutate(density = get_density(x, y, n = 30)) %>% 
  dplyr::ungroup()

# Plot
sdens_sge_of = sdens_sge_df %>% 
  # take only states 1:4
  dplyr::filter(assay == "open field") %>% 
  ggplot() +
  geom_point(aes(x, y, colour = density),
             alpha = 0.1, size = 0.2) +
  facet_grid(cols = vars(assay), rows = vars(test_fish)) +
  #colorspace::scale_color_continuous_sequential(palette = "Mako", rev = F) +
  scale_colour_viridis_c(option = "rocket") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_x_continuous(breaks = c(0,150,300)) +
  scale_y_continuous(breaks = c(0,150,300)) +
  xlab("x") +
  ylab("y")

# NO

sdens_sge_no = sdens_sge_df %>% 
  dplyr::filter(assay == "novel object") %>% 
  # Plot
  ggplot() +
  geom_point(aes(x, y, colour = density),
             alpha = 0.1, size = 0.2) +
  #coord_polar() +
  facet_grid(cols = vars(assay), rows = vars(test_fish)) +
  #colorspace::scale_color_continuous_sequential(palette = "Mako", rev = F) +
  scale_colour_viridis_c(option = "rocket") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_x_continuous(breaks = c(0,150,300)) +
  scale_y_continuous(breaks = c(0,150,300)) +
  xlab("x") +
  ylab("y")


############################
# Combine and save
############################

final_dge = cowplot::plot_grid(dge_tile_of +
                                 theme(strip.background.y = element_blank(),
                                       strip.text.y = element_blank(),
                                       axis.title.y = element_text(vjust=-5)),
                               time_dens_dge_of,
                               sdens_dge_of,
                               dge_tile_no +
                                 theme(strip.background.y = element_blank(),
                                       strip.text.y = element_blank(),
                                       axis.title.y = element_text(vjust=-5)),
                               time_dens_dge_no,
                               sdens_dge_no,
                               nrow = 2, ncol = 3,
                               rel_widths = c(1,1,0.6,1,1,0.6),
                               align = "hv",
                               labels = c('A', 'B', 'C', 'D', 'E', 'F'))

####################
# With polars coloured by significance
###################
final_dge_raw = cowplot::plot_grid(dge_tile_of +
                                 theme(strip.background.y = element_blank(),
                                       strip.text.y = element_blank(),
                                       axis.title.y = element_text(vjust=-5)),
                               time_dens_dge_of + 
                                 theme(strip.background.y = element_blank(),
                                       strip.text.y = element_blank()),
                               sdens_dge_of,
                               dge_tile_no +
                                 theme(strip.background.y = element_blank(),
                                       strip.text.y = element_blank(),
                                       axis.title.y = element_text(vjust=-5)),
                               time_dens_dge_no +
                                 theme(strip.background.y = element_blank(),
                                       strip.text.y = element_blank()),
                               sdens_dge_no,
                               nrow = 2, ncol = 3,
                               rel_widths = c(1,1,0.6,1,1,0.6),
                               align = "hv",
                               labels = c('B', 'C', 'D', 'E', 'F', 'G'))

final_dge_with_polar = ggdraw() +
  cowplot::draw_image(POLAR_ALL_DGE_SIG_OF,
                      x = 0, y = 0.8,
                      width = 1, height = 0.2) +
  cowplot::draw_plot(final_dge_raw,
                     x = 0, y = 0.2,
                     width = 1, height = 0.6) +
  cowplot::draw_image(POLAR_ALL_DGE_SIG_NO,
                      x = 0, y = 0,
                      width = 1, height = 0.2) +
  cowplot::draw_plot_label(c('A', 'H'),
                           x = c(0,0), y = c(1, 0.2),
                           size = 14)

ggsave(OUT_DGE,
       final_dge_with_polar,
       device = "png",
       width = 11.5,
       height = 20,
       units = "in",
       dpi = 400)

final_sge_raw = cowplot::plot_grid(sge_tile_of +
                                     theme(strip.background.y = element_blank(),
                                           strip.text.y = element_blank(),
                                           axis.title.y = element_text(vjust=-5)),
                                   time_dens_sge_of + 
                                     theme(strip.background.y = element_blank(),
                                           strip.text.y = element_blank()),
                                   sdens_sge_of,
                                   sge_tile_no +
                                     theme(strip.background.y = element_blank(),
                                           strip.text.y = element_blank(),
                                           axis.title.y = element_text(vjust=-5)),
                                   time_dens_sge_no +
                                     theme(strip.background.y = element_blank(),
                                           strip.text.y = element_blank()),
                                   sdens_sge_no,
                                   nrow = 2, ncol = 3,
                                   rel_widths = c(1,1,0.6,1,1,0.6),
                                   align = "hv",
                                   labels = c('B', 'C', 'D', 'E', 'F', 'G'))

final_sge_with_polar = ggdraw() +
  cowplot::draw_image(POLAR_ALL_SGE_SIG_OF,
                      x = 0, y = 0.8,
                      width = 1, height = 0.2) +
  cowplot::draw_plot(final_sge_raw,
                     x = 0, y = 0.2,
                     width = 1, height = 0.6) +
  cowplot::draw_image(POLAR_ALL_SGE_SIG_NO,
                      x = 0, y = 0,
                      width = 1, height = 0.2) +
  cowplot::draw_plot_label(c('A', 'H'),
                           x = c(0,0), y = c(1, 0.2),
                           size = 14)

ggsave(OUT_SGE,
       final_sge_with_polar,
       device = "png",
       width = 11.5,
       height = 20,
       units = "in",
       dpi = 400)


