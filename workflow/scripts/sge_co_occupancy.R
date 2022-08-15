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
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out/0.08/dist_angle/15.csv"
N_STATES = 15
INTERVAL = 0.08
LINE_COLS = here::here("config/line_colours/line_colours_0.08.csv")

OUT_COOC_ALL = here::here("book/figs/sge_F0/co-occupancy/dist_angle/0.08_15_cooc_all.png")

## True
IN = snakemake@input[["data"]]
LINE_COLS = snakemake@input[["line_cols"]]
N_STATES = snakemake@params[["n_states"]]
OUT_HEAT = snakemake@output[["heatmap"]]
OUT_BOX_ALL = snakemake@output[["boxplot_all"]]
OUT_PER_STATE = snakemake@output[["box_and_heat_per_state"]]


# Read in line colours
line_cols = readr::read_csv(LINE_COLS)
line_vec = line_cols$line
pal = line_cols$colour; names(pal) = line_cols$line


# Get number of rows (for plotting) based on number of states

if (N_STATES == 15){
  N_ROWS = 5
} else if (N_STATES == 12 | N_STATES == 16) {
  N_ROWS = 4
} else if (N_STATES == 17 | N_STATES == 18){
  N_ROWS = 6
} else if (N_STATES == 14){
  N_ROWS = 2
}

# Read in file

df = readr::read_csv(IN) %>% 
  # recode angle to sit between 0 and 360
  dplyr::mutate(angle_recode = ifelse(angle < 0,
                                      180 + (180 + angle),
                                      angle)) %>% 
  # convert `time` to character and add a 0 if only 3 characters
  dplyr::mutate(time = as.character(time),
                time = dplyr::if_else(nchar(time) == 3,
                                      paste("0", time, sep = ""),
                                      time))

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

# Get SGE df

sge_df = df %>% 
  # filter for F0
  dplyr::filter(dataset == "F0") %>% 
  # rename outbred
  dplyr::mutate(test_fish = dplyr::if_else(stringr::str_detect(test_fish, "outbred"),
                                           "Kiyosu CC",
                                           test_fish)) %>% 
  ## add `line`
  dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                        fish == "test" ~ test_fish)) %>% 
  # order `test_fish`
  dplyr::mutate(test_fish = factor(test_fish, levels = line_vec)) %>% 
  # rename and reorder assay
  dplyr::mutate(assay = stringr::str_replace(assay, "_", " "),
                assay = factor(assay, levels = c("open field", "novel object"))) %>% 
  # add `run` column
  tidyr::unite(date, time, quadrant,
               col = "run",
               sep = "_",
               remove = F) %>% 
  # add `indiv` column
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv",
               sep = "_",
               remove = F)

#######################
# Co-occupancy boxplot -- all states
#######################

cooc = sge_df %>% 
  # pivot wider to get cols for ref and test
  tidyr::pivot_wider(id_cols = c("run", "assay", "test_fish", "seconds"),
                     names_from = fish,
                     values_from = state_recode) %>% 
  # group by run and assay
  dplyr::group_by(run, assay, test_fish) %>% 
  summarise(TOTAL_ROWS = n(),
            TOTAL_CONC = sum(ref == test, na.rm = T),
            FREQ_CONC = TOTAL_CONC / TOTAL_ROWS) %>% 
  dplyr::ungroup()

# Get KW stat
kw_all = cooc %>% 
  dplyr::group_by(assay) %>% 
  rstatix::kruskal_test(FREQ_CONC ~ test_fish) %>% 
  rstatix::adjust_pvalue(method = "fdr") %>% 
  rstatix::add_significance(p.col = "p.adj") %>% 
  #dplyr::mutate(p.adj = as.character(signif(p.adj, digits = 3))) %>% 
  # paste p-value and significance together
  dplyr::mutate(p_final = dplyr::case_when(p.adj.signif == "ns" ~ paste("p =", scales::scientific(p.adj, digits = 2)),
                                           TRUE ~ paste("p =", scales::scientific(p.adj, digits = 2), p.adj.signif)))

# Plot

box_all = cooc %>% 
  #dplyr::mutate(test_fish = factor(test_fish, levels = line_vec)) %>% 
  ggplot() +
  geom_boxplot(aes(test_fish, FREQ_CONC, fill = test_fish)) +
  scale_fill_manual(values = pal) +
  facet_grid(cols = vars(assay)) +
  geom_text(data = kw_all,
            aes(x = "33-1", y = 0.45, label = p_final),
            size = 4) +
  cowplot::theme_cowplot() +
  guides(fill = "none") +
  xlab("test fish") +
  ylab("frequency of state co-occupancy") +
  labs(fill = "test fish") +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(cooc$test_fish)))
  
  



#######################
# Co-occupancy median -- all states
#######################

col_all = cooc %>% 
  dplyr::group_by(assay, test_fish) %>% 
  summarise(median_cooc = median(FREQ_CONC)) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
    geom_col(aes(test_fish, median_cooc, fill = test_fish)) +
    scale_fill_manual(values = pal) +
    facet_grid(cols = vars(assay)) +  
    xlab("test fish") +
    ylab("line median for\nfrequency of state co-occupancy") +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(cooc$test_fish))) +
    cowplot::theme_cowplot() +
    guides(fill = "none") 


cooc_final = cowplot::plot_grid(box_all,
                                col_all +
                                  theme(axis.title.y = element_blank()),
                                align = "hv",
                                ncol = 2,
                                labels = c("A", "B"),
                                label_size = 16)

ggsave(OUT_COOC_ALL,
       cooc_final,
       device = "png",
       width = 13,
       height = 11,
       units = "in",
       dpi = 400)

#######################
# Co-occupancy heatmap
#######################

cooc_heat = sge_df %>% 
  # pivot wider to get cols for ref and test
  tidyr::pivot_wider(id_cols = c("run", "assay", "test_fish", "seconds"),
                     names_from = fish,
                     values_from = state_recode) %>% 
  dplyr::group_by(assay, test_fish) %>% 
  dplyr::count(assay, ref, test) %>% 
  dplyr::add_count(assay, test_fish, wt = n, name = "nn") %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(FREQ_COOC = n / nn) 

# Choose lines with high/low co-occupancy for heatmaps
heat_lines = c("iCab", "8-2", "18-2", "139-4", "50-2", "38-2", "21-2", "40-1")

# OF
cooc_heat_plot = cooc_heat %>% 
  dplyr::filter(test_fish %in% heat_lines) %>% 
  # recode NAs as character
  dplyr::mutate(across(c(ref, test),
                       ~as.character(.))) %>% 
  # replace NA with character
  dplyr::mutate(across(c(ref, test),
                       ~tidyr::replace_na(., "NA"))) %>% 
  # convert to factor for order
  dplyr::mutate(across(c(ref, test),
                       ~factor(., levels = c(seq(1:N_STATES), "NA")))) %>% 
  ggplot() +
  geom_tile(aes(ref, test, fill = FREQ_COOC)) +
  facet_grid(rows = vars(test_fish),
             cols = vars(assay)) +
  cowplot::theme_cowplot(font_size = FONT_SIZE) +
  theme(aspect.ratio = 1) +
  #scale_x_di(breaks = unique(cooc_heat$ref)) +
  #scale_y_di(breaks = unique(cooc_heat$test))  +
  labs(fill = "Frequency\nof state\nco-occupancy\nwithin\nline-pairing") +
  scale_fill_viridis_c(option = "plasma") +
  xlab("reference fish state") +
  ylab("test fish state")

ggsave(OUT_HEAT,
       cooc_heat_plot,
       device = "png",
       width = 18,
       height = 7,
       units = "in",
       dpi = 400)

#######################
# Final figure
#######################

# Put together
final = cowplot::ggdraw() +
  cowplot::draw_plot(polar,
                     x = 0, y= 0.7,
                     width = 1, height = 0.3) +
  cowplot::draw_plot(box_per_state_of,
                     x = 0, y = 0.3,
                     width = 0.5, height = 0.4) +
  cowplot::draw_plot(box_per_state_no +
                       theme(axis.text.y=element_blank(),
                             axis.ticks.y=element_blank(),
                             axis.line.y = element_blank()) +
                       ylab(NULL),
                     x = 0.5, y = 0.3,
                     width = 0.5, height = 0.4) +
  cowplot::draw_plot(cooc_heat_plot,
                     x = 0, y = 0,
                     width = 1, height = 0.3) +
  cowplot::draw_plot_label(c("A", "B", "C"),
                           x = c(0,0,0),
                           y = c(1, 0.7, 0.3))

ggsave(OUT_PER_STATE,
       final,
       device = "png",
       width = 18,
       height = 20,
       units = "in",
       dpi = 400)


