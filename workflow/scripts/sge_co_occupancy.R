# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out/0.08/dist_angle/15.csv"
LINE_COLS = here::here("config/line_colours/line_colours_0.08.csv")
HEAT_LINES = as.list(c("iCab", "8-2", "18-2", "139-4", "50-2", "38-2", "21-2", "40-1")) %>% 
  unlist
BOX_PNG = here::here("book/figs/sge_F0/co-occupancy/dist_angle/0.08_15_cooc_all.png")
BOX_PDF = here::here("book/figs/sge_F0/co-occupancy/dist_angle/0.08_15_cooc_all.pdf")
HEAT_PNG = here::here("book/figs/sge_F0/co-occupancy/dist_angle/0.08_15_cooc_heat.png")
HEAT_PDF = here::here("book/figs/sge_F0/co-occupancy/dist_angle/0.08_15_cooc_heat.pdf")

## True
IN = snakemake@input[["data"]]
LINE_COLS = snakemake@input[["line_cols"]]
N_STATES = snakemake@params[["n_states"]]
HEAT_LINES = snakemake@params[["selected_lines"]] %>% 
  unlist()
BOX_PNG = snakemake@output[["box_png"]]
BOX_PDF = snakemake@output[["box_pdf"]]
HEAT_PNG = snakemake@output[["heat_png"]]
HEAT_PDF = snakemake@output[["heat_pdf"]]


#######################
# Read in files
#######################

# Read in line colours

line_cols = readr::read_csv(LINE_COLS)
line_vec = line_cols$line
pal = line_cols$colour; names(pal) = line_cols$line


# Read in data

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

ggsave(BOX_PNG,
       cooc_final,
       device = "png",
       width = 13,
       height = 11,
       units = "in",
       dpi = 400)

ggsave(BOX_PDF,
       cooc_final,
       device = "pdf",
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
HEAT_LINES = c(HEAT_LINES, "139-4", "14-2")

# OF
cooc_heat_plot = cooc_heat %>% 
  dplyr::filter(test_fish %in% HEAT_LINES) %>% 
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
  cowplot::theme_cowplot() +
  theme(aspect.ratio = 1) +
  #scale_x_di(breaks = unique(cooc_heat$ref)) +
  #scale_y_di(breaks = unique(cooc_heat$test))  +
  labs(fill = "Frequency\nof state\nco-occupancy\nwithin\nline-pairing") +
  scale_fill_viridis_c(option = "plasma") +
  xlab("reference fish state") +
  ylab("test fish state") +
  theme(strip.text = element_text(face = "bold"),
        axis.text = element_text(size = 8))

ggsave(HEAT_PNG,
       cooc_heat_plot,
       device = "png",
       width = 7,
       height = 20.5,
       units = "in",
       dpi = 400)

ggsave(HEAT_PDF,
       cooc_heat_plot,
       device = "pdf",
       width = 7,
       height = 20.5,
       units = "in",
       dpi = 400)

