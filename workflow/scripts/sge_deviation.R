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
HEAT_LINES = as.list(c("iCab", "8-2", "18-2", "50-2", "38-2", "21-2", "40-1")) %>% 
  unlist
N_STATES = 15
PNG = here::here("book/figs/sge_F0/co-occupancy/dist_angle/0.08_15_cooc_heat.png")
PDF = here::here("book/figs/sge_F0/co-occupancy/dist_angle/0.08_15_cooc_heat.pdf")

## True
IN = snakemake@input[["data"]]
LINE_COLS = snakemake@input[["line_cols"]]
N_STATES = snakemake@params[["n_states"]]
PNG = snakemake@output[["png"]]
PDF = snakemake@output[["pdf"]]


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
  # take only iCabs
  dplyr::filter(line == "iCab") %>% 
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
# iCab control boxplot
#######################

control_df = sge_df %>% 
  # take only iCabs paired with another iCab
  dplyr::filter(test_fish == "iCab") %>% 
  # count the number of observations per state for each individual/assay combination
  dplyr::count(assay, indiv, state_recode) %>% 
  # add a column with total count for each assay/individual
  dplyr::add_count(assay, indiv, wt = n, name = "nn") %>% 
  dplyr::ungroup() %>% 
  # get a frequency for each state
  dplyr::mutate(FREQ = n / nn) 


control_box = control_df %>% 
  dplyr::mutate(state_recode = factor(state_recode, levels = 1:N_STATES)) %>% 
  ggplot() + 
    geom_boxplot(aes(state_recode, FREQ, fill = state_recode)) +
    facet_grid(rows = vars(assay)) +
    cowplot::theme_cowplot() +
    xlab("HMM state") +
    ylab("frequency of time spent in HMM state\n(iCab-iCab pairings)") +
    scale_fill_viridis_d(option = "plasma") +
    guides(fill = "none")

#######################
# Deviation from control median
#######################

# Calculate median frequency for iCab controls

#control_median_df = control_df %>% 
#  dplyr::group_by(assay, state_recode) %>% 
#  dplyr::summarise(CONTROL_MED_FREQ = median(FREQ, na.rm = T)) %>% 
#  dplyr::ungroup()

# Calculate median frequency when with other lines

freq_df = sge_df %>% 
  # count the number of observations per state for each individual/assay combination
  dplyr::count(assay, indiv, test_fish, state_recode) %>% 
  # add a column with total count for each assay/individual
  dplyr::add_count(assay, indiv, wt = n, name = "nn") %>% 
  dplyr::ungroup() %>% 
  # get a frequency for each state
  dplyr::mutate(FREQ = n / nn) 


# Split into list by test fish and assay
freq_list = freq_df %>% 
  split(., f = ~ test_fish + assay)

# For each state, run a Mann-Whitney test between
LINE_ASSAYS = names(freq_list)[-which(grepl("iCab", names(freq_list)))]; names(LINE_ASSAYS) = LINE_ASSAYS

mw_out = purrr::map_dfr(LINE_ASSAYS, function(LINE_ASSAY){
  # split LINE_ASSAY
  TEST_LINE = stringr::str_split(LINE_ASSAY, "\\.", simplify = T) %>%
    .[[1]]
  ASSAY = stringr::str_split(LINE_ASSAY, "\\.", simplify = T) %>%
    .[[2]]
  # Get name of iCab list element
  ICAB_NAME = paste("iCab", ASSAY, sep = ".")
  # Run loop over each state for Mann-U test
  purrr::map_dfr(1:N_STATES, function(STATE){
    # Run Mann-U test
    man_out = wilcox.test(freq_list[[ICAB_NAME]] %>% 
                  dplyr::filter(state_recode == STATE) %>% 
                  dplyr::pull(FREQ),
                freq_list[[LINE_ASSAY]] %>% 
                  dplyr::filter(state_recode == STATE) %>% 
                  dplyr::pull(FREQ))
    # Welch's t-test
    ttest_out = t.test(freq_list[[ICAB_NAME]] %>% 
                         dplyr::filter(state_recode == STATE) %>% 
                         dplyr::pull(FREQ),
                       freq_list[[LINE_ASSAY]] %>% 
                         dplyr::filter(state_recode == STATE) %>% 
                         dplyr::pull(FREQ))
    # Pull out statistics into tibble
    out = tibble::tibble(ASSAY = ASSAY,
                         TEST_LINE = TEST_LINE,
                         MANN_U = man_out$statistic,
                         MANN_P = man_out$p.value,
                         T_STAT = ttest_out$statistic,
                         T_P = ttest_out$p.value)
  }, .id = "STATE")
})

# Sum statistics within combinations of assay and test line
mw_final = mw_out %>% 
  dplyr::group_by(ASSAY, TEST_LINE) %>% 
  dplyr::summarise(SUM_MW = sum(MANN_U, na.rm = T),
                   SUM_P = sum(MANN_P, na.rm = T),
                   SUM_T = sum(abs(T_STAT))) %>% 
  dplyr::ungroup() %>% 
  # order TEST_LINE
  dplyr::mutate(TEST_LINE = factor(TEST_LINE, levels = line_vec))

# Plot

## Remove non-MIKK lines from `line_vec`
MIKK_VEC = line_vec[-which(line_vec %in% c("iCab", "Kiyosu CC", "luzonensis"))]

## Summed T-statistic
t_fig = mw_final %>% 
  # remove Kiyosu CC and luzonensis
  dplyr::filter(!TEST_LINE %in% c("Kiyosu CC", "luzonensis")) %>% 
  # order assay
  dplyr::mutate(ASSAY = factor(ASSAY, levels = c("open field", "novel object"))) %>% 
  dplyr::mutate(TEST_LINE = factor(TEST_LINE, levels = rev(MIKK_VEC))) %>% 
  ggplot() + 
  geom_col(aes(TEST_LINE, SUM_T, fill = TEST_LINE)) +
  scale_fill_manual(values = pal) +
  facet_grid(cols = vars(ASSAY))   +
  #scale_x_discrete(limits = rev(levels(MIKK_VEC))) +
  coord_flip() +
  cowplot::theme_cowplot() +
  guides(fill = "none")   +
  ylab("absolute t-statistic comparing state frequency for\niCab-iCab vs iCab-MIKK pairings,\nsummed across all states") +
  xlab("test fish line")

#######################
# Compile and save fig
#######################

final = cowplot::plot_grid(control_box,
                           t_fig,align = "hv",axis = "tblr", rel_widths = c(0.45,0.55), labels = c("A", "B"),label_size = 16)


ggsave(PNG,
       final,
       device = "png",
       width = 12,
       height = 12,
       units = "in",
       dpi = 400)

ggsave(PDF,
       final,
       device = "pdf",
       width = 12,
       height = 12,
       units = "in",
       dpi = 400)

