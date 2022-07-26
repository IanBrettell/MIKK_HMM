# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out/0.05/dist_angle/15.csv"
LINE_COLS = here::here("config/line_colours/line_colours_0.08.csv")
N_STATES = 15
SHEET_ID = "118HNbI7ch_0rgRSaDB1b73mpRo5Z2g7uQvNhjZUD7WI"
DGE_HIST = here::here("book/figs/state_freq_F0_select/dist_angle/0.05_15_state_freq_F0_select_dge.png")
SGE_HIST = here::here("book/figs/state_freq_F0_select/dist_angle/0.05_15_state_freq_F0_select_sge.png")
SELECTED_LINES = list("8-2", "18-2", "21-2", "38-2", "40-1", "50-2") %>% 
  unlist()

## True
IN = snakemake@input[["data"]]
LINE_COLS = snakemake@input[["line_cols"]]
SELECTED_LINES = snakemake@input[["selected_lines"]] %>% 
  unlist()
N_STATES = snakemake@params[["n_states"]]
SHEET_ID = snakemake@params[["sheet_id"]]
DGE_HIST = snakemake@output[["dge_hist"]]
SGE_HIST = snakemake@output[["dge_hist"]]


# Deauthorise google sheets so that it doesn't ask for prompt
#googlesheets4::gs4_deauth()

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


# FILTER FOR F0 and selected lines
df = df %>% 
  dplyr::filter(dataset == "F0") %>% 
  dplyr::filter(test_fish %in% SELECTED_LINES)

## DGE

# Get proportion of time each fish spent in each state
df_dge = df %>% 
  # remove iCab when paired with a different test fish
  dplyr::filter(!(fish == "ref" & test_fish != "iCab")) %>% 
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

### Calculate variance explained

aov_dge = df_dge %>% 
  dplyr::group_by(assay, state_recode) %>% 
  tidyr::nest() %>%
  dplyr::mutate(model = purrr::map(data, ~aov(
    state_freq_invnorm ~ date + time + quadrant + tank_side + line,
    data = .))) %>%
  select(-data) %>% 
  dplyr::mutate(model_tidy = purrr::map(model, broom::tidy)) %>%
  tidyr::unnest(model_tidy) %>% 
  rstatix::adjust_pvalue(p.col = "p.value", method = "fdr") %>% 
  rstatix::add_significance(p.col = "p.value.adj") %>% 
  # reduce to 3 digits
  dplyr::mutate(dplyr::across(c("sumsq", "meansq", "statistic", "p.value", "p.value.adj"),
                              ~signif(.x, digits = 3)))


### Save to Google Drive

# Open field
dge_tbl_of = aov_dge %>% 
  dplyr::filter(assay == "open field") %>% 
  # add variance explained
  dplyr::group_by(assay, state_recode) %>% 
  dplyr::mutate(tot_ss = sum(sumsq)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(var_expl_perc = (sumsq / tot_ss) * 100 ) %>%
  # select and rename key columns
  dplyr::select(State = state_recode,
                Variable = term,
                Statistic = statistic,
                `p-value` = p.value,
                `p-value FDR-adj` = p.value.adj,
                `Significance (p-value FDR-adj)` = p.value.adj.signif,
                `Variance explained (%)` = var_expl_perc) %>% 
  # show only 2 decimals
  dplyr::mutate(dplyr::across(c("Statistic", 
                                #"p-value", 
                                "Variance explained (%)"),
                              ~ format(round(.x, 2), nsmall = 2)))

## Write to google sheet
googlesheets4::write_sheet(
  data = dge_tbl_of,
  ss = SHEET_ID,
  sheet = "DGE_OF")

# Novel object
dge_tbl_no = aov_dge %>% 
  dplyr::filter(assay == "novel object") %>% 
  # add variance explained
  dplyr::group_by(assay, state_recode) %>% 
  dplyr::mutate(tot_ss = sum(sumsq)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(var_expl_perc = (sumsq / tot_ss) * 100 ) %>%
  # select and rename key columns
  dplyr::select(State = state_recode,
                Variable = term,
                Statistic = statistic,
                `p-value` = p.value,
                `p-value FDR-adj` = p.value.adj,
                `Significance (p-value FDR-adj)` = p.value.adj.signif,
                `Variance explained (%)` = var_expl_perc) %>% 
  # show only 2 decimals
  dplyr::mutate(dplyr::across(c("Statistic", 
                                #"p-value", 
                                "Variance explained (%)"),
                              ~ format(round(.x, 2), nsmall = 2)))

## Write to google sheet
googlesheets4::write_sheet(
  data = dge_tbl_no,
  ss = SHEET_ID,
  sheet = "DGE_NO")


## SGE


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

### Calculate variance explained

aov_sge = df_sge %>% 
  dplyr::group_by(assay, state_recode) %>% 
  tidyr::nest() %>%
  dplyr::mutate(model = purrr::map(data, ~aov(
    state_freq_invnorm ~ date + time + quadrant + tank_side + test_fish,
    data = .))) %>%
  select(-data) %>% 
  dplyr::mutate(model_tidy = purrr::map(model, broom::tidy)) %>%
  tidyr::unnest(model_tidy) %>% 
  rstatix::adjust_pvalue(p.col = "p.value", method = "fdr") %>% 
  rstatix::add_significance(p.col = "p.value.adj") %>% 
  # reduce to 3 digits
  dplyr::mutate(dplyr::across(c("sumsq", "meansq", "statistic", "p.value", "p.value.adj"),
                              ~signif(.x, digits = 3)))

DT::datatable(aov_sge %>% 
                dplyr::select(-model),
              options = list(pageLength = 20))


### Save to Google Drive


# Open field
sge_tbl_of = aov_sge %>% 
  dplyr::filter(assay == "open field") %>% 
  dplyr::select(-model) %>% 
  # add variance explained
  dplyr::group_by(assay, state_recode) %>% 
  dplyr::mutate(tot_ss = sum(sumsq)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(var_expl_perc = (sumsq / tot_ss) * 100 ) %>%
  # select and rename key columns
  dplyr::select(State = state_recode,
                Variable = term,
                Statistic = statistic,
                `p-value` = p.value,
                `p-value FDR-adj` = p.value.adj,
                `Significance (p-value FDR-adj)` = p.value.adj.signif,
                `Variance explained (%)` = var_expl_perc) %>% 
  # show only 2 decimals
  dplyr::mutate(dplyr::across(c("Statistic", 
                                #"p-value", 
                                "Variance explained (%)"),
                              ~ format(round(.x, 2), nsmall = 2)))

## Write to Google sheet
googlesheets4::write_sheet(
  data = sge_tbl_of,
  ss = SHEET_ID,
  sheet = "SGE_OF")

# Novel object
sge_tbl_no = aov_sge %>% 
  dplyr::filter(assay == "novel object") %>% 
  dplyr::select(-model) %>% 
  # add variance explained
  dplyr::group_by(assay, state_recode) %>% 
  dplyr::mutate(tot_ss = sum(sumsq)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(var_expl_perc = (sumsq / tot_ss) * 100 ) %>%
  # select and rename key columns
  dplyr::select(State = state_recode,
                Variable = term,
                Statistic = statistic,
                `p-value` = p.value,
                `p-value FDR-adj` = p.value.adj,
                `Significance (p-value FDR-adj)` = p.value.adj.signif,
                `Variance explained (%)` = var_expl_perc) %>% 
  # show only 2 decimals
  dplyr::mutate(dplyr::across(c("Statistic", 
                                #"p-value", 
                                "Variance explained (%)"),
                              ~ format(round(.x, 2), nsmall = 2)))

## Write to Google sheet
googlesheets4::write_sheet(
  data = sge_tbl_no,
  ss = SHEET_ID,
  sheet = "SGE_NO")


## Write final table with only significant variables

### DGE

final_dge = dplyr::bind_rows(
  list(
    "open field" = dge_tbl_of,
    "novel_object" = dge_tbl_no
    ),
  .id = "Assay") %>% 
  # filter for significant rows
  dplyr::filter(`p-value FDR-adj` < 0.05) %>% 
  # filter for states where `line` was significant
  dplyr::filter(Variable == "line") %>% 
  # remove p-value
  dplyr::select(-`p-value`) %>% 
  # convert adjusted p-value to character in scientific notation
  dplyr::mutate(`p-value FDR-adj` = as.character(scales::scientific(`p-value FDR-adj`, digits = 3))) %>% 
  # remove underscores from values
  dplyr::mutate(dplyr::across(c("Assay", "Variable"),
                              ~stringr::str_replace(., pattern = "_", " "))) %>% 
  # rename columns
  dplyr::rename("p-value (FDR-adjusted)" = "p-value FDR-adj",
                "Significance" = "Significance (p-value FDR-adj)") %>% 
  # select columns
  dplyr::select(Assay, State, `Variance explained (%)`, `p-value (FDR-adjusted)`)

## Write to Google sheet
googlesheets4::write_sheet(
  data = final_dge,
  ss = SHEET_ID,
  sheet = "DGE_FINAL")


### SGE


final_sge = dplyr::bind_rows(
  list(
    "open field" = sge_tbl_of,
    "novel_object" = sge_tbl_no
    ),
  .id = "Assay") %>% 
  # filter for significant rows
  dplyr::filter(`p-value FDR-adj` < 0.05) %>% 
  # filter for states where `test_fish` was significant
  dplyr::filter(Variable == "test_fish") %>% 
  # remove p-value
  dplyr::select(-`p-value`) %>% 
  # convert adjusted p-value to character in scientific notation
  dplyr::mutate(`p-value FDR-adj` = as.character(scales::scientific(`p-value FDR-adj`, digits = 3))) %>% 
  # remove underscores from values
  dplyr::mutate(dplyr::across(c("Assay", "Variable"),
                              ~stringr::str_replace(., pattern = "_", " "))) %>% 
  # rename columns
  dplyr::rename("p-value (FDR-adjusted)" = "p-value FDR-adj",
                "Significance" = "Significance (p-value FDR-adj)") %>% 
  # select columns
  dplyr::select(Assay, State, `Variance explained (%)`, `p-value (FDR-adjusted)`)

## Write to Google sheet
googlesheets4::write_sheet(
  data = final_sge,
  ss = SHEET_ID,
  sheet = "SGE_FINAL")


