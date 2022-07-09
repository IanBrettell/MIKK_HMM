# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)

# Variables

## Debug
IN = "/hps/nobackup/birney/users/ian/MIKK_F0_tracking/merged/0.08.csv"
INTERVAL = 0.08

## True
IN = snakemake@input[[1]]
INTERVAL = snakemake@params[["interval"]] %>% 
  as.numeric()
OUT_FIG = snakemake@output[["fig"]]
OUT_CSV = snakemake@output[["csv"]]

# Read in file

raw = readr::read_csv(IN)

# Process

df = raw %>% 
  # exclude iCabs when paired with a different test fish
  dplyr::filter(!(fish == "ref" & test_fish != "iCab")) %>% 
  # rename outbred
  dplyr::mutate(test_fish = dplyr::if_else(stringr::str_detect(test_fish, "outbred"),
                                           "outbred",
                                           test_fish)) %>% 
  # create `indiv` column
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv",
               sep = "_",
               remove = F) %>% 
  # group by individual
  dplyr::group_by(indiv, test_fish) %>% 
  # calculate mean speed
  dplyr::summarise(mean_speed = mean(distance))

# Rank by line

df_rank = df %>% 
  dplyr::group_by(test_fish) %>% 
  dplyr::summarise(median_speed = median(mean_speed)) %>% 
  dplyr::arrange(median_speed)

# Get palette for all MIKK lines

## How many lines excluding iCab, outbred, and luzonensis
EXCL = c("iCab", "luzonensis", "outbred")
MIKK_RANK = df_rank %>% 
  dplyr::filter(!test_fish %in% EXCL) %>% 
  dplyr::pull(test_fish)
N_MIKK = length(MIKK_RANK)
PAL = rev(scales::hue_pal()(N_MIKK))
names(PAL) = MIKK_RANK
# Add greys to PAL
EXTRA = c("#424B4D", "#C8CED0", "#718084")
names(EXTRA) = EXCL
PAL = c(EXTRA["iCab"], PAL, EXTRA["outbred"], EXTRA["luzonensis"])

# Convert `test_fish` to factor to order
df = df %>% 
  dplyr::mutate(test_fish = factor(test_fish, levels = names(PAL)))

# Plot

out_plot = df %>% 
  ggplot(aes(test_fish, mean_speed, fill = test_fish)) + 
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(size = 0.25) +
  scale_fill_manual(values = PAL) +
  coord_flip() +
  guides(fill = "none") +
  cowplot::theme_cowplot() +
  scale_x_discrete(limits = rev(levels(df$test_fish))) +
  xlab("line") +
  ylab(paste("mean speed (pixels per", INTERVAL, "seconds)", sep = " "))
  
ggsave(OUT_FIG,
       out_plot,
       device = "png",
       width = 8,
       height = 11,
       units = "in",
       dpi = 400)

# Save colours to .csv
pal_df = tibble::tibble("line" = names(PAL),
                        "colour" = PAL)
readr::write_csv(pal_df, OUT_CSV)
