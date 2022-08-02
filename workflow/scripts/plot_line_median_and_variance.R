# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)

# Variables

## Debug
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out_split/0.05/dist_angle/15/F0.csv"
INTERVAL = 0.05
SELECTED_LINES = list("8-2", "18-2", "21-2", "38-2", "40-1", "50-2") %>% 
  unlist()
LINE_COLS = "config/line_colours/line_colours_0.05.csv"

## True
IN = snakemake@input[["dat"]]
INTERVAL = snakemake@params[["interval"]] %>% 
  as.numeric()
SELECTED_LINES = snakemake@params[["selected_lines"]] %>% 
  unlist()
LINE_COLS = snakemake@input[["line_cols"]]
PNG_ALL = snakemake@output[["png_all"]]
PDF_ALL = snakemake@output[["pdf_all"]]
PNG_SELECT = snakemake@output[["png_select"]]
PDF_SELECT = snakemake@output[["pdf_select"]]

# Read in file

raw = readr::read_csv(IN) %>% 
  # convert `time` to character and add a 0 if only 3 characters
  dplyr::mutate(time = as.character(time),
                time = dplyr::if_else(nchar(time) == 3,
                                      paste("0", time, sep = ""),
                                      time))

# Read in line colours
line_cols = readr::read_csv(LINE_COLS)

# Process

df = raw %>% 
  # exclude iCabs when paired with a different test fish
  dplyr::filter(!(fish == "ref" & test_fish != "iCab")) %>% 
  # rename outbred
  dplyr::mutate(test_fish = dplyr::if_else(stringr::str_detect(test_fish, "outbred"),
                                           "Kiyosu CC",
                                           test_fish)) %>% 
  # create `indiv` column
  tidyr::unite(date, time, quadrant, fish,
               col = "indiv",
               sep = "_",
               remove = F) %>% 
  # group by individual
  dplyr::group_by(indiv, test_fish) %>% 
  # calculate mean speed
  dplyr::summarise(mean_speed = mean(distance)) %>% 
  dplyr::ungroup() %>% 
  # get each line's median and variance for individual mean speed
  dplyr::group_by(test_fish) %>% 
  dplyr::summarise(median_speed = median(mean_speed),
                   variance_speed = var(mean_speed))

# Get palette for all MIKK lines


PAL = line_cols$colour
names(PAL) = line_cols$line

# Plot all lines

fig_all = df %>% 
  ggplot(aes(median_speed, variance_speed, colour = test_fish)) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = test_fish),
                           size = 3
  ) +
  scale_colour_manual(values = PAL) +
  coord_flip() +
  guides(colour = "none") +
  cowplot::theme_cowplot() +
  xlab(paste("line median for individual mean speed\n(pixels per", INTERVAL, "seconds)", sep = " ")) +
  ylab(paste("line variance for individual mean speed\n(pixels per", INTERVAL, "seconds)", sep = " ")) 

ggsave(PNG_ALL,
       fig_all,
       device = "png",
       width = 7,
       height = 7,
       units = "in",
       dpi = 400)

ggsave(PDF_ALL,
       fig_all,
       device = "pdf",
       width = 7,
       height = 7,
       units = "in",
       dpi = 400)
  
# Highlighting selected lines

# Add column based on whether it was selected for the F2 cross
df = df %>% 
  dplyr::left_join(line_cols,
                   by = c("test_fish" = "line")) %>% 
  dplyr::mutate(selected = dplyr::if_else(test_fish %in% SELECTED_LINES,
                                          colour,
                                          "dimgrey"))

## Plot

fig_select = df %>% 
  ggplot(aes(median_speed, variance_speed)) + 
  geom_point(colour = df$selected) +
  ggrepel::geom_text_repel(aes(label = test_fish),
                           size = 3,
                           colour = df$selected
  ) +
  coord_flip() +
  guides(colour = "none") +
  cowplot::theme_cowplot() +
  xlab(paste("line median for individual mean speed\n(pixels per", INTERVAL, "seconds)", sep = " ")) +
  ylab(paste("line variance for individual mean speed\n(pixels per", INTERVAL, "seconds)", sep = " ")) 


ggsave(PNG_SELECT,
       fig_select,
       device = "png",
       width = 7,
       height = 7,
       units = "in",
       dpi = 400)

ggsave(PDF_SELECT,
       fig_select,
       device = "pdf",
       width = 7,
       height = 7,
       units = "in",
       dpi = 400)

