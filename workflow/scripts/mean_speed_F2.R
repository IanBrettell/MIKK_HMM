# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out_split/0.05/dist_angle/15/F2.csv"
DGE_SGE = "dge"
INTERVAL = "0.05"

## True
IN = snakemake@input[["data"]]
OUT_CSV_NT = snakemake@output[["csv_notrans"]]
OUT_CSV_IV = snakemake@output[["csv_invnorm"]]
OUT_HIST = snakemake@output[["hist"]]
DGE_SGE = snakemake@params[["dge_sge"]]

# Get output directories for .csv files

OUT_DIR_NT_OF = dirname(OUT_CSV_NT[["notrans_of"]])
OUT_DIR_NT_NO = dirname(OUT_CSV_NT[["notrans_no"]])
OUT_DIR_IV_OF = dirname(OUT_CSV_IV[["invnorm_of"]])
OUT_DIR_IV_NO = dirname(OUT_CSV_IV[["invnorm_no"]])

# Read data

raw = readr::read_csv(IN)

# Add inverse-normalisation function
invnorm = function(x) {
  res = rank(x)
  res = qnorm(res/(length(res)+0.5))
  return(res)
}

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
  # convert `date` to factor
  dplyr::mutate(date = factor(date))

#######################
## DGE
#######################

if (DGE_SGE == "dge"){
  # Get proportion of time each fish spent in each state
  df_dge = df %>% 
    # remove iCab 
    dplyr::filter(fish == "test") %>% 
    # group by individual
    dplyr::group_by(indiv, assay, test_fish) %>% 
    # calculate mean speed
    dplyr::summarise(mean_speed = mean(distance)) %>% 
    # ungroup
    dplyr::ungroup()
  
  # Write to file
  readr::write_csv(df_dge, OUT_PATH)
  
  # Split by assay
  
  dge_hist_pretrans = df_dge %>% 
    ggplot() + 
    geom_histogram(aes(mean_speed),
                   bins = 40,
                   fill = "#EAB464",
                   colour = "#A7754D") +
    facet_grid(cols = vars(assay)) +
    cowplot::theme_cowplot() +
    scale_fill_viridis_d() +
    guides(fill = "none")+
    xlab(paste("mean distance in pixels\nbetween",
               INTERVAL,
               "second intervals"))
  
  ### Inverse-normalise
  
  df_dge = df_dge %>% 
    dplyr::group_by(assay) %>% 
    dplyr::mutate(mean_speed = invnorm(mean_speed)) %>% 
    dplyr::ungroup()
  
  # Write to file
  readr::write_csv(DF, OUT_PATH)

  
  # Plot
  
  dge_hist_posttrans = df_dge %>% 
    ggplot() + 
    geom_histogram(aes(mean_speed),
                   bins = 40,
                   fill = "#EAB464",
                   colour = "#A7754D") +
    facet_grid(cols = vars(assay)) +
    cowplot::theme_cowplot() +
    scale_fill_viridis_d() +
    guides(fill = "none")+
    xlab(paste("mean distance in pixels\nbetween",
               INTERVAL,
               "second intervals"))
  
  # Compile into single plot
  
  dge_hist = cowplot::plot_grid(dge_hist_pretrans,
                                dge_hist_posttrans,
                                align = "hv",axis = "tblr",
                                labels = c("A", "B"))
  
  ggsave(OUT_HIST,
         dge_hist,
         device = "png",
         width = 11,
         height = 6,
         units = "in",
         dpi = 400)
}

#######################
## SGE
#######################

if (DGE_SGE == "sge"){
  # Get proportion of time each fish spent in each state
  df_sge = df %>% 
    # remove iCab 
    dplyr::filter(fish == "ref") %>% 
    # group by individual
    dplyr::group_by(indiv, assay, test_fish) %>% 
    # calculate mean speed
    dplyr::summarise(mean_speed = mean(distance)) %>% 
    # ungroup
    dplyr::ungroup()
  
  # Write to file
  readr::write_csv(DF, OUT_PATH)

  
  sge_hist_pretrans = df_sge %>% 
    ggplot() + 
    geom_histogram(aes(mean_speed),
                   bins = 40,
                   fill = "#934B00",
                   colour = "#690500") +
    facet_grid(cols = vars(assay)) +
    cowplot::theme_cowplot() +
    scale_fill_viridis_d() +
    guides(fill = "none")+
    xlab(paste("mean distance in pixels\nbetween",
               INTERVAL,
               "second intervals"))
  
  
  ### Inverse-normalise
  
  df_sge = df_sge %>% 
    dplyr::group_by(assay) %>% 
    dplyr::mutate(mean_speed = invnorm(mean_speed)) %>% 
    dplyr::ungroup()
  
  # Write to file
  readr::write_csv(DF, OUT_PATH)

  sge_hist_posttrans = df_sge %>% 
    ggplot() + 
    geom_histogram(aes(mean_speed),
                   bins = 40,
                   fill = "#934B00",
                   colour = "#690500") +
    facet_grid(cols = vars(assay)) +
    cowplot::theme_cowplot() +
    scale_fill_viridis_d() +
    guides(fill = "none")+
    xlab(paste("mean distance in pixels\nbetween",
               INTERVAL,
               "second intervals"))
  
  # Compile into single plot
  
  sge_hist = cowplot::plot_grid(sge_hist_pretrans,
                                sge_hist_posttrans,
                                align = "hv",axis = "tblr")
  
  ggsave(OUT_HIST,
         sge_hist,
         device = "png",
         width = 11,
         height = 6,
         units = "in",
         dpi = 400)
}



