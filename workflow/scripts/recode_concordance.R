# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)

# Get variables

## Debug
#IN_FILE = "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_out/0.05/dist_angle/17.csv"
#N_STATES =17
#INTERVAL = 0.05

## True
IN_FILE = snakemake@input[[1]]
OUT_RDS = snakemake@output[["rds"]]
OUT_CONFM_PNG = snakemake@output[["confmat_png"]]
OUT_CONFM_PDF = snakemake@output[["confmat_pdf"]]
OUT_POLAR_PNG = snakemake@output[["polar_png"]]
OUT_POLAR_PDF = snakemake@output[["polar_pdf"]]
N_STATES = snakemake@params[["n_states"]] %>% 
  as.numeric()
INTERVAL = snakemake@params[["interval"]] %>% 
  as.numeric()


#######################
# Read in data
#######################

# Get number of rows (for plotting) based on number of states

N_ROWS = N_STATES / 5

# Read in file

df = readr::read_csv(IN_FILE) %>% 
  # recode angle to sit between 0 and 360
  dplyr::mutate(angle_recode = ifelse(angle < 0,
                                      180 + (180 + angle),
                                      angle))

#######################
# Confusion matrix
#######################

# Rank states by highest overlap between train_self and train_other (confusion matrix)

# Split by group

df_list = df %>% 
  split(., ~ group)

# Match states and plot confusion matrix

conf_list = purrr::imap(df_list, function(GROUP, NAME){
  # set up list to direct output
  res = list()
  
  # Get order of most populous states
  self_states_rank = GROUP %>% 
    dplyr::count(train_self) %>% 
    dplyr::arrange(dplyr::desc(n)) %>% 
    # add column with ranked state
    dplyr::mutate(train_self_rc = 0:(n() - 1))
  
  self_state_rc = self_states_rank$train_self_rc
  names(self_state_rc) = self_states_rank$train_self
  
  # recode train_self in GROUP to ranked state
  GROUP = GROUP %>% 
    dplyr::mutate(train_self = dplyr::recode(train_self, !!!self_state_rc))
  
  # Loop over most populous states
  ## Create empty vector for output
  vec = double()
  STATES = 0:(N_STATES - 1)
  out = purrr::imap(STATES, function(STATE, COUNTER){
    # Set initial target row as 1
    target_row = 1
    
    res = GROUP %>% 
      dplyr::filter(train_self == STATE) %>% 
      dplyr::count(train_self, train_other) %>% 
      # sort
      dplyr::arrange(dplyr::desc(n)) %>% 
      # pull out first row
      dplyr::slice(target_row) %>% 
      # get train_other 
      dplyr::pull(train_other)
    
    while (res %in% vec){
      target_row = target_row + 1
      res = GROUP %>% 
        dplyr::filter(train_self == STATE) %>% 
        dplyr::count(train_self, train_other) %>% 
        # sort
        dplyr::arrange(dplyr::desc(n)) %>% 
        # pull out first row
        dplyr::slice(target_row) %>% 
        # get train_other 
        dplyr::pull(train_other)
      
      # If there are no more states with overlapping states that have not already been taken,
      # Take the first unused state
      if (length(res) == 0){
        res = STATES[which(!STATES %in% vec)][1]
      }
      
    }
    
    vec[COUNTER] <<- res
    
  })
  
  
  # combine into recode vector
  recode_vec = self_state_rc
  names(recode_vec) = vec
  
  # recode
  
  GROUP = GROUP %>% 
    dplyr::mutate(train_other = dplyr::recode(train_other, !!!recode_vec))
  
  # concordance
  
  concordance = length(which(GROUP$train_self == GROUP$train_other)) / nrow(GROUP)
  
  # Plot confusion matrix
  
  conf_mat = GROUP %>% 
    dplyr::count(train_self, train_other) %>% 
    ggplot() +
    geom_tile(aes(train_self, train_other, fill = n)) +
    scale_fill_viridis_c(option = "plasma") +
    theme_bw() +
    theme(aspect.ratio = 1) +
    xlab("self-trained state") +
    ylab("other-trained state (recoded)") +
    ggtitle(paste("Group: ",
                  NAME,
                  "\nInterval: ",
                  INTERVAL,
                  " seconds",
                  "\nN states: ",
                  N_STATES,
                  "\nConcordance: ",
                  concordance,
                  sep = ""))
  
  res[["data"]] = GROUP
  res[["recode_vec"]] = recode_vec
  res[["concordance"]] = concordance
  res[["conf_mat"]] = conf_mat
  
  return(res)
  
})

# Compile plots

final_cmats = cowplot::plot_grid(conf_list$A$conf_mat,
                                 conf_list$B$conf_mat,
                                 nrow = 1,
                                 ncol = 2,
                                 axis = "tblr")

# Save .png
ggsave(OUT_CONFM_PNG,
       final_cmats,
       device = "png",
       width = 9.708,
       height = 6,
       units = "in",
       dpi = 400)

# Save .pdf
ggsave(OUT_CONFM_PDF,
       final_cmats,
       device = "pdf",
       width = 9.708,
       height = 6,
       units = "in",
       dpi = 400)

#######################
# Eye of Cabs
#######################

# Plot

## Split by group/training_data combinations
df_split = purrr::map(conf_list, "data") %>% 
  # bind rows
  dplyr::bind_rows(.) %>% 
  # pivot longer
  tidyr::pivot_longer(cols = c(train_self, train_other),
                      names_to = "training_data",
                      names_prefix = "train_",
                      values_to = "state") %>% 
  # order `training_data`
  dplyr::mutate(training_data = factor(training_data, levels = c("self", "other"))) %>% 
  split(., ~ group + training_data)


out_plots = purrr::imap(df_split, function(DATA, NAME){
  
  # get plot title 
  title = NAME %>% 
    stringr::str_split(pattern = "\\.", simplify = T)
  title = paste("Group: ", title[1], "\nTraining data: ", title[2], sep = "")
  # Plot
  out_plot = DATA %>% 
    # select random sample of 1e5 rows
    dplyr::slice_sample(n = 1e5) %>% 
    ggplot() +
    geom_point(aes(angle_recode, log10(distance), colour = state),
               alpha = 0.3, size = 0.2) +
    coord_polar() +
    facet_wrap(~state, nrow = N_ROWS) +
    theme_dark(base_size = 8) +
    scale_x_continuous(labels = c(0, 90, 180, 270),
                       breaks = c(0, 90, 180, 270)) +
    scale_color_viridis_c() +
    guides(colour = "none") +
    xlab("angle of travel") +
    ylab(expression(log[10]("distance travelled in pixels"))) +
    ggtitle(title)
  
  return(out_plot)
})

# Compile 4 plots into single plot

final_plot = cowplot::plot_grid(plotlist = out_plots,
                                nrow = 2,
                                ncol = 2,
                                axis = "tblr")

# Save .png
ggsave(OUT_POLAR_PNG,
       final_plot,
       device = "png",
       width = 19.416,
       height = 12,
       units = "in",
       dpi = 400)

# Save .pdf
ggsave(OUT_POLAR_PDF,
       final_plot,
       device = "pdf",
       width = 19.416,
       height = 12,
       units = "in",
       dpi = 400)


#######################
# Save concordance metrics
#######################

out_df = tibble(INTERVAL = INTERVAL,
                N_STATES = N_STATES,
                GROUP = purrr::imap_chr(conf_list, ~.y),
                CONCORDANCE = purrr::map_dbl(conf_list, "concordance"))

saveRDS(out_df, OUT_RDS)
