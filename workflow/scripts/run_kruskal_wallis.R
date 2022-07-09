# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
#IN = list("/hps/nobackup/birney/users/ian/pilot/hmm_out/0.5/dist_angle/10.csv",
#          "/hps/nobackup/birney/users/ian/pilot/hmm_out/0.5/dist_angle/15.csv")

## True
IN = snakemake@input
OUT = snakemake@output[[1]]

# Split parameters from file names

filenames = tibble(FILE = unlist(IN)) %>% 
  tidyr::separate(.,
                  FILE,
                  into = c(rep(NA, 8), "INTERVAL", "VARIABLES", "N_STATES"),
                  sep = "/") %>% 
  # remove file ext
  dplyr::mutate(N_STATES = stringr::str_remove(N_STATES, ".csv"))

# Read in files

kw_tests = purrr::map(IN, readr::read_csv) %>% 
  # loop over each param combination
  purrr::map(., function(PARAM_DF){
    out = list()
    
    # Get frequencies of each state per fish
    
    out[["STATE_FREQ_DF"]] = PARAM_DF %>% 
      # remove iCab when paired with a different test fish
      dplyr::filter(!(fish == "ref" & ref_fish == "icab" & test_fish != "icab")) %>% 
      # add `line` %>% 
      dplyr::mutate(line = dplyr::case_when(fish == "ref" ~ ref_fish,
                                            fish == "test" ~ test_fish)) %>% 
      # get proportions of time spent in each state
      ## count rows per fish per state
      dplyr::count(assay, date, time, quadrant, fish, line, state) %>% 
      # add total row count per fish
      dplyr::add_count(assay, date, time, quadrant, fish, line, wt = n, name = "nn") %>% 
      # get proportion of time fish spent in each state
      dplyr::mutate(state_freq = n / nn)
    
    # Run KS test per state
    
    out[["KW_PER_STATE"]] = out[["STATE_FREQ_DF"]] %>% 
      dplyr::group_split(state) %>% 
      purrr::map_dfr(., function(DF) {
        broom::tidy(kruskal.test(x = DF$state_freq, g = DF$line))
      }) 
    
    # Sum statistic
    
    out[["SUM_KW_STAT"]] = sum(out[["KW_PER_STATE"]]$statistic)
    
    return(out)
  })


# Bind into single DF

out = dplyr::bind_cols(filenames,
                       SUM_KW_STAT = purrr::map_dbl(kw_tests, "SUM_KW_STAT"))

# Save to file

saveRDS(out, OUT)
