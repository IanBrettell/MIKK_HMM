# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

IN_PED = list.files("/hps/nobackup/birney/users/ian/MIKK_HMM/peds_contigs/F2/hdrr/5000/0.8",
                    full.names = T,
                    pattern = ".ped")

IN_MAP = list.files("/hps/nobackup/birney/users/ian/MIKK_HMM/peds_contigs/F2/hdrr/5000/0.8",
                    full.names = T,
                    pattern = ".map")

## True
IN_PED = snakemake@input[["ped"]] %>% 
  unlist()
IN_MAP = snakemake@input[["map"]] %>% 
  unlist()

# Ensure the are in numeric order

chroms_p = IN_PED %>% 
  basename(.) %>% 
  stringr::str_remove(".ped") %>% 
  as.numeric()

IN_PED = IN_PED[order(chroms_p)]

chroms_m = IN_MAP %>% 
  basename(.) %>% 
  stringr::str_remove(".map") %>% 
  as.numeric()

IN_MAP = IN_MAP[order(chroms_m)]

# Read in

counter = 0
out_p = lapply(IN_PED[c(22,23,24)], function(PED){
  out = as.matrix(readr::read_tsv(PED,col_names = F))
  out = out[-1, ]
})

