# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(cowplot)

# Set variables

## Debug
IN = list(list.files(here::here("book/figs/gcta/hdrr/5000/0.8/dge/invnorm/open_field"),
                     pattern = "None.png",
                     full.names = T)) %>% 
  unlist()
SHEET_ID = "118HNbI7ch_0rgRSaDB1b73mpRo5Z2g7uQvNhjZUD7WI"
ASSAY = "open_field"
DGE_SGE = "dge" %>% 
  toupper()
COVARS = "None"

############################
# Get significant states
############################

if (DGE_SGE == "DGE"){
  if (ASSAY == "open_field"){
    SIGS = googlesheets4::read_sheet(SHEET_ID, sheet = "DGE_OF") %>% 
      dplyr::filter(`p-value FDR-adj` < 0.05 & `Variable` == "line") %>% 
      dplyr::pull(State) %>% 
      as.integer()
  }
  else if (ASSAY == "novel_object"){
    SIGS = googlesheets4::read_sheet(SHEET_ID, sheet = "DGE_NO") %>% 
      dplyr::filter(`p-value FDR-adj` < 0.05 & `Variable` == "line") %>% 
      dplyr::pull(State)%>% 
      as.integer()
  }
} else if (DGE_SGE == "SGE"){
  if (ASSAY == "open_field"){
    SIGS = googlesheets4::read_sheet(SHEET_ID, sheet = "SGE_OF") %>% 
      dplyr::filter(`p-value FDR-adj` < 0.05 & `Variable` == "test_fish") %>% 
      dplyr::pull(State)%>% 
      as.integer()
  } else if (ASSAY == "novel_object"){
    SIGS = googlesheets4::read_sheet(SHEET_ID, sheet = "SGE_NO") %>% 
      dplyr::filter(`p-value FDR-adj` < 0.05 & `Variable` == "test_fish") %>% 
      dplyr::pull(State)%>% 
      as.integer()
  }
}

############################
# Filter files for significant states
############################

in_files = tibble::tibble(PATH = IN) %>% 
  dplyr::mutate(NAME = basename(IN) %>% 
                  stringr::str_remove(".png")) %>% 
  tidyr::separate(col = NAME,
                  into = c("STATE", "COVAR"),
                  sep = "_") %>% 
  # convert STATE to numeric
  dplyr::mutate(STATE = as.numeric(STATE)) %>% 
  # filter for significant states
  dplyr::filter(STATE %in% SIGS) %>% 
  # arrange
  dplyr::arrange(STATE) %>% 
  # pull out paths to files
  dplyr::pull(PATH)


############################
# Compile into single file
############################

img = magick::image_read(in_files) # read from vector of paths
img2 = magick::image_append(img, stack = TRUE) # places pics above one another
test = magick::image_read(in_files[1])
magick::image_write(test, format = "pdf", path = here::here("tmp.pdf"))
