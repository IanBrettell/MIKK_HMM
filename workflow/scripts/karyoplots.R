# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(karyoploteR)
library(ggbeeswarm)

# Get variables

## Debug
BIN_LENGTH = as.numeric("5000")
IN_FILE = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out/F2/hdrr/hmmlearn_true/5000/0.8.csv"
COV = as.numeric("0.8")
LINE_COLS = here::here("config/line_colours/line_colours_0.08.csv")

## True
IN_FILE = snakemake@input[["data"]]
LINE_COLS = snakemake@input[["line_cols"]]
COV = snakemake@params[["cov"]]
BIN_LENGTH = as.numeric(snakemake@params[["bin_length"]])
KARYOPLOT_NOMISS = snakemake@output[["karyoplot_no_missing"]]
KARYOPLOT_WIMISS = snakemake@output[["karyoplot_with_missing"]]


######################
# Number of states
######################

N_STATES = 3


######################
# Read in data
######################

df = readr::read_csv(IN_FILE,
                     col_types = "ciiiiidi") %>% 
  # add key variables
  dplyr::mutate(BIN_START = (BIN * BIN_LENGTH) + 1,
                BIN_END = ((BIN + 1) * BIN_LENGTH)) %>% 
  # split `SAMPLE` into `SAMPLE`, `PAT` and `MAT`
  tidyr::separate(col = "SAMPLE",
                  into = c("LINE_ID", "PAT", "MAT"),
                  sep = "_",
                  remove = F) %>% 
  # rename `SAMPLE` and `LINE_ID`
  dplyr::rename(SAMPLE_PAT_MAT = SAMPLE,
                SAMPLE = LINE_ID) %>% 
  # convert `SAMPLE` to numeric and sort
  dplyr::mutate(SAMPLE = as.integer(SAMPLE)) %>% 
  dplyr::arrange(SAMPLE, CHROM, BIN)

# Recode vector

recode_vec = c(`0` = "Homozygous pat",
               `1` = "Heterozygous",
               `2` = "Homozygous mat")

# Recode `STATE`

df = df %>% 
  dplyr::mutate(STATE_RECODE = dplyr::recode(STATE,
                                             !!!recode_vec))

# Read in line colours
line_cols = readr::read_csv(LINE_COLS)

# Add `COLOUR` column based on genotype and pat/mat line

df = df %>% 
  dplyr::left_join(line_cols %>% 
              dplyr::rename(PAT_COL = colour),
            by = c("PAT" = "line")) %>% 
  dplyr::left_join(line_cols %>% 
              dplyr::rename(MAT_COL = colour),
            by = c("MAT" = "line")) %>% 
  dplyr::mutate(COLOUR = dplyr::case_when(STATE_RECODE == "Homozygous pat" ~ PAT_COL,
                                          STATE_RECODE == "Homozygous mat" ~ MAT_COL,
                                          STATE_RECODE == "Heterozygous" ~ "#2B303A"))


######################
# Palette and plotting params
######################

# Set states to loop over

states = 0:(N_STATES - 1)

# Read in total medaka genome count

## Get chromosome lengths
med_chr_lens = readr::read_csv(here::here("config/hdrr_chrom_lengths.csv"),
                               col_names = c("chr", "end"))
## Add start
med_chr_lens$start = 1
## Reorder
med_chr_lens = med_chr_lens %>% 
  dplyr::select(chr, start, end) %>% 
  # remove MT
  dplyr::filter(chr != "MT")
## Total HdrR sequence length
total_hdrr_bases = sum(med_chr_lens$end)

# Get number of samples (for setting the height of the Karyoplots)

N_SAMPLES = unique(df$SAMPLE) %>% 
  length()

#######################
# Karyoplot no missing
#######################

# Create custom genome 

med_genome = regioneR::toGRanges(as.data.frame(med_chr_lens))

# Convert data to list of block boundaries for each LANE

block_bounds_list = df %>% 
  # loop over LANE
  split(., f = .$SAMPLE) %>% 
  purrr::map(., function(SAMPLE){
    # loop over CHR
    SAMPLE %>% 
      split(., f = .$CHROM) %>% 
      purrr::map(., function(CHROM){
        # Get lengths of each contiguous state
        cont_len = rle(CHROM$STATE)
        
        # Get cumulative sum of those lengths
        cum_blocks = cumsum(cont_len$lengths)
        
        # Get rows that correspond to block changes
        block_bounds = CHROM[cum_blocks, ] %>% 
          # Add end of previous block
          dplyr::mutate(END_PREV = dplyr::lag(BIN_END)) %>% 
          # Replace the NA in the first row with `1`
          dplyr::mutate(END_PREV = tidyr::replace_na(END_PREV, 1)) 
        
      }) %>% 
      dplyr::bind_rows()
    
  })

# Extract y cutoff points for each lane

lane_cutoffs = cut(0:1, breaks = length(block_bounds_list), dig.lab = 7) %>% 
  levels(.) %>% 
  data.frame(lower = as.numeric( sub("\\((.+),.*", "\\1", .) ),
             upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", .) )) %>% 
  dplyr::arrange(dplyr::desc(lower))

# Plot karyoplot WITH NO missing blocks

png(file=KARYOPLOT_NOMISS,
    width=7800,
    height=round(37.6*N_SAMPLES),
    units = "px",
    res = 400)

# Plot ideogram
kp = karyoploteR::plotKaryotype(med_genome, plot.type = 5)

# Plot title
karyoploteR::kpAddMainTitle(kp,
                            paste("Emission (co)variances: ", 
                                  COV,
                                  "\nBin length: ",
                                  BIN_LENGTH,
                                  sep = ""),
                            cex=4)

# Add rectangles in loop
counter = 0
purrr::map(block_bounds_list, function(SAMPLE){
  # Add to counter_B
  counter <<- counter + 1
  # Add rectangles
  karyoploteR::kpRect(kp,
                      chr = SAMPLE$CHROM,
                      x0 = SAMPLE$END_PREV,
                      x1 = SAMPLE$BIN_END,
                      y0 = lane_cutoffs[counter, ] %>% 
                        dplyr::pull(lower),
                      y1 = lane_cutoffs[counter, ] %>% 
                        dplyr::pull(upper),
                      col = SAMPLE$COLOUR,
                      border = NA)
  # Add axis label
  karyoploteR::kpAddLabels(kp, labels = unique(SAMPLE$SAMPLE_PAT_MAT),
                           r0 = lane_cutoffs[counter, ] %>% 
                             dplyr::pull(lower),
                           r1 = lane_cutoffs[counter, ] %>% 
                             dplyr::pull(upper),
                           cex = 0.5)
})


dev.off()

#######################
# Karyoplot WITH missing
#######################

## Process
block_bounds_list = df %>% 
  # loop over LANE
  split(., f = .$SAMPLE) %>% 
  purrr::map(., function(SAMPLE){
  
    STRAIN = unique(SAMPLE$SAMPLE)
    # Create list of possible bins
    poss_bins = purrr::map(med_chr_lens$chr, function(CHROM){
      # Get chr end
      CHR_END = med_chr_lens %>% 
        dplyr::filter(chr == CHROM) %>% 
        dplyr::pull(end) %>% 
        as.numeric()
      # Get bin starts
      out = tibble::tibble(CHROM = as.numeric(CHROM),
                           BIN_START = seq(from = 1, to = CHR_END, by = BIN_LENGTH),
                           BIN_END = BIN_START + BIN_LENGTH - 1
      )
      # Adjust final bin end 
      out[nrow(out), "BIN_END"] = CHR_END
      
      return(out)
    }) %>% 
      dplyr::bind_rows()
  
    
    # Bind DF
    new_df = dplyr::left_join(poss_bins,
                              SAMPLE %>% 
                                dplyr::select(SAMPLE_PAT_MAT, CHROM, BIN_START, BIN_END, STATE, COLOUR),
                              by = c("CHROM", "BIN_START", "BIN_END")) %>% 
      # replace NAs with `UNCLASSIFIED`
      dplyr::mutate(STATE = as.character(STATE),
                    STATE = STATE %>% 
                      tidyr::replace_na("UNCLASSIFIED"),
                    # add STRAIN
                    SAMPLE = STRAIN)
  
            
  })

# Plot karyoplot

png(file=KARYOPLOT_WIMISS,
    width=7800,
    height=round(37.6*N_SAMPLES),
    units = "px",
    res = 400)

# Plot ideogram
kp = karyoploteR::plotKaryotype(med_genome, plot.type = 5)

# Plot title
karyoploteR::kpAddMainTitle(kp,
                            paste("Emission (co)variances: ", 
                                  COV,
                                  "\nBin length: ",
                                  BIN_LENGTH,
                                  sep = ""),
                            cex=4)


# Add rectangles in loop
counter = 0
purrr::map(block_bounds_list, function(SAMPLE){
  # Add to counter
  counter <<- counter + 1
  # Add rectangles
  karyoploteR::kpRect(kp,
                      chr = SAMPLE$CHROM,
                      x0 = SAMPLE$BIN_START,
                      x1 = SAMPLE$BIN_END,
                      y0 = lane_cutoffs[counter, ] %>% 
                        dplyr::pull(lower),
                      y1 = lane_cutoffs[counter, ] %>% 
                        dplyr::pull(upper),
                      col = SAMPLE$COLOUR,
                      border = NA)
  # Add axis label
  karyoploteR::kpAddLabels(kp, labels = unique(SAMPLE$SAMPLE_PAT_MAT),
                           r0 = lane_cutoffs[counter, ] %>% 
                             dplyr::pull(lower),
                           r1 = lane_cutoffs[counter, ] %>% 
                             dplyr::pull(upper),
                           cex = 0.5)
})


dev.off()  

