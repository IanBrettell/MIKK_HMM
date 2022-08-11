# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(biomaRt)

# Get variables

## Debug
SIGS = "/hps/nobackup/birney/users/ian/MIKK_HMM/sig_snps/hdrr/dist_angle/dist_angle/15/5000/0.8/invnorm/open_field_sigs.csv"

## True
RES = snakemake@input[["res"]]
MIN_P = snakemake@input[["min_p"]]
SITES = snakemake@input[["sites"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
OUT_BINS_ALL = snakemake@output[["bins_all"]]
OUT_BINS_UNQ = snakemake@output[["bins_unq"]]
OUT_SNPS = snakemake@output[["snps"]]
#OUT_SNPS_UNQ = snakemake@output[["snps_unq"]]

# Read in files

sigs = readr::read_csv(SIGS,
                       col_types = c("cciiciccdddd"))


# Get annotations

## Select dataset
olat_mart = biomaRt::useEnsembl(biomart = "ensembl", dataset = "olatipes_gene_ensembl", mirror = "uswest")

olat_attr = biomaRt::listAttributes(olat_mart)

olat_genes = biomaRt::getBM(attributes = c("chromosome_name",
                                           "start_position",
                                           "end_position",
                                           "ensembl_gene_id",
                                           "hgnc_symbol",
                                           "ensembl_exon_id",
                                           "description",
                                           "strand",
                                           "transcript_start",
                                           "transcript_end"),
                            mart = olat_mart) 

olat_genes_r = olat_genes %>% 
  # change strand characters
  dplyr::mutate(strand = dplyr::recode(.$strand,
                                       `1` = "+",
                                       `-1` = "-")
  ) %>% 
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "chromosome_name",
                                          start.field = "start_position",
                                          end.field = "end_position")


# SNPs

## convert hits to genomic ranges
snp_loc_r = sigs %>% 
  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "Chr",
                                          start.field = "bp",
                                          end.field = "bp",
                                          ignore.strand = T)


## find overlaps
olaps_snps = GenomicRanges::findOverlaps(snp_loc_r, olat_genes_r)

# Pull out data frame of hits

snp_hits_all = dplyr::bind_cols(sigs %>% 
                                  dplyr::slice(olaps_snps@from),
                                olat_genes[olaps_snps@to, ])

bin_hits_unq = olat_genes[unique(olaps_snps@to), ]


## SNPs
#
### convert hits to genomic ranges
#snp_loc_r = sites %>% 
#  GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "CHROM",
#                                          start.field = "POS",
#                                          end.field = "POS",
#                                          ignore.strand = T)
#
#
### find overlaps
#olaps_snps = GenomicRanges::findOverlaps(snp_loc_r, olat_genes_r)
#
## Pull out data frame of hits
#
#snp_hits_all = dplyr::bind_cols(sites %>% 
#                                  dplyr::slice(olaps_snps@from),
#                                olat_genes[olaps_snps@to, ])
#
#snp_hits_unq = olat_genes[unique(olaps_snps@to), ]


# Write to files

readr::write_csv(bin_hits_all, OUT_BINS_ALL)
readr::write_csv(bin_hits_unq, OUT_BINS_UNQ)
readr::write_csv(sites, OUT_SNPS)
#readr::write_csv(snp_hits_unq, OUT_SNPS_UNQ)
