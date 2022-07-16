# Send log

import sys,os
import logging, traceback
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )
# Install exception handler
sys.excepthook = handle_exception


# Import libraries

import pandas as pd
import numpy as np

# Get variables
DP4 = snakemake.input.dp4[0]
SITES = snakemake.input.sites_file
OUT = snakemake.output[0]

# Read in dp4 file

dp4_cols = ["CHR", "POS", "REF", "TOTAL", "A", "C", "G", "T", "N"]
dp4 = pd.read_csv(DP4, sep = "\t", names = dp4_cols, index_col = ["CHR", "POS"])

# Read in sites file

sites_cols = ["CHR", "POS", "POS_2", "REF", "ALT", "F0_1", "F0_2"]
sites = pd.read_csv(SITES, sep = "\t", names = sites_cols, index_col = ["CHR", "POS"])

# Add columns with F0 and F1 alleles

sites = sites.assign(F0_1_ALLELE =lambda x: np.where(x["F0_1"] == "0/0", x["REF"], x["ALT"]))
sites = sites.assign(F0_2_ALLELE =lambda x: np.where(x["F0_2"] == "0/0", x["REF"], x["ALT"]))

# Filter dp4 for necessary columns

dp4 = dp4[["A", "C", "G", "T"]]

# Join dp4 to sites

joined = dp4.join(sites)

# Create new columns with counts for F0_1 and F0_2 alleles

def rules(row, column):
    if row[column] == "A":
        return row["A"]
    elif row[column] == "C":
        return row["C"]
    elif row[column] == "G":
        return row["G"]            
    elif row[column] == "T":
        return row["T"]               
    else:
        return None
joined['F0_1_COUNT'] = joined.apply(lambda x: rules(x, "F0_1_ALLELE"), 1)
joined['F0_2_COUNT'] = joined.apply(lambda x: rules(x, "F0_2_ALLELE"), 1)

# Drop NAs

joined = joined.dropna()

# Convert F0_1_COUNT and F0_2_COUNT to integers

joined['F0_1_COUNT'] = joined['F0_1_COUNT'].astype(int)
joined['F0_2_COUNT'] = joined['F0_2_COUNT'].astype(int)

# Write to file

joined[["F0_1_ALLELE", "F0_1_COUNT", "F0_2_ALLELE", "F0_2_COUNT"]].to_csv(OUT, sep = "\t", header = False)
