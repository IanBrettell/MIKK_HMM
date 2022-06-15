# Send stdout and stderr to log file
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

# Get variables

## Debug
IN = [
    "/hps/nobackup/birney/users/ian/MIKK_all/0.08/F0.csv",
    "/hps/nobackup/birney/users/ian/MIKK_all/0.08/F2.csv",
    "/hps/nobackup/birney/users/ian/MIKK_all/0.08/Kiyosu_CC.csv"
     ]

## True
IN = snakemake.input
