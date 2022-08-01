# Notes

## 1 August 2022

Discovered that GCTA outputs a different GRM matrix than the one computed manually, and looks incorrect.

The .ped with all SNPs contains 20,796,377 SNPs.

## 29 July 2022

SNP count of .ped with all HOM SNPs: 
    - 21,294,644
SNP count of .ped with all HOM SNPs, excluding SNPs with any missing calls:
    - 21,294,644

## 18 July 2022

Best combination of HMM parameters is:
    - 0.05 interval
    - 15 states

## 15 July 2022

* `config/F2_samples_converted.csv` was copied from the repo/path `MIKK_F2_tracking/config/samples_converted.csv`.

* `config/mikk_vcf_cram2line-id_dupes-ed.txt` copied from https://github.com/brettellebi/mikk_genome/blob/master/data/20200305_cram2line_full_dupes-edited.txt, then edited manually to replace the underscores in the line IDs with hyphens.

