# Notes

## 20 August 2022

SNP hunting:

* For state 1, the SNP at ENSORLG00000009385 is *VAV1*: vav guanine nucleotide exchange factor 1. Cancer phenotypes in human, immune system phenotypes in mouse.

* ENSORLG00000028050 is a novel gene.

* ENSORLG00000024879 is a novel gene. Phenotypes: bone mineral density in human, and PR interval (cardiovascular) in mouse.

* ENSORLG00000013962 is echo-ADP-ribosyltransferase 5-like. No phenotypes.

* ENSORLG00000029574 is a novel gene. No phenotypes.

* 9:9802754: ENSORLG00000024663: grid2 (Projected ZFIN): glutamate ionotropic receptor delta type subunit 2. Phenotypes include ataxia in humans, and neurological defects in mouse.

* 10:15319434: ENSORLG00000024866: protocadherin alpha-C2-like, axon branching receptor (famous).

* 10:18537719: ENSORLG00000006464: nlgn3a (neuroligin 3a). synapse related protein involved in Autism. Phenotypes include autism in human, and extensive neurological phenotypes in mouse. 

## 10 August 2022

Something very strange:

We have duplicated rows for the samples collected on 20191113 in the right tank, starting from 20191113_1247

SOMETHING HAPPENED DURING THE MOVEMENT METRICS STEP

Due to the low frame rate in those videos.

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

