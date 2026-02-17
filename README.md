
[![DOI](https://zenodo.org/badge/1072837401.svg)](https://doi.org/10.5281/zenodo.18199468)

# MauE_supplementary_code
This is a repository for the code used in dataset construction for the scientific publication center in MauE as a disulfide bond formation protein.

All code used is in this github. 
supplementary data files can be gound at https://doi.org/10.5281/zenodo.18197764

the notebooks listed are sequential and are set up so they can be repeated by others to reproduce the data analysis in the listed order
additionally taxanomic metadata was annotated using an NCBI taxa dump from 2025, which can be recreated in the notebooks or found in the supplementary
data at https://doi.org/10.5281/zenodo.18197764 under the phylogeny tables as a "pkl" file. This work uses the ENTREZ tool in Biopython for anything not found in the
NCBI 2025 taxa dump.

This data was constructed before in 10.2024 and there has since been updates to Uniprot and Interpro that will change output if starting from 
raw API pulls. To recreate our exact materials use the initial data tables in the directory "uncombined_rawdata" as a starting point for notebook
1, and continue from there.

