# MixInfect2

### Calling mixed microbial samples from variant calling on whole genome sequencing data

This tool is an update of a previously published method for detecing mixed infection in TB samples (github.com/bensobkowiak/Mixinfect) - [link to paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4988-z).

A new preprint detailing the performance of this tool is available [here](https://www.biorxiv.org/content/10.1101/2024.04.26.591283v1).

This tool can be run using the MixInfect2.R script. The "MaskedRegions.csv" file can be used with the 'maskFile' parameter to mask regions identified as pe/ppe genes, known antimicrobial resistance-conferring genes, and regions of excess repeats.

The MixInfect2_PL.R script can be run on command line using Rscript and allows for multi-threading
