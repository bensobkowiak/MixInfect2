# MixInfect2

### Calling mixed microbial samples from variant calling on whole genome sequencing data

This tool is an update of a previously published method for detecing mixed infection in TB samples (github.com/bensobkowiak/Mixinfect) - [link to paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4988-z).

A new preprint detailing the performance of this tool is available [here](https://www.biorxiv.org/content/10.1101/2024.04.26.591283v1).

This page contains two main scripts to first detect mixed infections from (MixInfect.R) and then reconstruct the constituent strains (reconstructConstitutents.R) using two methods - using consensus allele frequencies and matching to the genetically closest non-mixed strain in the dataset. Please see the pre-print above for recommendations on which method provides the most accurately representation of a consensus sequence in mixes.

Both scripts are best run in the command line with Rscript as follows:

```bash
Rscript MixInfect2.R --VCFfile input.vcf --prefix output --maskFile MaskedRegions.csv ...
```

```bash
Rscript reconstructConstituents.R --VCFfile input.vcf --outputprefix output --MixInfect2Result output_MixSampleSummary.csv --maskFile MaskedRegions.csv ...
```

