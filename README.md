# MixInfect2

### Detecting multi-strain microbial samples from short-read whole genome sequencing data

This tool is an update of a previously published method for detecing mixed infection in TB samples (github.com/bensobkowiak/Mixinfect) - [link to paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4988-z).

A new preprint detailing the performance of MixInfect2 is available [here](https://www.biorxiv.org/content/10.1101/2024.04.26.591283v1).

This page contains two main scripts that are executable with Rscript in the command line using VCF files from a variant calling software (e.g., GATK, BCFtools etc.), along with test data to run the tools.

First MixInfect2.R estimates the number of strains present in each sample by using a Gaussian Mixture Model to cluster read frequencies at mixed loci identified per sample from the VCF file. 

The next script (reconstructConstitutents.R) then uses the same VCF file along with the output file from MixInfect2.R to reconstuct the constituent sequences of samples that were found to contain more than one possible strain. Constituent sequences in mixed samples are predicted using two methods: 1. assigning strain sequences using consensus read frequencies and 2. matching to the genetically closest non-mixed sample in the VCF file. Please see the pre-print above for recommendations on which method provides the most accurately representation of a consensus sequence in mixes. 

_Note. this tool has only been tested on mixes of two strains and from short-read Illumina data but we are in the process of updating it to incorporate more than two strain mixes and assemblies of long-read sequence data._

### Running MixInfect2.R and reconstructConstituents.R

Both scripts are best run in the command line with Rscript and can be run with default parameters as follows:

```bash
Rscript MixInfect2.R --VCFfile input.vcf --prefix output --maskFile MaskedRegions.csv 
```

```bash
Rscript reconstructConstituents.R --VCFfile input.vcf --outputprefix output --MixInfect2Result output_MixSampleSummary.csv --maskFile MaskedRegions.csv 
```

