# MixInfect2

### Detecting multi-strain microbial samples from short-read whole genome sequencing data

This tool is an update of a previously published method for detecing mixed infection in TB samples (github.com/bensobkowiak/Mixinfect) - [link to paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4988-z).

A new preprint detailing the performance of MixInfect2 is available [here](https://www.biorxiv.org/content/10.1101/2024.04.26.591283v1).

This page contains two main scripts that are executable with Rscript in the command line using VCF files from a variant calling software (e.g., GATK, BCFtools etc.), along with test data to run the tools.

- **MixInfect2.R** estimates the number of strains present in each sample by using a Gaussian Mixture Model to cluster read frequencies at mixed loci identified per sample from the VCF file.

- **reconstructConstitutents.R**  uses the same VCF file along with the output file from MixInfect2.R to reconstuct the constituent sequences of samples that were found to contain more than one possible strain. Constituent sequences in mixed samples are predicted using two methods:
    - Assigning strain sequences using consensus read frequencies
    - Matching to the genetically closest non-mixed sample in the VCF file.

Please see the pre-print above for recommendations on which method provides the most accurately representation of a consensus sequence in mixes. 

_Note. this tool has only been tested on mixes of two strains and from short-read Illumina data but we are in the process of updating it to incorporate more than two strain mixes and assemblies of long-read sequence data._

<br>

## Running MixInfect2 and reconstructing constituent sequences of mixed samples

### MixInfect2.R

A VCF file is the main input to run MixInfect2 and is specified with the ```--VCFfile``` option.

A user-defined prefix for all output files is specified with the ```--prefix``` option.

A CSV file containing genes to mask in the VCF file that are not considered when estimating mixed samples e.g., AMR-conferring genes or known repeats, is is specified with the ```--maskFile``` option. (recommended) 

#### Usage

MixInfect2 can then be run using Rscript in the command line as follows:

```bash
Rscript MixInfect2.R --VCFfile input.vcf --prefix output --maskFile MaskedRegions.csv 
```

The following options can be specified:

| Option                 | Type       | Default | Description                                                                                                                             |
|------------------------|------------|---------|-----------------------------------------------------------------------------------------------------------------------------------------|
| `--VCFfile`            | character  |     | VCF file - must have GT and AD fields                                                                                                   |
| `--prefix`             | character  | output  | Output file prefix                                                                                                                      |
| `--maskFile`           | character  | NULL    | CSV file with regions to mask, with start position in column 1 and end position in column 2                                             |
| `--useFilter`          | logical    | TRUE    | Use the 'FILTER' column in VCF file to filter SNPs                                                                                      |
| `--minQual`            | numeric    | 20      | Minimum per loci quality                                                                                                                |
| `--LowCov`             | numeric    | 10      | Minimum read depth at site to call either a cSNP or hSNP allele frequency                                                               |
| `--minDepth`             | integer    | 5      | Minimum read depth of minor frequency allele for a mixed call                                                                                             |
| `--popFreq_threshold`  | numeric    | 1       | Remove hSNPs found in greater than this proportion of sequences in VCF (set as 1 for single sample VCF)                                                                                 |
| `--SNPwindow`          | numeric    | 100     | Take the median of hSNP allele frequencies within this distance on the genome                                                           |
| `--n_threads`          | numeric    | 4       | Number of threads to use                                                                                                                |


#### Output

MixInfect2.R produces two output files

- _prefix_MixSampleSummary.csv_ - A CSV file showing whether each sample is mixed or non-mixed, the number of hSNPs and cSNPs, proportion of hSNPs/cSNPs, the number of estimated strains, and, for mixed samples, the estimated major strain proportion.
- _prefix_BICvalues.csv_ - For each mixed sample, the corresponding BIC values for the hSNP read frequencies at cluster of size 2, 4 and 6. 

<br>

### reconstructConstituents.R

The same VCF file as used to run MixInfect2 is specified with the ```--VCFfile``` option.

A user-defined prefix for all output files is specified with the ```--outputprefix``` option.

The main output CSV file from the MixInfect2.R script that classifies samples as mixed or non-mixed (with the suffix "_MixSampleSummary.csv") is specified with the ```--MixInfect2Result``` option.

A CSV file containing genes to mask in the VCF file that are not considered when estimating mixed samples e.g., AMR-conferring genes or known repeats, is is specified with the ```--maskFile``` option. (recommended) 

#### Usage

reconstructConstituents.R can then be run using Rscript in the command line as follows:

```bash
Rscript reconstructConstituents.R --VCFfile input.vcf --prefix output --MixInfect2Result output_MixSampleSummary.csv --maskFile MaskedRegions.csv 
```

The following options can be specified:

| Option                       | Type       | Default | Description                                                                                                                           |
|------------------------------|------------|---------|---------------------------------------------------------------------------------------------------------------------------------------|
| `-v, --VCFfile`              | character  |         | Input VCF file (same as for MixInfect2)                                                                                                                       |
| `-o, --prefix`         | character  | output  | Prefix for output files                                                                                                               |
| `-r, --MixInfect2Result`     | character  |         | output CSV file generated using MixInfect2                                                                                                   |
| `-f, --maskFile`             | character  | NULL    |CSV file with regions to mask, with start position in column 1 and end position in column 2                |
| `-q, --minQual`              | integer    | 20      | Minimum per loci quality                                                                                              |
| `-c, --closestStrain`        | logical    | TRUE    | Reconstruct constituent strains using closest strain method (set to FALSE for single sample VCFs)                                                                          |
| `-x, --maxDistance`          | integer    | 5000    | Maximum distance to closest non-mixed strain (If -c is TRUE)                                                                                         |
| `-p, --popFreq_threshold`    | numeric    | 1       | Remove hSNPs found in greater than this proportion of sequences in VCF (set as 1 for single sample VCF)                                                              |
| `--LowCov`             | numeric    | 10      | Minimum read depth at site to call either a cSNP or hSNP allele frequency                                                               |
| `-d, --minDepth`             | integer    | 5      | Minimum read depth of minor frequency allele for a mixed call                                                                                             |
| `-m, --mixProp`              | numeric    | 0.9     | Proportion allele frequency at hSNPs to assign call in non-mixed strains                                                              |


#### Output

reconstructConstituents.R produces two output files and one optional output file:

- _prefix_constituents.fasta_ - A FASTA file containing a concatenated SNP alignment for all non-mixed strains and estimated reconstructed sequences of the mixed samples.
- _prefix_SNPindex.txt_ - A SNP index of the sites retained in the FASTA file, corresponding to their position in the reference strain as determined in the input VCF file. 

(optional if closestStrain = TRUE):

- _prefix_closestStrainSummary.csv_ - A CSV file containing, for each mixed sample, the names and genetic distance to the closest non-mixed strain in the VCF file for both the major and minor constituent strains.

## Citation

Please cite the following paper if you use MixInfect2 for your analysis, [Sobkowiak et. al. 2024](https://www.biorxiv.org/content/10.1101/2024.04.26.591283v1).
