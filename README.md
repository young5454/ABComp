# ABComp : Assembly Polishing and Bacterial Whole-genome Comparison Pipeline for Multi-group Clinical Isolates

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)

<img src="docs/figures/ABComp.v2.png" width="1000px" align="center" />

# Getting Started
## Installation

Clone and enter repository to create your workspace directory :

```bash
git clone https://github.com/young5454/ABComp.git
cd ABComp && mkdir {your_workspace}
```

⚠️ Make sure your workspace name is correctly specified in the configuration file (please see section [Configuration file](#configuration-file) for more detail).

## Directory structure and file names
ABComp starts with two main type of input data, the draft assembly of clinical isolates and Illumina paired-end short reads. In order for ABComp to run smoothly, the input data should be saved in a proper format of both directory structure and file names. 

All of the genomes, reads, and reference genbank file should be saved inside the workspace directory as below :

```
0.Assembly
├── {group}_{strain}
│   ├── genome
│   │   └── {group}_{strain}.fasta
│   └── reads
│       ├── {group}_{strain}_1.fastq.gz
│       └── {group}_{strain}_2.fastq.gz
└── {ref}
    ├── genome
    │	└── {ref}.fasta
    └── {ref}.gb
```

⚠️ Make sure all the file extensions strictly follow the names (`.fasta`, `.fastq.gz`, `.gb`) provided above. Also, the fastq read files should be gzipped beforehand. The fixed directory names should **NOT** be altered.

⚠️ Make sure your `{group}` and `{strain}` names do **NOT** contain a underscore (_) inside. The name formats are intended to have the unique underscore so that during the Snakemake run, the group and strain wildcard info are automatically extracted and saved. Also, the file names should perfectly match the names provided in the configuration file (please see section [Configuration file](#configuration-file) for more detail).

The assembly polishing module will create a directory `trimmed/` inside each `{group}_{strain}/` folder after running adapter trimming, and save the trimmed reads. Also, the polished assemblies will be saved inside the `genome/` directory with a name format of `{group}_{strain}_polished.fasta`. 

⚠️ If you have your assemblies ready and plan to run only the downstream comparative analysis module, please save all of your genome assemblies with the name format of `{group}_{strain}_polished.fasta` inside the `genome/` directory.

## Configuration file
ABComp requires two configuration files for running the pipeline. These yaml files can be found in the `config/` directory. 

`config.yml` is a default configuration setting for the overall Snakemake run. Make sure you specify the correct parameters and directory names of your preference.

`groups_original.yml` is a configuration file for the complete group-strain information of your clinical isolates. Below is an example yaml file of a 2-group, 6-strain setting :

```YAML
IPMR:
    - BFO17
    - BFO18
    - BFS01
IPMS:
    - BFO42
    - BFO67
    - BFO85
```

⚠️ Make sure that the group and strain names perfectly match with the draft assembly name format : `{group}_{strain}.fasta`. For example, the BFO17 strain should have the directory and file name of `IPMR_BFO17`.

## Conda environments
ABComp suggests single conda environment for a single software for avoiding potential dependency errors. Before running the pipeline, make sure you have all of the environments set with the provided yaml files. The yaml file of each individual environment can be found in `workflow/envs/`.

ABComp provides a single bash script `create_envs.sh` that creates all of the required conda environments at once. You can run the script with the following command :

```
bash create_envs.sh 
```

## Pre-downloaded databases and files
ABComp requires pre-downloaded databases for running **BUSCO** assembly assessment and **EggNOG-mapper**. Please refer to the website or Github of each software for more information :

### BUSCO lineage datasets
[BUSCO](https://busco.ezlab.org/busco_userguide.html#lineage-datasets/)

BUSCO lineage datasets contain sets of highly conserved single-copy genes specific to different evolutionary lineages. These datasets are used to evaluate the completeness of genome assemblies by comparing the presence of these key genes in the assembly. Pre-download the appropriate lineage dataset according to the lineage of the input strain. 

⚠️ Make sure the lineage datasets are saved inside a directory named `busco_downloads/`, inside your workspace directory.

### EggNOG annotation databases
[EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)

EggNOG-mapper databases consist of precomputed orthologous groups and functional annotations derived from a vast range of organisms. These databases enable the functional annotation of genes in your genome by mapping them to known orthologs and their associated functional data. Pre-downloading the necessary EggNOG-mapper databases ensures that the pipeline can efficiently perform functional annotation. 

⚠️ Please carefully read through the **Setup** section in the [EggNOG-mapper wiki](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#user-content-Setup). ABComp assumes you have downloaded all the annotation databases with running the script provided by EggNOG-mapper with a **default option**:

```bash
download_eggnog_data.py
```

That is, a directory named `data/` must be created within your EggNOG-mapper directory, then the above script is run to download all databases provided by EggNOG-mapper.

### Trimmomatic adapter FASTAs
For running trimmomatic, please save the required adapter sequence FASTAs you need to specify in the `adapters/` directory inside your workspace. The sequence FASTA files can be acquired from the [Trimmomatic github](https://github.com/usadellab/Trimmomatic).

# Running the Pipeline
Once you have followed the Getting Started section and setup everything, you can start the Snakemake run!

## Assembly polishing
The below code will start the assembly polishing pipeline with using up to 4 cores. Make sure you activate the default conda environment before your Snakemake run.

```bash
conda activate env_default

snakemake \
	--snakefile workflow/Snakefile_assembly_polishing.py \
	--use-conda \
	--cores 4
```

## Comparative genomics
The below code will start the comparative analysis pipeline with using up to 4 cores. This task contains gene annotation, pangenome analysis, antibiotic resistance gene finding and more. Make sure you activate the default conda environment before your Snakemake run.

```bash
conda activate env_default

snakemake \
	--snakefile workflow/Snakefile_downstream.py \
	--use-conda \
	--cores 4
```

# Pathogenic marker discovery
<img src="docs/figures/downstream.png" width="1000px" align="center" />
The pathogenic marker discovery module operates on the result of commandline BLASTp, therefore the user should first manually conduct the alignment process. In this example, let group B’s core set is queried to the total annotation set of a single strain X.

First, create the subject database :

```bash
# Creating subject DB
makeblastdb \
	-in strain_X_annot.fasta \
	-dbtype prot \
	-out strain_X_annnot
```

Next, run BLASTp to query all the core entries to the subject database :

```bash
# Advanced: Coverage option
blastp \
	-query group_B_core.fasta \
	-db strain_X_annnot \
	-out core_B_to_strain_X.csv \
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
```

⚠️ Make sure to use the output format specified above, since query coverage is a fixed metric used for parsing the BLASTp results.
Finally, running `pairwise_analysis_plus.py` will do the following :

1. Parse raw BLASTp output table with percent identity, E-value and query coverage threshold for determining presence of query entry within subject set
2. For the passed queries, the single best query-to-subject match (Best Matches) will be saved as a BLASTp output table
3. The best match queries and subjects will be saved as a FASTA file
4. Not-found queries will be saved as a BLASTp output table with all results
5. Not-found queries will be saved as a FASTA file

Below is a run with setting the percent identity threshold as 35, E-value threshold as 0.001, and query coverage threshold as 80 :

```bash
python pairwise_analysis_plus.py \
	-b core_B_to_strain_X.csv \         # Path to raw BLASTp output table
	-qf core_B.fasta \                  # Path to query set FASTA
	-sf strain_X.fasta \                # Path to subject set FASTA
	-p 35 \                             # Percent identity threshold
	-e 0.001 \                          # E-value threshold
	-c 80 \                             # Query coverage threshold
	-s Pairwise.Plus/ \                 # Path to save all final data
	-bmt best_match.csv \               # File name of best match BLASTp output table
	-nft not_founds.csv \               # File name of not-founds BLASTp output table
	-pqf pass_query.fasta \             # File name of best match query FASTA
	-psf pass_subject.fasta \           # File name of best match subject FASTA
	-nff not_founds.fasta               # File name of not-founds query FASTA
```