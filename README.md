# Repeat-Aware_phylogenetic_distance_estimator

This tool efficiently estimates phylogenetic distances and constructs a Neighbor-Joining (NJ) tree among multiple reference or assembly genomes in parallel. It is based on the [Repeat-Aware Substitution Rate Estimator](https://github.com/medvedevgroup/Repeat-Aware_Substitution_Rate_Estimator) and [Repeat-Aware Mutation Rate Estimator](https://github.com/Wu-Haonan/Repeat-Aware_mutation_rate_estimator), leveraging the power of **KMC** for k-mer counting and **sourmash** for sketching.

## Installation

### Git Clone with Submodules

```bash 
- Method 1: Clone with submodules in one step (recommended) 
git clone --recurse-submodules git@github.com:Wu-Haonan/Repeat-Aware_phylogenetic_distance_estimator.git
- Method 2: Clone first, then initialize submodules 
git clone git@github.com:Wu-Haonan/Repeat-Aware_phylogenetic_distance_estimator.git
cd Repeat-Aware_phylogenetic_distance_estimator 
git submodule update --init --recursive
```

### Install required packages using Conda

We recommend using Conda to set up the environment:


Install all the packages based on `requirements.txt`

```bash
conda create -n <env_name> python=3.10 --file requirements.txt -c conda-forge -c bioconda -y
conda activate <env_name>
```

### Compile `kmc_histogram.cpp`

```bash
cd ./cpp/
make
```


## Quick Start

We provide a small example under `example_data/`:

```bash
python Phylogenetic_distance_estimator.py \
    --config example_data/config.csv \
    --k 31 \
    --theta 0.01 \
    --cores 4 \
    --outdir results
```

## Usage

### Input Configuration

Create a CSV configuration file with columns: `label,filepath`

Example `config.csv`:

```csv
label,filepath
Genome_A,/path/to/genomeA.fasta
Genome_B,/path/to/genomeB.fasta
Genome_C,/path/to/genomeC.fasta
```

### Basic Usage

```bash
python Phylogenetic_distance_estimator.py \
    --config config.csv \
    --k 31 \
    --theta 0.01 \
    --cores 4 \
    --outdir results
```

### Distance Computation Modes

- **`--mode half`** (default): Compute unique pairwise distances (like, A-B, A-C, B-C)
- **`--mode all`**: Compute all directional distances (like, A-B, B-A, A-C, C-A, B-C, C-B)

Example for all-vs-all comparison:

```bash
python Phylogenetic_distance_estimator.py \
    --config config.csv \
    --mode all \
    --k 31 \
    --theta 0.01 \
    --cores 8 \
    --outdir results
```

## Parameters

### Required Parameters

- **`--config`**: CSV file with `label,filepath` columns specifying genome names and file paths

### Optional Parameters

- **`--mode`**: Distance computation mode (`half` or `all`, default: `half`)
- **`--k`**: K-mer size (default: `31`). Must be between 1-256 due to KMC limitations
- **`--theta`**: FracMinHash sketching rate (default: `0.01`). Controls sampling rate for acceleration
- **`--cores`**: Number of CPU cores for parallel processing (default: `4`)
- **`--kmc-threads`**: Threads per KMC task (default: `1`)
- **`--outdir`**: Output directory (default: `results`)
- **`--cleanup`**: Remove intermediate files `./middle_results` after completion

## Output Files

The tool generates several output files in the specified output directory:

1. **`pairwise_distances.csv`**: Pairwise distance matrix in CSV format

2. **`nj_tree.nwk`**: Neighbor-joining phylogenetic tree in Newick format

3. Intermediate files (if  `--cleanup` not specified):

    - KMC databases: `kmc_{label}_k{k}/`
   - Sourmash signatures: `{label}_k{k}_sig.sig`
   - K-mer histograms: `{label}_k{k}_hist.csv`
   - Individual distance files: `{labelA}___{labelB}.csv`

## FracMinHash Sketching

The `--theta` parameter controls FracMinHash sketching to approximate intersection sizes:

**Note**: Sourmash requires a minimum scaled parameter of 100, corresponding to `--theta` $\leq 0.01$. Values above 0.01 will be automatically adjusted by sourmash.

## Performance Considerations

### Parallelization

The tool supports multi-level parallelization:

- **Process-level**: Multiple genomes processed simultaneously
- **Thread-level**: KMC uses multiple threads per task

### Memory Usage

- KMC temporary files require significant disk space
- Memory usage scales with genome size and k-mer diversity
- Consider using `--cleanup` for large datasets to save disk space

## Algorithm Overview

1. **K-mer Counting**: Uses KMC to count k-mers in each genome
2. **Sketching**: Applies FracMinHash via sourmash for efficient intersection estimation
3. **Histogram Generation**: Creates k-mer frequency distributions
4. **Distance Calculation**: Estimates mutation rates using repeat-aware methods
5. **Tree Construction**: Builds neighbor-joining phylogenetic tree using BioPython

