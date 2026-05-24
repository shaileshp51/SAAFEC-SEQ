# SAAFEC-SEQ

SAAFEC-SEQ predicts the effect of missense mutations on protein thermodynamic stability using a sequence-based workflow that combines local sequence context, mutation descriptors, and PSSM-derived features.

## Supported platforms

SAAFEC-SEQ is intended for Linux and Linux-like systems, including macOS.

The package expects:

- a local PSI-BLAST installation
- a locally indexed UniRef100 database
- a writable user configuration file at `~/.saafec-seqrc`

Windows is not a supported target.

## Command name

After installation, the public command is:

```bash
saafec-seq
```

## Requirements

### Tested runtime stack

The following versions were validated in the worker environment used for deployment:

- Python 3.11
- `numpy==1.26.4`
- `scipy==1.14.1`
- `pandas==2.2.3`
- `scikit-learn==1.6.0`
- `xgboost==0.82`
- `requests==2.32.3`


## External tools and data: BLAST+ and UniRef100 setup

- BLAST+ `psiblast`
- a local UniRef100 PSI-BLAST database

SAAFEC-SEQ requires a local BLAST+ installation and a local UniRef100 database. Install BLAST+ from the official [NCBI BLAST+ download page](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html), where NCBI provides the command-line tools and source, and obtain UniRef100 from UniProt’s [downloads page](https://www.uniprot.org/help/downloads) or [UniRef page](https://www.uniprot.org/uniref). :contentReference[oaicite:0]{index=0}

UniRef100 is the 100% clustered UniProt reference set used to reduce redundancy while retaining sequence coverage. :contentReference[oaicite:1]{index=1}

If you already have a formatted BLAST database, point SAAFEC-SEQ at the database prefix. Run `makeblastdb` only when starting from a FASTA file and creating the local BLAST database yourself. :contentReference[oaicite:1]{index=1}

For example, if you are creating a local UniRef100 protein database from a FASTA file, a typical command is:
```bash
makeblastdb -in uniref100.fasta -dbtype prot -out /path/to/uniref100/uniref100
```

Now, update your local SAAFEC-SEQ configuration so it points to the `psiblast` binary directory and the UniRef100 database directory.

If you deviate from the tested stack, keep the same major workflow but expect to verify compatibility in your own environment.

## Installation

### Recommended: Miniforge or Miniconda

The project is tested against a conda-based Python stack. That keeps the Python version and compiled dependencies aligned with the validated environment.

Create and activate a dedicated environment:

```bash
conda create -n saafec-seq python=3.11 pip
conda activate saafec-seq
```

Install the pinned runtime packages from `conda-forge`:

```bash
conda install -y -c conda-forge \
  numpy==1.26.4 \
  scipy==1.14.1 \
  pandas==2.2.3 \
  scikit-learn==1.6.0
```

Install the remaining Python packages with pip:

```bash
python -m pip install --no-cache-dir xgboost==0.82 --no-deps
python -m pip install --no-cache-dir requests==2.32.3
```

### Install the package

After the runtime stack above is available, install SAAFEC-SEQ from source:

```bash
pip install .
```

For development:

```bash
pip install -e .
```

## Configuration

Configuration is stored in a plain text file:

```text
~/.saafec-seqrc
```

The file uses simple `KEY=VALUE` lines. Example:

```text
PSIBLAST_BIN_DIR=/opt/shared/blast
UNIREF_DB_DIR=/opt/shared/uniref100
UNIREF_DB_NAME=uniref100
BLAST_NUM_THREADS=4
```

### Configuration order

On startup, SAAFEC-SEQ uses this order:

1. `~/.saafec-seqrc`
2. prompt the user for values if the file is missing or invalid
3. save the working values back to `~/.saafec-seqrc`

If the saved configuration points to a missing `psiblast` executable or an invalid database location, the program prompts again and rewrites the file.

## First run behavior

On the first run, or whenever the saved config is invalid, SAAFEC-SEQ asks for:

- path to the PSI-BLAST binary directory
- path to the UniRef100 database directory
- UniRef100 database name
- number of BLAST threads

The answers are written to `~/.saafec-seqrc` for later use.

## Usage

### Single mutation

```bash
saafec-seq \
  -i example.fasta \
  -p 42 \
  -w A \
  -m V \
  -o output.out
```

### Multiple mutations

```bash
saafec-seq \
  -i example.fasta \
  -f mutations.txt \
  -o output.out
```

### Arguments

- `-i`, `--input-sequence`: FASTA file containing the protein sequence
- `-o`, `--output`: output file name, default: `output.out`
- `-v`, `--verbose`: print verbose messages

#### Single-mutation mode

- `-p`, `--position`: 1-based mutation position
- `-w`, `--wild-type`: wild-type amino acid one-letter code
- `-m`, `--mutant`: mutant amino acid one-letter code

#### Multiple-mutation mode

- `-f`, `--mutation-list-file`: plain text file with one mutation per line

## Input format

### FASTA sequence

The input sequence must contain only standard amino acids:

`A C D E F G H I K L M N P Q R S T V W Y`

Header lines beginning with `>` are allowed.

### Multiple mutation list

One mutation per line:

```text
A 42 V
G 101 D
L 155 P
```

Blank lines and comment lines beginning with `#` are ignored.

## Output format

Example output:

```text
Position Wild Mutant ddG Type
42 A V 1.23 Destabilizing
```

Possible result types:

- `Neutral`
- `Stabilizing`
- `Destabilizing`

## Sequence database setup

PSI-BLAST requires a properly formatted local database.

Example:

```bash
makeblastdb -in uniref100.fasta -dbtype prot -out /data/uniref100/uniref100
```

Then place the matching values into `~/.saafec-seqrc`:

```text
PSIBLAST_BIN_DIR=/opt/shared/blast
UNIREF_DB_DIR=/data/uniref100
UNIREF_DB_NAME=uniref100
BLAST_NUM_THREADS=6
```

## Troubleshooting

### `saafec-seq: command not found`

Check that the package is installed in the active conda environment and that the environment is activated.

### `psiblast` not found

Check the configured BLAST directory and ensure the executable exists.

### UniRef100 database not found

Check the configured database directory and prefix.

### Invalid FASTA input

The input sequence must contain only standard amino acids.

### Mutation mismatch

The wild-type residue must match the residue at the chosen position in the input sequence.

### Sequence too short

Very short sequences may not produce reliable features or may fail validation.

## For unusual environments

If your system has a non-standard software layout, edit `~/.saafec-seqrc` after first run or let the program prompt you again by removing that file.

This project is intended for Linux and Linux-like systems, including macOS. Windows is not a supported target.

## Citation

If you use SAAFEC-SEQ in published work, please cite:

Li, G.; Panday, S.K.; Alexov, E. *SAAFEC-SEQ: A Sequence-Based Method for Predicting the Effect of Single Point Mutations on Protein Thermodynamic Stability*. Int. J. Mol. Sci. 2021, 22, 606.

## License

SAAFEC-SEQ is licensed under the GNU Affero General Public License v3.0 or later (AGPL-3.0-or-later).
