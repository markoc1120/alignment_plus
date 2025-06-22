# Alignment Plus

A simple command-line tool for performing pairwise and multiple sequence alignment (MSA). Pairwise alignment supports both linear and affine gap costs and it gives us all of the optimal alignments, while the MSA only supports linear gap costs and only gives back a single optimal alignment.

## Features

- **Pairwise Alignment**: Global alignment between two sequences
- **Multiple Sequence Alignment (MSA)**: Alignment of more than two sequences
- **Customizable Scoring**: User-defined scoring matrices and gap penalties
- **Flexible Output**: Results can be saved to file or displayed in terminal

## Installation

### Using Conda (Recommended)

```bash
# Create environment from the provided environment.yaml
conda env create -f environment.yaml

# Activate the environment
conda activate alignment_plus
```

### Manual Installation

Required dependencies:

- Python 3.13+
- NumPy 2.2.2+
- Pandas 2.2.3+
- BioPython 1.85+

```bash
pip install numpy pandas biopython
```

## Usage

### Basic Command Structure

```bash
python main.py -i <input_file> -at <alignment_type> -s <settings> [OPTIONS]
```

### Required Arguments

- `-i, --input`: Name of the input FASTA file containing sequences to align
- `-at, --alignment-type`: Type of alignment (`pair` for pairwise, `multi` for MSA)
- `-s, --settings`: Name of the settings JSON file for custom parameters

### Optional Arguments

- `-o, --output`: Output path to save alignment results (FASTA format)
- `-gm, --gap-model`: Gap penalty model (`linear` or `affine`, default: `linear`)

### Examples

#### Pairwise Alignment with Linear Gap Model

```bash
python main.py -i pairwise_input.fasta -at pair -s settings.json -o output/pairwise_result.fasta
```

#### Pairwise Alignment with Affine Gap Model

```bash
python main.py -i pairwise_input.fasta -at pair -s settings.json -gm affine -o output/pairwise_affine_result.fasta
```

#### Multiple Sequence Alignment with Linear Gap Model

```bash
python main.py -i msa_input.fasta -at multi -s settings.json -o output/msa_result.fasta
```

## Input Format

### FASTA Files

Input sequences must be in FASTA format:

```fasta
>sequence1
ATGGATTTATCTGCTCTTCG
>sequence2
GTTCCGAAAGGCTAGCGCTAGGCGCC
```

#### Requirements:

- **Pairwise alignment**: Exactly 2 sequences
- **MSA**: More than 2 sequences
- Supported nucleotides: A, C, G, T (case-insensitive)

### Settings File

The settings file should be in JSON format and can include:

```json
{
  "gapopen": 5,
  "gapextend": 5,
  "score_matrix_file": "score_matrix.csv"
}
```

#### Parameters:

- `gapopen`: Gap opening penalty (required)
- `gapextend`: Gap extension penalty (required for affine model)
- `score_matrix_file`: Name of the scoring matrix CSV file

### Scoring Matrix

The scoring matrix should be a CSV file with nucleotide scores:

```csv
A, C, G, T
0, 5, 2, 5
5, 0, 5, 2
2, 5, 0, 5
5, 2, 5, 0
```

## Algorithms

### Pairwise Alignment

- **Linear**: Dynamic programming with linear gap penalties
- **Affine**: Three-matrix approach (S, D, I matrices) for affine gap penalties

### Multiple Sequence Alignment

- **Linear**: Multi-dimensional dynamic programming
- Currently supports linear gap model only
