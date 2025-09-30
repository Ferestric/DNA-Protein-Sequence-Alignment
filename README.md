# DNA-Protein Sequence Alignment

Align two DNA or protein sequences using Global, Semi-Local, or Local sequence alignment methods.

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
  - [Command Line / Script](#command-line--script)
  - [Examples](#examples)
- [Input & Output](#input--output)
- [Scoring Matrices](#scoring-matrices)
- [Project Structure](#project-structure)
- [Testing](#testing)
- [Contributing](#contributing)
- [License](#license)

## Introduction
This project implements classical sequence alignment algorithms (Global, Semi-Local, Local) to align DNA or protein sequences. You can compute optimal alignments, score them, and examine the alignment path. It supports both nucleic acid and amino acid sequences.

## Features
- Global (Needleman–Wunsch) alignment  
- Local (Smith–Waterman) alignment  
- Semi-local alignment (hybrid)  
- Support for DNA (nucleotides) and protein sequences  
- Customizable scoring matrices / substitution tables  
- Outputs the aligned sequences and scores  
- Example test files and output included  

## Prerequisites
- Python 3.x (tested on 3.7+)  
- No external libraries (pure Python implementation)  

## Installation
Clone the repository:
```bash
git clone https://github.com/Ferestric/DNA-Protein-Sequence-Alignment.git
cd DNA-Protein-Sequence-Alignment
```

You may optionally set up a virtual environment:
```bash
python3 -m venv venv
source venv/bin/activate    # macOS / Linux
venv\Scripts\activate       # Windows
```

## Usage

### Command Line / Script
Run the main script:
```bash
python align.py <seq1> <seq2> --mode [global|local|semilocal] [--matrix MATRIX_FILE] [--gap-open GO] [--gap-extend GE]
```

**Options**
- `--mode` — alignment mode (global, local, semilocal)  
- `--matrix` — (optional) substitution matrix file (for proteins)  
- `--gap-open` — gap opening penalty  
- `--gap-extend` — gap extension penalty  

### Examples
Example run with a protein substitution matrix:
```bash
python align.py seqA.fasta seqB.fasta --mode local --matrix BLOSUM45
```

See `out.txt` in the repository for a sample alignment result.

## Input & Output
- **Input**: sequences (as strings or FASTA files) and optional substitution matrix  
- **Output**: aligned sequences (with gaps), alignment score, and alignment path  
- Example results are in `out.txt`  

## Scoring Matrices
- Protein alignments: substitution matrices such as **BLOSUM45** are supported  
- DNA alignments: simple match/mismatch scheme defined internally  
- Additional matrices (e.g., BLOSUM62, PAM250) can be added to the repo and passed with `--matrix`  

## Project Structure
```
.
├── align.py             # Main alignment script
├── BLOSUM45             # Example substitution matrix file
├── dnaMatrix            # DNA scoring definitions
├── dnaMatrix2           # Alternative DNA matrix file
├── out.txt              # Example output
├── .gitignore
├── LICENSE
└── README.md
```

## Testing
Test with included files:
```bash
python align.py test_seq1.fasta test_seq2.fasta --mode global
```

Check the output against expected alignments.

## Contributing
Contributions are welcome! You can:
- Add support for more substitution matrices  
- Optimize performance  
- Add a GUI or web interface  
- Expand test cases  

Fork the repo, commit your changes, and open a pull request.

## License
This project is licensed under the **MIT License**.
