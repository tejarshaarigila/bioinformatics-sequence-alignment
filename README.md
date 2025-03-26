# Sequence Alignment Algorithms

This repository contains Python implementations of two classical sequence alignment algorithms: **Needleman-Wunsch** (Global Alignment) and **Smith-Waterman** (Local Alignment). These algorithms are essential in bioinformatics for comparing and aligning biological sequences like DNA, RNA, and proteins.

## Features

- **Needleman-Wunsch Algorithm**: Performs global sequence alignment using dynamic programming.
- **Smith-Waterman Algorithm**: Performs local sequence alignment using dynamic programming with affine gap penalties.

## Algorithms Overview

- **Needleman-Wunsch (Global Alignment)**:
  - Aligns the entire sequences from start to end.
  - Computes an optimal alignment for the entire sequence length, considering all possible pairings.
  
- **Smith-Waterman (Local Alignment)**:
  - Aligns the most similar subsequences of two sequences.
  - Optimizes local matches with gaps, making it suitable for finding regions of high similarity within larger sequences.

## Requirements

- Python 3.x

## Usage

### Global Alignment (Needleman-Wunsch)
To perform global alignment, run the following command:

```bash
python AlignGlobal.py <sequence_file1> <sequence_file2>
```

### Local Alignment (Smith-Waterman)
To perform local alignment, run the following command:

```bash
python AlignLocal.py <sequence_file1> <sequence_file2>
```

## Input Format

The sequence files should be in **FASTA format**. Each file should contain a sequence starting with a header line (e.g., `>sequence_name`).

Example of a FASTA file:

```
>seq1
AGCTGTA
```
```
>seq2
AGCTGCA
```

## Example Output

**Global Alignment (Needleman-Wunsch)**:
```
Alignment Score: <score>
Alignment:
AGCTGTA
|||||||
AGCTGCA
```

**Local Alignment (Smith-Waterman)**:
```
Score: <score>
AGCT
||||
AGCT
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
