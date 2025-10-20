# catching_repeats.sh

A SLURM-based pipeline for identifying and analyzing repetitive DNA sequences using BLAST and Tandem Repeats Finder (TRF).

## Description

This pipeline takes a database of genomic sequences and a set of repeat sequences, identifies matches through BLAST alignment, and processes the results to generate consensus sequences for tandem repeats using TRF, MAFFT, and HMMER.

## Prerequisites

- SLURM job scheduler
- BLAST+ suite (blastn, makeblastdb)
- Tandem Repeats Finder (TRF)
- MAFFT (multiple sequence alignment)
- HMMER (hmmbuild, hmmemit)
- gfastats
- Standard Unix tools (awk, sort, uniq)

## Usage
```bash
sbatch catching_repeats.sh <database.fasta> <repeat.fasta> <names_map_file>
```

### Parameters

1. **database.fasta**: FASTA file containing genomic sequences to search against
2. **repeat.fasta**: FASTA file containing repeat sequences to search for  
3. **names_map_file**: Tab-delimited file mapping sequence identifiers (col1 → col2)

## Pipeline Steps

### 1. Input Validation

- Verifies all three input files exist
- Exits with error message if any file is missing

### 2. Database Preparation

- Creates a dimerized version of the database: `<database>.dimer.fasta`
- Builds BLAST database from the dimerized sequences

### 3. BLAST Search

- Runs blastn alignment of repeat sequences against dimerized database
- Outputs results in tabular format (outfmt 6) to `<repeat>.hits.tsv`

### 4. Hit Processing

- Converts GCF accessions to GCA format → `<repeat>.hits.GCA.tsv`
- Maps accessions using names_map file → `<repeat>.hits.tolid.tsv`  
- Deduplicates hits by target sequence → `<repeat>.hits.tolid.dedup.tsv`

### 5. Sequence Extraction

- Updates database FASTA headers to GCA then ToLID format → `<database>.GCA.tolid.fasta`
- Extracts sequences for each unique hit
- Groups sequences by ToLID prefix (text before first underscore)
- Outputs grouped sequences to `<TolidPrefix>.blast.fasta` files

### 6. Tandem Repeat Analysis

For each extracted sequence file:

- Creates subdirectory and copies sequence file
- Runs TRF with parameters: `trf <file> 2 7 7 75 10 40 240 -d -h`
- Extracts consensus monomers from TRF output
- Sorts monomers by length
- Performs multiple sequence alignment with MAFFT
- Builds HMM profile and generates consensus sequence
- Appends consensus to `all_consensus_sequences.fa`

## Output Files

- `<repeat>.hits.tsv` - Raw BLAST results
- `<repeat>.hits.tolid.dedup.tsv` - Processed and deduplicated hits
- `<database>.GCA.tolid.fasta` - Database with updated headers
- `<TolidPrefix>.blast.fasta` - Extracted sequences per species/group
- `<species>/<base>.monomers.sorted.aln.fa` - Aligned monomers
- `<species>/<base>.consensus.fa` - Consensus sequences
- `all_consensus_sequences.fa` - Combined consensus sequences

## Example
```bash
# Run the pipeline
sbatch catching_repeats.sh genome_db.fasta known_repeats.fasta accession_map.txt

