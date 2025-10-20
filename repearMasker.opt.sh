#!/bin/bash
#SBATCH --job-name=bird_repeatmasker
#SBATCH --time=3-00:00:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --partition=vgl_a
#SBATCH --account=vgl_condo_bank
#SBATCH --output=bird_repeatmasker_%j.out
#SBATCH --error=bird_repeatmasker_%j.err

set -euo pipefail

# Input arguments
genome_list="$1"
lib="$2"

# Export the repeat library for use in subprocesses
export lib

process_chromosome() {
    local chromosome="$1"
    local fasta="$2"
    local genome="$3"

    echo "[INFO] Processing chromosome: $chromosome"

    out_dir="$genome/out/second_pass/chromosomes"
    mkdir -p "$out_dir"


    chr_fasta="${genome}/out/chromosomes/${chromosome}.fasta"
#    echo "[INFO] Extracting chromosome with gfastats: $chromosome"
#    gfastats "$fasta" "$chromosome" -o "$out_fasta"

    echo "[INFO] Running RepeatMasker on $chr_fasta"
    RepeatMasker -dir "$out_dir" -lib "$lib" "$chr_fasta" -nolow -no_is -pa 8 -gff
    echo "[INFO] Completed: $chromosome"
}

export -f process_chromosome

process_genome() {
    local genome="$1"
    echo "[INFO] Processing genome: $genome"

    fasta=$(find "$genome" -maxdepth 1 -name '*.fna' | head -n 1)
    chr_list=$(find "$genome" -maxdepth 1 -name '*chromosomes.lst' | head -n 1)

    if [[ ! -f "$fasta" || ! -f "$chr_list" ]]; then
        echo "[WARNING] Missing fasta or chromosome list in $genome. Skipping."
        return
    fi

    export fasta
    export genome

    # 4 concurrent chromosomes per genome, each using 8 threads
    cat "$chr_list" | parallel --jobs 4 "process_chromosome {} \"$fasta\" \"$genome\""
}

export -f process_genome

echo "[INFO] Starting genome batch processing..."

# 4 concurrent genome folders, each using 32 threads
awk '{print $1}' "$genome_list" | parallel --jobs $SLURM_NTASKS --link --joblog parallel.log "srun --exclusive -N1 -n1 -c32 bash -c 'process_genome {}'"

echo "[INFO] All genome processing completed."
