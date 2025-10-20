#!/bin/bash

#SBATCH --output="catchinrepeats.log"  
#SBATCH --job-name="tgut191alike"  
#SBATCH --time=8:00:00 # walltime  
#SBATCH --cpus-per-task=2 # number of cores  
#SBATCH --mem-per-cpu=8G # memory per CPU core


# Input parameters
database=$1
repeat=$2
names_map=$3

# Validate inputs
if [[ -z "$database" || -z "$repeat" || -z "$names_map" ]]; then
    echo "Usage: $0 <database.fasta> <repeat.fasta> <names_map_file>" >&2
    exit 1
fi

if [[ ! -f "$database" || ! -f "$repeat" || ! -f "$names_map" ]]; then
    echo "Error: One or more input files do not exist" >&2
    exit 1
fi

# Extract base names
db_name=$(basename "$database" .fasta)
repeat_name=$(basename "$repeat" .fasta)

# Output file names
dimer_fasta="${db_name}.dimer.fasta"
dimer_tolid_fasta="${db_name}.dimer.tolid.fasta"
hits_file="${repeat_name}.hits.tsv"
hits_gca="${repeat_name}.hits.GCA.tsv"
hits_tolid="${repeat_name}.hits.tolid.tsv"
hits_dedup="${repeat_name}.hits.tolid.dedup.tsv"

echo "Creating dimer fasta..."
# Create dimer fasta
awk '/^>/ {
    if (header) print header ORS seq seq;
    header=$0;
    seq="";
    next
}
{
    seq=seq $0
}
END {
    print header ORS seq seq
}' "$database" > "$dimer_fasta"

echo "Building BLAST database..."
# Build database from dimer
makeblastdb -in "$dimer_fasta" -dbtype nucl
if [[ $? -ne 0 ]]; then
    echo "Error: makeblastdb failed" >&2
    exit 1
fi

echo "Running BLASTN..."
# Use the dimer database in blastn
blastn \
    -query "$repeat" \
    -db "$dimer_fasta" \
    -strand both \
    -evalue 1e-1 \
    -word_size 7 \
    -reward 1 \
    -penalty -2 \
    -gapopen 4 \
    -gapextend 2 \
    -dust no \
    -soft_masking true \
    -max_hsps 2 \
    -max_target_seqs 5000 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore qlen slen qcovs stitle" \
    -num_threads 8 \
    -out "$hits_file"

if [[ $? -ne 0 ]]; then
    echo "Error: blastn failed" >&2
    exit 1
fi

echo "Processing BLAST results..."
# Replace GCF with GCA
awk -v OFS="\t" '{gsub("GCF", "GCA", $2); print}' "$hits_file" > "$hits_gca"

# Map accessions using names_map
awk 'BEGIN {FS=OFS="\t"}
NR==FNR {
    map[$1]=$2
    next
}
{
    for (gca in map) {
        if (index($2, gca)) {
            gsub(gca, map[gca], $2)
            break
        }
    }
    print
}' "$names_map" "$hits_gca" > "$hits_tolid"

# Deduplicate by column 2
awk '!seen[$2]++' "$hits_tolid" > "$hits_dedup"

echo "Replacing headers in original db fasta..."
# Replace headers in original db fasta using names_map
awk '/^>/ {gsub(/GCF/,"GCA")} {print}' $database > ${db_name}.GCA.fasta
awk 'NR==FNR {map[$1]=$2; next}
/^>/ {
    for (key in map) {
        if (index($0, key)) {
            sub(key, map[key])
            break
        }
    }
}
{print}' "$names_map" "${db_name}.GCA.fasta" > "${db_name}.GCA.tolid.fasta"

echo "Extracting sequences..."
# Extract sequences for each unique hit
cut -f2 "$hits_dedup" | while read -r line; do
    # Extract Tolid prefix for naming
    name=$(echo "$line" | cut -d "_" -f1)
    # Extract the sequence from the FASTA using gfastats
    gfastats ${db_name}.GCA.tolid.fasta "$line" -o fa >> "${name}.blast.fasta"
done

echo "Processing TRF for each species..."
# Process each input FASTA that ends with "blast.fasta"
for file in *blast.fasta; do
    [[ ! -f "$file" ]] && continue
    
    base="${file%.blast.fasta}"
    dir="$base"
    
    echo "Processing $base..."
    mkdir -p "$dir"
    cp -f "$file" "$dir/"
    
    (
        cd "$dir" || exit 1
        
        # 1) Run TRF (sensitive preset)
        echo "  Running TRF..."
        trf "$file" 2 7 7 75 10 40 240 -d -h
        
        # 2) Locate the TRF .dat file
        dat=""
        for d in "$file".*.dat; do
            if [[ -f "$d" ]]; then
                dat="$d"
                break
            fi
        done
        
        if [[ -z "$dat" || ! -s "$dat" ]]; then
            echo "  WARNING: no TRF .dat produced for $file" >&2
            exit 0
        fi
        
        # 3) Extract consensus monomers from the TRF .dat (field 14)
        echo "  Extracting monomers..."
        awk '
            /^Sequence:/ {sid=$2; next}
            /^[0-9]/ {
                printf(">%s_start%s_end%s_len%s\n%s\n", sid, $1, $2, $3, $14)
            }
        ' "$dat" > "${base}.monomers.fa"
        
        # 4) Sort monomers by length (descending)
        echo "  Sorting monomers..."
        awk '
            /^>/ { if (seq) {print hdr "\t" length(seq) "\t" seq}; hdr=$0; seq=""; next }
            { seq=seq $0 }
            END { if (seq) print hdr "\t" length(seq) "\t" seq }
        ' "${base}.monomers.fa" \
        | sort -k2,2nr \
        | awk -F'\t' '{print $1"\n"$3}' > "${base}.monomers.sorted.fa"
        
        # 5) Multiple alignment (auto mode, reverse-complement correction)
        echo "  Running MAFFT alignment..."
        mafft --thread -1 --auto --adjustdirectionaccurately \
            "${base}.monomers.sorted.fa" > "${base}.monomers.sorted.aln.fa" 2>/dev/null
        
        # 6) Build HMM and emit consensus sequence
        echo "  Building HMM and generating consensus..."
        hmmbuild "${base}.monomers.sorted.hmm" "${base}.monomers.sorted.aln.fa" >/dev/null 2>&1
        hmmemit -c "${base}.monomers.sorted.hmm" > "${base}.consensus.fa"
        
        echo "  Done: $PWD"

        #7 Concatenate consensus sequences into a single file
        cat "${base}.consensus.fa" >> ../all_consensus_sequences.fa  
    )
done

find 

echo "Pipeline complete!"

