#!/bin/bash

# Script to count repeats in sex chromosomes within 3Mb from both ends
# Input: sex_chr_with_lengths.tsv
# Output: repeat_sex_count_per_accession.txt

# Create output file with header
echo -e "species\taccession\tsex_chr\tstart\tend\ttotal" > repeat_sex_count_per_accession.txt

echo "Processing accessions..."

# Process each line of the TSV file
tail -n +2 sex_chr_with_lengths.tsv | while IFS=$'\t' read -r tolid gca sex_chr acc gff status count found length; do 
    window=3000000
    
    # Check if GFF file exists
    if [[ -f "$gff" ]]; then
        # Count repeats at start (0-3Mb) and end (length-3Mb to length)
        grep "$acc" "$gff" | awk -v window="$window" -v len="$length" \
            -v species="$tolid" -v acc="$acc" -v chr="$sex_chr" '
        {
            at_start = ($4 <= window || $5 <= window)
            at_end = ($4 >= len - window || $5 >= len - window)
            
            if (at_start) count_start++
            if (at_end) count_end++
        }
        END {
            total = count_start + count_end
            printf "%s\t%s\t%s\t%d\t%d\t%d\n", 
                species, acc, chr, 
                count_start+0, 
                count_end+0,
                total+0
        }' >> repeat_sex_count_per_accession.txt
        
        echo "$tolid - $acc processed"
    else
        echo "Warning: File not found: $gff" >&2
    fi
done

echo "Done! Results saved to repeat_sex_count_per_accession.txt"
