#!/usr/bin/env bash
set -euo pipefail

# Usage: ./check_sex_chr_in_gff.sh sex_chr_list.tsv [base_dir] [out_tsv]
list="${1:-sex_chr_list.tsv}"
basedir="${2:-.}"
out="${3:-sex_chr_in_gff.tsv}"

printf "tolid\tchr\taccession\tgff_path\tstatus\tmatch_count\n" > "$out"

# Accept tabs or spaces; ignore extra trailing columns if present
while IFS=$' \t' read -r tolid chr acc _; do
  # Skip blank/comment lines
  [[ -z "${tolid:-}" || "$tolid" =~ ^# ]] && continue

  gff="${basedir}/${tolid}/out/second_pass/${tolid}.merged.gff"

  if [[ ! -f "$gff" ]]; then
    printf "%s\t%s\t%s\t%s\tfile_missing\t0\n" "$tolid" "$chr" "$acc" "$gff" >> "$out"
    continue
  fi

  # Check if accession is the seqid (column 1) in the GFF
  seqid_count=$(awk -F'\t' -v a="$acc" '($0!~/^#/ && $1==a){c++} END{print c+0}' "$gff")

  if (( seqid_count > 0 )); then
    printf "%s\t%s\t%s\t%s\tseqid\t%d\n" "$tolid" "$chr" "$acc" "$gff" "$seqid_count" >> "$out"
    continue
  fi

  # Otherwise, look anywhere in the file (attributes, etc.)
  any_count=$(grep -F -c -- "$acc" "$gff" || true)
  if (( any_count > 0 )); then
    printf "%s\t%s\t%s\t%s\tfound_elsewhere\t%d\n" "$tolid" "$chr" "$acc" "$gff" "$any_count" >> "$out"
  else
    printf "%s\t%s\t%s\t%s\tabsent\t0\n" "$tolid" "$chr" "$acc" "$gff" >> "$out"
  fi
done < "$list"

echo "Wrote: $out"
