#!/bin/sh


#SBATCH --job-name=downloding    # Job name
#SBATCH --output=%x_%j.slurm


file=$1


awk -F'\t' '{print $1 "\t" $2}' birds.txt | while read -r tolid accession;do
echo "folder $tolid created"
mkdir $tolid;
echo "downloading data for $tolid"
(cd $tolid &&
 datasets download genome accession $accession --include genome &&
 unzip ncbi_dataset.zip &&
 mv ncbi_dataset/data/$accession/* . &&
 rm -r ncbi_dataset)
done
