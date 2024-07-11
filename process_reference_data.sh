#!/bin/bash

# Written by Jesper Vang @DTU 2024
# uses https://savannah.gnu.org/projects/parallel/ 
# install with: brew install parallel
# enter once in terminal to remove the citation request: parallel --citation

#	1.	Confirm that multiallelic sites (if present) in your reference data files are decomposed. If they are not, the script will split the multiallelic sites into biallelic records.
#	2.	Renaming chromosomes in the VCF files. Confirms that the chromosome notation in your reference data files follows the GRCh38/h38 notations as 'chr#' for autosomal,
#		'chrX' for chromosome 23 and 'chrM' for mitochondrial sites.
#	3.	Generate a tab-delimited file of the reference data allele frequencies, one line per variant, with columns CHR, SNP (as
#		CHR_POS_REF_ALT), REF, ALT, AF (including the header line).
#	4.	Querying the required fields from each VCF file and appending them to the ref_data.frq file.

# The final ref_data.frq file will be a single tab-delimited file containing the allele frequencies for all chromosomes.

# Define the reference file path
REFERENCE_FILE="/Users/jespervang/Projects/MTX/GWAS_paper/data/ref_37/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz"

# Decompose multiallelic sites in parallel
echo "Decomposing multiallelic sites..."
seq 1 22 | parallel -j 12 'bcftools norm -m -any '"${REFERENCE_FILE}"' -Oz -o ref_data_split_multiallelic_chr{}.vcf.gz'
bcftools norm -m -any "${REFERENCE_FILE}" -Oz -o ref_data_split_multiallelic_chrX.vcf.gz

# Generate a chromosome renaming file
echo "Renaming chromosomes..."
for CHR in {1..22} X; do
    echo ${CHR} chr${CHR}
done > chr_names.txt

# Rename chromosomes in VCF files in parallel
seq 1 22 | parallel -j 12 'bcftools annotate --rename-chrs chr_names.txt ref_data_split_multiallelic_chr{}.vcf.gz -Oz -o ref_data_chr{}.vcf.gz'
bcftools annotate --rename-chrs chr_names.txt ref_data_split_multiallelic_chrX.vcf.gz -Oz -o ref_data_chrX.vcf.gz

# Calculate AF for each chromosome VCF file in parallel
echo "Generating allele frequency file..."
seq 1 22 | parallel -j 12 'bcftools +fill-tags ref_data_chr{}.vcf.gz -Oz -o ref_data_AF_chr{}.vcf.gz -- -t AF'
bcftools +fill-tags ref_data_chrX.vcf.gz -Oz -o ref_data_AF_chrX.vcf.gz -- -t AF

# Generate a tab-delimited header
echo -e 'CHR\tSNP\tREF\tALT\tAF' > ref_data.frq

# Query the required fields from the reference VCF files, append to the allele frequency file in parallel
seq 1 22 | parallel -j 12 'bcftools query -f "%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\n" ref_data_chr{}.vcf.gz >> ref_data.frq'
bcftools query -f "%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\n" ref_data_chrX.vcf.gz >> ref_data.frq

echo "Process completed."