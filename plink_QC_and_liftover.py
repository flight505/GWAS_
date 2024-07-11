import collections
import logging
import os
import subprocess

import matplotlib.pyplot as plt
import pandas as pd
import yaml

# 1 Update build - using strand files - should all be TOP and GRCh37.
# 2 Update the rsid's - using Illumina rsid files (The batches have ids from Illumina). Running bim_build_and_chip_check.py should confirm
# 3 Update SEX of the IIDs (using the NOPHO sex file)
# 4 Identify duplicate IIDs in the files, write to file and remove with plink2 (Maybe keep first)
# 5 Run a simple QC
# 6 Create VCF files from plink files.
# 7 liftover VCF files to GRCh38
# 8 Create frequency files from liftover VCF
# 9


# TODO
# does the org plink files have a reference? different versions GRCh37/b37 and Hg19  https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19
# if so, use that for liftover
# if not, use the hs37-1kg reference panel for liftover

# We likly need to add some plink commands:
# --keep-allele-order Use this EVERY SINGLE TIME you call a plink command, otherwise the order of Allele1 and Allele2 may (or probably will) be flipped in your data. \
# --allow-no-sex PLINK will default to removing individuals that have unassigned sex, use this to force it to keep them. \
# --snps-only Removes indels from your variant data and keeps only snps \
# --biallelic-only Removes sites with 2+ alleles \

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def run_command(command):
    """Runs a shell command and logs the output."""
    try:
        result = subprocess.run(
            command, shell=True, check=True, text=True, capture_output=True
        )
        logging.info(result.stdout)
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        logging.error(
            f"Command '{e.cmd}' returned non-zero exit status {e.returncode}. Error: {e.stderr}"
        )
        return False, e.stderr


def create_directory(directory):
    """Creates a directory if it does not exist."""
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
            logging.info(f"Created directory: {directory}")
        else:
            logging.info(f"Directory already exists: {directory}")
    except Exception as e:
        logging.error(f"Failed to create directory {directory}: {e}")


def update_build(batch_prefix, strand_file, output_dir):
    """Updates build using strand file and PLINK."""
    create_directory(output_dir)
    temp_dir = os.path.join(output_dir, "temp")
    create_directory(temp_dir)

    temp_prefix = os.path.join(temp_dir, "TEMP_FILE_XX72262628_")
    chr_file = f"{strand_file}.chr"
    pos_file = f"{strand_file}.pos"
    flip_file = f"{strand_file}.flip"

    run_command(f"cut -f 1,2 {strand_file} > {chr_file}")
    run_command(f"cut -f 1,3 {strand_file} > {pos_file}")
    run_command(
        f"awk '{{if ($5==\"-\") print $0}}' {strand_file} | cut -f 1 > {flip_file}"
    )

    commands = [
        f"plink --bfile {batch_prefix} --update-map {chr_file} --update-chr --make-bed --out {temp_prefix}1",
        f"plink --bfile {temp_prefix}1 --update-map {pos_file} --make-bed --out {temp_prefix}2",
        f"plink --bfile {temp_prefix}2 --flip {flip_file} --make-bed --out {temp_prefix}3",
        f"plink --bfile {temp_prefix}3 --extract {pos_file} --make-bed --out {os.path.join(output_dir, 'updated')}",
    ]

    for cmd in commands:
        run_command(cmd)

    # Cleanup temporary files
    for file in [f"{temp_prefix}{i}" for i in range(1, 4)]:
        for ext in [".bed", ".bim", ".fam", ".log", ".nosex"]:
            try:
                os.remove(file + ext)
            except FileNotFoundError:
                pass

    logging.info("Process completed successfully.")
    return os.path.join(output_dir, "updated")


def update_rsids(input_prefix, output_prefix, rsid_data):
    """Updates RS IDs using PLINK."""
    output_dir = os.path.dirname(output_prefix)
    create_directory(output_dir)  # Ensure the directory exists
    success, _ = run_command(
        f"plink2 --bfile {input_prefix} --update-name {rsid_data} --make-bed --out {output_prefix}"
    )
    if success and all(
        os.path.exists(f"{output_prefix}.{ext}") for ext in ["bed", "bim", "fam"]
    ):
        logging.info(f"RS IDs updated successfully for {output_prefix}")
        return output_prefix
    else:
        logging.error(f"Failed to update RS IDs for {output_prefix}")
        return None


def identify_duplicate_iids(input_prefix):
    fam_file = f"{input_prefix}.fam"
    duplicates_file = f"{input_prefix}_duplicate_iids.txt"

    if not os.path.exists(fam_file):
        logging.error(f"File not found: {fam_file}")
        return None

    try:
        with open(fam_file, "r") as file:
            family_data = [line.strip().split()[:2] for line in file.readlines()]
    except Exception as e:
        logging.error(f"Error reading file {fam_file}: {e}")
        return None

    iid_to_fids = collections.defaultdict(list)
    for fid, iid in family_data:
        iid_to_fids[iid].append(fid)

    duplicates = [
        (fid, iid) for iid, fids in iid_to_fids.items() for fid in fids if len(fids) > 1
    ]

    if duplicates:
        try:
            with open(duplicates_file, "w") as file:
                file.write("#FID\tIID\n")
                for fid, iid in duplicates:
                    file.write(f"{fid}\t{iid}\n")
            logging.info(
                f"Duplicate FID-IID pairs identified and written to {duplicates_file}"
            )
        except Exception as e:
            logging.error(f"Error writing duplicates file {duplicates_file}: {e}")
        return duplicates_file
    else:
        logging.warning(f"No duplicate IIDs found in {input_prefix}.fam")
        return None


def remove_duplicate_samples(input_prefix, output_prefix):
    duplicates_file = identify_duplicate_iids(input_prefix)
    if duplicates_file:
        output_dir = os.path.dirname(output_prefix)
        create_directory(output_dir)  # Ensure the directory exists
        remove_command = f"plink2 --bfile {input_prefix} --remove {duplicates_file} --make-bed --out {output_prefix}"
        success, message = run_command(remove_command)
        if success:
            logging.info(f"Duplicate samples removed successfully for {output_prefix}")
            return output_prefix
        else:
            logging.error(
                f"Failed to remove duplicate samples for {output_prefix}: {message}"
            )
            return None
    else:
        logging.error(f"No duplicate IIDs identified for {input_prefix}")
        return None


def check_duplicate_samples(fam_file):
    """Checks for duplicate sample names in the .fam file."""
    with open(fam_file, "r") as file:
        samples = [line.strip().split()[1] for line in file]
    duplicates = [
        item for item, count in collections.Counter(samples).items() if count > 1
    ]
    return duplicates


def plink_bed_to_vcf(
    bed_file,
    bim_file,
    fam_file,
    fasta_file=None,
    out_prefix=None,
    snps_only=True,
    chr_prefix=True,
):
    """Converts PLINK .bed files to VCF format."""
    if out_prefix is None:
        out_prefix = os.path.splitext(bed_file)[0]

    logging.info(f"Converting BED to VCF: {bed_file}, {bim_file}, {fam_file}")

    output_dir = os.path.dirname(out_prefix)
    create_directory(output_dir)  # Ensure the directory exists

    duplicates = check_duplicate_samples(fam_file)
    if duplicates:
        logging.error(f"Duplicate sample names found: {duplicates}")
        return None

    snps_only_option = "--snps-only 'just-acgt'" if snps_only else ""
    chr_prefix_option = "--output-chr chrM" if chr_prefix else ""
    fasta_option = f"--ref-from-fa --fa {fasta_file}" if fasta_file else ""

    commands = [
        f"plink2 --bed {bed_file} --bim {bim_file} --fam {fam_file} --make-pgen --merge-x --sort-vars {snps_only_option} --out sorted",
        f"plink2 --pfile sorted --export vcf id-paste=iid bgz {fasta_option} --out {out_prefix} {chr_prefix_option}",
        "rm sorted.*",
    ]

    for cmd in commands:
        run_command(cmd)

    vcf_path = f"{out_prefix}.vcf.gz"
    if os.path.exists(vcf_path):
        logging.info(f"VCF file generated successfully: {vcf_path}")
        logging.info(f"Starting indexing {vcf_path}")
        run_command(f"bcftools index --threads 12 {vcf_path}")
        logging.info(f"Indexing completed for VCF file: {vcf_path}")
    else:
        logging.error(f"Failed to generate VCF file at {vcf_path}")
        return None

    return vcf_path


def liftover_to_hg38(
    input_vcf, output_bcf, reject_bcf, src_fasta, ref_fasta, chain_file
):
    """Performs liftover to hg38 using bcftools."""
    create_directory(os.path.dirname(output_bcf))  # Ensure output directory exists
    create_directory(os.path.dirname(reject_bcf))  # Ensure reject directory exists

    command = (
        f"bcftools +liftover --no-version -Ou {input_vcf} -- "
        f"-s {src_fasta} -f {ref_fasta} -c {chain_file} --reject {reject_bcf} --write-src | "
        f"bcftools sort -o {output_bcf} -Ob --write-index"
    )
    success, message = run_command(command)
    if success:
        logging.info(f"Liftover completed successfully for {input_vcf}")
        return output_bcf
    else:
        logging.error(f"Liftover failed for {input_vcf}: {message}")
        return None


# Note: Chromosome notation should follow the GRCh38/hg38 notations ('chr#' for autosomal chromosomes and 'chrX' for
# chromosome 23). BUT you are using hg37 for this part as you havent lifted the build . perhaps run this after liftover as a check
def create_frequency_files(vcf_path, dataset_prefix):
    """Creates frequency files from VCF using bcftools."""
    output_vcf = f"{dataset_prefix}_AF.vcf.gz"
    freq_file = f"{dataset_prefix}.frq"
    logging.info(f"Creating frequency file for: {vcf_path}")

    commands = [
        f"bcftools +fill-tags {vcf_path} -Oz --threads 12 -o {output_vcf} -- -t AF",
        f"bcftools query -f '%CHROM\\t%CHROM_%POS_%REF_%ALT\\t%REF\\t%ALT\\t%INFO/AF\\n' {output_vcf} | sed '1i\\\nCHR\\tSNP\\tREF\\tALT\\tAF' > {freq_file}",
    ]

    for cmd in commands:
        success, message = run_command(cmd)
        if not success:
            logging.error(f"Command failed: {cmd}")
            logging.error(message)
            return None

    return freq_file


def analyze_frequency_files(
    freq_file, ref_freq_file, dataset_prefix, af_diff_limit=0.1
):
    """Analyzes frequency files to compare allele frequencies."""
    try:
        # Ensure the delimiter is correctly specified, assuming tab-delimited files
        freq_data = pd.read_csv(freq_file, sep="\t")
        ref_data = pd.read_csv(ref_freq_file, sep="\t")

        merged_data = pd.merge(
            freq_data, ref_data, on="SNP", suffixes=(".chip", ".ref")
        )
        merged_data["AF_diff"] = abs(merged_data["AF.chip"] - merged_data["AF.ref"])
        high_diff_data = merged_data[merged_data["AF_diff"] > af_diff_limit]

        plt.figure(figsize=(10, 6))
        plt.scatter(
            merged_data["AF.ref"],
            merged_data["AF.chip"],
            c="blue",
            label="Within Threshold",
        )
        plt.scatter(
            high_diff_data["AF.ref"],
            high_diff_data["AF.chip"],
            c="red",
            label="Above Threshold",
        )
        plt.xlabel("Reference Allele Frequency")
        plt.ylabel("Chip Allele Frequency")
        plt.title(f"Allele Frequency Comparison for {dataset_prefix}")
        plt.legend()
        plt.grid(True)

        output_image_path = f"{dataset_prefix}_allele_frequency_comparison.png"
        plt.savefig(output_image_path)
        logging.info(f"Plot saved to {output_image_path}")

        # Display the plot if running in an interactive environment
        plt.show()

        return high_diff_data
    except pd.errors.ParserError as e:
        logging.error(f"Error parsing CSV file: {e}")
        return None
    except Exception as e:
        logging.error(f"Error in analyzing frequency files: {e}")
        return None


def load_config(config_path) -> dict:
    """Loads configuration from a YAML file."""
    with open(config_path, "r") as file:
        return yaml.load(file, Loader=yaml.FullLoader)


def process_batches(config):
    """Main processing workflow for genomic data."""
    batches = config["batches"]
    strand_files = config["strand_files"]
    rsid_files = config["rsid_files"]
    human_g1k_v37_path = config["human_g1k_v37_path"]
    GRCh37_path = config["GRCh37_path"]
    GRCh38_path = config["GRCh38_path"]
    chain_file = config["chain_file"]

    updated_batches = {
        batch: update_build(path, strand_files[batch], f"data/1_updated/{batch}")
        for batch, path in batches.items()
    }

    updated_rsids = {
        batch: update_rsids(
            updated_path, f"data/2_updated_rsids/{batch}/{batch}", rsid_files[batch]
        )
        for batch, updated_path in updated_batches.items()
        if updated_path
    }

    removed_duplicates = {
        batch: remove_duplicate_samples(
            updated_path, f"data/3_removed_duplicates/{batch}/{batch}_no_duplicates"
        )
        for batch, updated_path in updated_rsids.items()
        if updated_path
    }

    vcf_outputs = {
        batch: plink_bed_to_vcf(
            f"{path}.bed",
            f"{path}.bim",
            f"{path}.fam",
            fasta_file=GRCh37_path,
            out_prefix=f"data/4_vcf_outputs/{batch}/{batch}",
        )
        for batch, path in removed_duplicates.items()
        if path
    }

    for batch, vcf_file in vcf_outputs.items():
        if vcf_file:
            create_frequency_files(vcf_file, f"data/4_vcf_outputs/{batch}")

    frequency_plots = {
        batch: analyze_frequency_files(
            vcf_file,
            f"data/reference_frequency/{batch}.frq",
            f"data/4_vcf_outputs/{batch}",
            0.1,
        )
        for batch, vcf_file in vcf_outputs.items()
        if vcf_file
    }

    frequency_plots = {
        batch: plot_data
        for batch, plot_data in frequency_plots.items()
        if plot_data is not None
    }
    logging.info(f"Frequency plots data: {frequency_plots}")

    liftover_outputs = {
        batch: liftover_to_hg38(
            vcf_file,
            os.path.join(f"data/5_liftover/{batch}", f"{batch}.bcf"),
            os.path.join(f"data/5_liftover/{batch}", f"{batch}_reject.bcf"),
            human_g1k_v37_path,
            GRCh38_path,
            chain_file,
        )
        for batch, vcf_file in vcf_outputs.items()
        if vcf_file
    }
    logging.info(f"Liftover outputs: {liftover_outputs}")


def main():
    config = load_config("config.yaml")
    process_batches(config)


if __name__ == "__main__":
    main()
