import argparse
import logging
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests
from tqdm import tqdm

"""
BIM Build and Genotyping Chip Check Tool

This script reads a .bim file, fetches SNP data from Ensembl, 
cross-checks the SNP positions against GRCh38 and GRCh37 builds, 
and analyzes genotyping chips.

Usage Instructions
1. Save the Script: Save the script as bim_build_and_chip_check.py.
2. Run the Script: Use the command line to run the script with the required .bim file path
and optional number of SNP samples.

Command Line Examples
python bim_build_and_chip_check.py /path/to/your/file.bim
Specifying the Number of SNP Samples:
python bim_build_and_chip_check.py /path/to/your/file.bim --num_samples 500

Output
Log File: Created in the same directory as the .bim file with the suffix _analysis.log.
Log Content: Includes total variants, 'rs' variants processed, number of SNPs sampled,
confidence levels for GRCh38 and GRCh37, and top genotyping chips.
"""


def setup_logging(bim_file_path):
    """
    Set up logging to log both to a file and to the console.

    Parameters:
    bim_file_path (str): Path to the .bim file. The log file will be created in the same directory with a suffix '_analysis.log'.
    """
    log_filename = bim_file_path.replace(".bim", "_analysis.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
        handlers=[
            logging.FileHandler(log_filename, mode="w"),  # Overwrite the log file
            logging.StreamHandler(),
        ],
    )


def read_bim_file(file_path):
    snp_data = []
    total_variants = 0
    with open(file_path, "r") as file:
        for line in file:
            parts = line.strip().split()
            total_variants += 1
            if len(parts) < 6 or not parts[1].startswith("rs"):
                continue
            chromosome, variant_id, morgan_position, bp_coordinate, allele1, allele2 = (
                parts
            )
            snp_data.append(
                {
                    "chromosome": chromosome,
                    "variant_id": variant_id,
                    "bp_coordinate": int(bp_coordinate),
                    "allele1": allele1,
                    "allele2": allele2,
                }
            )
    logging.info(f"Total variants loaded from .bim file: {total_variants}")
    logging.info(f"Total 'rs' variants processed: {len(snp_data)}")
    return pd.DataFrame(snp_data)


def fetch_snp_data(snp_id, servers, ext, headers):
    snp_info = {
        "rs_id": snp_id,
        "chromosome": None,
        "position_grch38": None,
        "position_grch37": None,
        "alleles": None,
        "genotyping_chips": [],
    }
    for assembly, server in servers.items():
        try:
            response = requests.get(
                f"{server}{ext}/{snp_id}?genotyping_chips=1", headers=headers
            )
            response.raise_for_status()  # Raises an HTTPError for bad responses
            data = response.json()
            mappings = data.get("mappings", [])
            if not mappings:
                logging.warning(f"No mapping data available for SNP {snp_id}")
            for mapping in mappings:
                if mapping["assembly_name"] == "GRCh38":
                    snp_info["position_grch38"] = mapping["start"]
                elif mapping["assembly_name"] == "GRCh37":
                    snp_info["position_grch37"] = mapping["start"]
                snp_info["chromosome"] = mapping["seq_region_name"]
                snp_info["alleles"] = mapping.get("alleles")
            snp_info["genotyping_chips"] = data.get("genotyping_chips", [])
        except requests.HTTPError as e:
            logging.error(f"HTTP error occurred: {e}")
        except KeyError as e:
            logging.error(f"Missing expected data key in response: {e}")
        except Exception as e:
            logging.error(f"An unexpected error occurred: {e}")

    return snp_info


def fetch_snp_data_from_ensembl(snp_ids):
    servers = {
        "GRCh37": "https://grch37.rest.ensembl.org",
        "GRCh38": "https://rest.ensembl.org",
    }
    ext = "/variation/human"
    headers = {"Content-Type": "application/json"}
    with ThreadPoolExecutor(max_workers=10) as executor:
        results = list(
            tqdm(
                executor.map(
                    lambda snp_id: fetch_snp_data(snp_id, servers, ext, headers),
                    snp_ids,
                ),
                total=len(snp_ids),
                desc="Fetching SNP data",
            )
        )
    return pd.DataFrame(results)


def cross_check_snp_data(plink_data, snp_data):
    grch38_matches = 0
    grch37_matches = 0
    mismatch_stats = {"strand_mismatches": 0, "potential_inversions": 0}
    plink_snps = set(plink_data["variant_id"].values)

    for snp in snp_data.itertuples():
        snp_id = str(snp.rs_id)
        if snp_id in plink_snps:
            plink_snp = plink_data[plink_data["variant_id"] == snp_id].iloc[0]
            plink_pos = plink_snp["bp_coordinate"]
            plink_alleles = {plink_snp["allele1"], plink_snp["allele2"]}
            fetched_alleles = set(snp.alleles.split("/") if snp.alleles else [])

            if plink_pos == snp.position_grch38:
                grch38_matches += 1
            if plink_pos == snp.position_grch37:
                grch37_matches += 1

            # Check for mismatches and potentially for inversion issues
            if plink_alleles != fetched_alleles:
                mismatch_stats["strand_mismatches"] += 1
                if plink_alleles != {a[::-1] for a in fetched_alleles}:
                    mismatch_stats["potential_inversions"] += 1
                logging.warning(
                    f"Potential strand mismatch or inversion for SNP {snp_id}"
                )

    return grch38_matches, grch37_matches, mismatch_stats


def analyze_genotyping_chips(snp_data):
    chip_counts = {}
    for snp in snp_data.itertuples():
        for chip in snp.genotyping_chips:
            if chip in chip_counts:
                chip_counts[chip] += 1
            else:
                chip_counts[chip] = 1
    sorted_chips = sorted(chip_counts.items(), key=lambda item: item[1], reverse=True)[
        :3
    ]
    total_chips = sum(chip_counts.values())
    top_chips_confidence = [
        (chip, count / total_chips * 100) for chip, count in sorted_chips
    ]
    return top_chips_confidence


def main(bim_file_path, num_samples):
    setup_logging(bim_file_path)
    logging.info("Starting SNP analysis script.")
    logging.info(
        "This script reads a .bim file, fetches SNP data from Ensembl, cross-checks the SNP positions, and analyzes genotyping chips."
    )
    logging.info("Steps involved:")
    logging.info("1. Read the .bim file and extract SNP data.")
    logging.info("2. Fetch SNP data from Ensembl using multithreading.")
    logging.info("3. Cross-check SNP positions with GRCh37 and GRCh38 assemblies.")
    logging.info("4. Analyze the top genotyping chips used for the SNPs.")

    plink_data = read_bim_file(bim_file_path)
    snps = plink_data["variant_id"].sample(n=num_samples).tolist()
    logging.info(f"Number of 'rs' SNPs sampled: {len(snps)}")

    snp_data = fetch_snp_data_from_ensembl(snps)
    if not snp_data.empty:
        grch38_matches, grch37_matches, mismatch_stats = cross_check_snp_data(
            plink_data, snp_data
        )
        logging.info(f"Confidence for GRCh38: {grch38_matches / len(snps) * 100:.2f}%")
        logging.info(f"Confidence for GRCh37: {grch37_matches / len(snps) * 100:.2f}%")
        logging.info(f"Total strand mismatches: {mismatch_stats['strand_mismatches']}")
        logging.info(
            f"Total potential inversions: {mismatch_stats['potential_inversions']}"
        )

        top_chips = analyze_genotyping_chips(snp_data)
        logging.info("Top Genotyping Chips and their confidence levels:")
        for chip, confidence in top_chips:
            logging.info(f"{chip}: {confidence:.2f}%")
    else:
        logging.info("No SNP data fetched from Ensembl.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SNP analysis script")
    parser.add_argument("bim_file_path", type=str, help="Path to the .bim file")
    parser.add_argument(
        "--num_samples",
        type=int,
        default=250,
        help="Number of SNP samples to check (default: 250)",
    )
    args = parser.parse_args()
    main(args.bim_file_path, args.num_samples)
