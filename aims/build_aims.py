#!/usr/bin/env python3
"""
Build a cross-region ancestry-informative marker (AIM) panel from PLINK .frq files.

The script identifies biallelic SNPs that:
  1. are present in all specified regions for both species,
  2. pass a minimum chromosome-count threshold in every region/species file,
  3. show the same allele-frequency direction in every region,
  4. exceed a minimum absolute allele-frequency difference in every region.

Expected input filename format:
  <REGION>_ADF_<SPECIES>.frq

Example usage:
  python build_aims.py \
    --regions AIS SIWSS FARG \
    --species-a sel \
    --species-b cin \
    --threshold 0.8 \
    --min-nchr 6 \
    --out-prefix aims_strict
"""

import argparse
import sys
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a cross-region AIM panel from PLINK .frq files using delta AF."
    )
    parser.add_argument("--regions", nargs="+", required=True, help="Region names, e.g. AIS SIWSS FARG")
    parser.add_argument("--species-a", required=True, help="Short code for species A in filenames, e.g. sel")
    parser.add_argument("--species-b", required=True, help="Short code for species B in filenames, e.g. cin")
    parser.add_argument("--threshold", type=float, default=0.8, help="Minimum absolute delta AF required in every region")
    parser.add_argument("--min-nchr", type=int, default=2, help="Minimum N_CHR required per species per region per SNP")
    parser.add_argument("--out-prefix", default="aims", help="Output prefix")
    return parser.parse_args()


def read_frq(filepath):
    """
    Read a PLINK .frq file with columns like:
      CHROM POS N_ALLELES N_CHR {ALLELE:FREQ}

    Returns a DataFrame with CHROM, POS, SNP_ID, N_ALLELES, N_CHR, FREQS.
    FREQS is a dictionary of allele -> frequency.
    """
    records = []

    with open(filepath, "r") as f:
        header = f.readline().strip().split()

        expected = {"CHROM", "POS", "N_ALLELES", "N_CHR"}
        if not expected.issubset(set(header)):
            raise ValueError(
                f"{filepath} is missing required columns. "
                f"Expected at least {sorted(expected)}, found {header}"
            )

        for line in f:
            if not line.strip():
                continue

            parts = line.strip().split()
            chrom = parts[0]
            pos = parts[1]

            try:
                n_alleles = int(parts[2])
                n_chr = int(parts[3])
            except ValueError:
                continue

            if n_chr == 0:
                continue

            freqs = {}
            for af in parts[4:]:
                if ":" not in af:
                    continue

                allele, freq = af.split(":", 1)

                if freq == "-nan":
                    continue

                try:
                    freqs[allele] = float(freq)
                except ValueError:
                    continue

            if not freqs:
                continue

            records.append({
                "CHROM": chrom,
                "POS": pos,
                "SNP_ID": f"{chrom}_{pos}",
                "N_ALLELES": n_alleles,
                "N_CHR": n_chr,
                "FREQS": freqs
            })

    return pd.DataFrame(records)


def main():
    args = parse_args()

    regions = args.regions
    spA = args.species_a
    spB = args.species_b
    threshold = args.threshold
    min_nchr = args.min_nchr
    out_prefix = args.out_prefix

    region_data = {}

    for region in regions:
        file_a = f"{region}_ADF_{spA}.frq"
        file_b = f"{region}_ADF_{spB}.frq"

        try:
            df_a = read_frq(file_a)
            df_b = read_frq(file_b)
        except FileNotFoundError as e:
            sys.exit(f"ERROR: file not found: {e.filename}")
        except Exception as e:
            sys.exit(f"ERROR reading files for region {region}: {e}")

        if df_a.empty:
            sys.exit(f"ERROR: {file_a} contains no usable records.")
        if df_b.empty:
            sys.exit(f"ERROR: {file_b} contains no usable records.")

        region_data[region] = {
            spA: df_a.set_index("SNP_ID"),
            spB: df_b.set_index("SNP_ID")
        }

    common_snps = None
    for region in regions:
        snps_here = set(region_data[region][spA].index) & set(region_data[region][spB].index)
        common_snps = snps_here if common_snps is None else common_snps & snps_here

    if not common_snps:
        sys.exit("ERROR: no common SNPs found across all regions.")

    passed = []
    total_checked = 0
    skipped_low_n = 0
    skipped_multiallelic = 0
    skipped_missing_ref = 0
    skipped_direction = 0
    skipped_threshold = 0

    first_region = regions[0]

    for snp in sorted(common_snps):
        total_checked += 1
        first_record = region_data[first_region][spA].loc[snp]

        if first_record["N_ALLELES"] != 2:
            skipped_multiallelic += 1
            continue

        ref_alleles = sorted(first_record["FREQS"].keys())
        if len(ref_alleles) != 2:
            skipped_multiallelic += 1
            continue

        ref_allele = ref_alleles[0]
        chrom = first_record["CHROM"]
        pos = first_record["POS"]

        per_region = []
        valid = True

        for region in regions:
            rec_a = region_data[region][spA].loc[snp]
            rec_b = region_data[region][spB].loc[snp]

            if rec_a["N_CHR"] < min_nchr or rec_b["N_CHR"] < min_nchr:
                skipped_low_n += 1
                valid = False
                break

            if rec_a["N_ALLELES"] != 2 or rec_b["N_ALLELES"] != 2:
                skipped_multiallelic += 1
                valid = False
                break

            alleles_a = set(rec_a["FREQS"].keys())
            alleles_b = set(rec_b["FREQS"].keys())

            if ref_allele not in alleles_a or ref_allele not in alleles_b:
                skipped_missing_ref += 1
                valid = False
                break

            af_a = rec_a["FREQS"][ref_allele]
            af_b = rec_b["FREQS"][ref_allele]
            daf = af_a - af_b

            per_region.append({
                "region": region,
                "af_a": af_a,
                "af_b": af_b,
                "daf": daf,
                "nchr_a": rec_a["N_CHR"],
                "nchr_b": rec_b["N_CHR"]
            })

        if not valid:
            continue

        dafs = [x["daf"] for x in per_region]
        all_positive = all(d > 0 for d in dafs)
        all_negative = all(d < 0 for d in dafs)

        if not (all_positive or all_negative):
            skipped_direction += 1
            continue

        min_abs_daf = min(abs(d) for d in dafs)
        mean_abs_daf = sum(abs(d) for d in dafs) / len(dafs)

        if min_abs_daf < threshold:
            skipped_threshold += 1
            continue

        row = {
            "CHROM": chrom,
            "POS": pos,
            "SNP_ID": snp,
            "REF_ALLELE": ref_allele,
            "MIN_ABS_DAF": round(min_abs_daf, 6),
            "MEAN_ABS_DAF": round(mean_abs_daf, 6),
            "DIRECTION": f"{spA}_gt_{spB}" if all_positive else f"{spB}_gt_{spA}"
        }

        for x in per_region:
            region = x["region"]
            row[f"AF_{spA}_{region}"] = round(x["af_a"], 6)
            row[f"AF_{spB}_{region}"] = round(x["af_b"], 6)
            row[f"DAF_{region}"] = round(x["daf"], 6)
            row[f"NCHR_{spA}_{region}"] = int(x["nchr_a"])
            row[f"NCHR_{spB}_{region}"] = int(x["nchr_b"])

        passed.append(row)

    passed_df = pd.DataFrame(passed)

    if not passed_df.empty:
        passed_df = passed_df.sort_values(
            by=["MIN_ABS_DAF", "MEAN_ABS_DAF"],
            ascending=False
        )

    tsv_out = f"{out_prefix}.tsv"
    snp_out = f"{out_prefix}.snps"

    passed_df.to_csv(tsv_out, sep="\t", index=False)

    with open(snp_out, "w") as f:
        if not passed_df.empty:
            for snp in passed_df["SNP_ID"]:
                f.write(f"{snp}\n")

    print("=== SUMMARY ===")
    print(f"Regions: {', '.join(regions)}")
    print(f"Species A: {spA}")
    print(f"Species B: {spB}")
    print(f"Threshold: {threshold}")
    print(f"Min N_CHR: {min_nchr}")
    print(f"Total SNPs checked: {total_checked}")
    print(f"Passed: {len(passed_df)}")
    print(f"Skipped low N_CHR: {skipped_low_n}")
    print(f"Skipped multiallelic: {skipped_multiallelic}")
    print(f"Skipped missing reference allele: {skipped_missing_ref}")
    print(f"Skipped inconsistent direction: {skipped_direction}")
    print(f"Skipped threshold: {skipped_threshold}")
    print(f"Output table: {tsv_out}")
    print(f"Output SNP list: {snp_out}")


if __name__ == "__main__":
    main()
