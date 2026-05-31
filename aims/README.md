# AIM panel builder

## Overview

This repository contains scripts for building ancestry-informative marker (AIM) panels from PLINK `.frq` files and extracting those AIMs from a VCF.

The primary goal is to identify species-diagnostic markers that are consistent across multiple geographic regions.

Rather than selecting markers from a single population pair, the approach requires allele-frequency differences to be replicated across all specified regions and in the same direction. This reduces the influence of local geographic structure and population-specific allele-frequency variation, enriching for markers that are broadly species-diagnostic.

The resulting AIM panel is therefore intended to capture species-level differentiation while minimizing the contribution of regional population structure.

## Files

- build_aims.py
- filter_vcf_to_aims.sh

## AIM selection criteria

A SNP is retained only if it:

1. is present in both species in all specified regions,
2. passes the minimum chromosome-count threshold in every region,
3. is biallelic,
4. shows the same allele-frequency difference direction across all regions,
5. exceeds the specified minimum absolute allele-frequency difference in every region.

By requiring consistency across regions, the resulting marker set is less likely to reflect geographic population structure and more likely to represent species-wide differentiation, useful for downstream analyses such as PCA and Admixture analyses.

## Input files

Expected format:

<REGION>_ADF_<SPECIES>.frq

Example:

AIS_ADF_sel.frq
AIS_ADF_cin.frq
SIWSS_ADF_sel.frq
SIWSS_ADF_cin.frq
FARG_ADF_sel.frq
FARG_ADF_cin.frq

## Example

python build_aims.py \
  --regions AIS SIWSS FARG \
  --species-a sel \
  --species-b cin \
  --threshold 0.8 \
  --min-nchr 6 \
  --out-prefix aims_strict

## Outputs

aims_strict.tsv
aims_strict.snps

## VCF filtering

bash filter_vcf_to_aims.sh YOUR_vcf_file.vcf.gz aims_strict

This converts the SNP list into a BED file, filters the VCF to AIM positions, compresses and indexes the filtered VCF, and reports the number of retained AIMs.

The resulting VCF can be used for downstream PCA, ADMIXTURE, PLINK conversion, or other ancestry analyses.
