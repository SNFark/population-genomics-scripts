#!/usr/bin/env bash

# Filter a VCF to AIMs produced by build_aims.py.
#
# Usage:
#   bash filter_vcf_to_aims.sh YOUR_vcf_file.vcf.gz aims_strict
#
# This expects aims_strict.snps to exist in the current directory.

set -euo pipefail

VCF_IN="${1:?ERROR: provide input VCF, e.g. YOUR_vcf_file.vcf.gz}"
AIM_PREFIX="${2:-aims_strict}"

SNP_LIST="${AIM_PREFIX}.snps"
BED_OUT="${AIM_PREFIX}_regions.bed"
VCF_PREFIX="${AIM_PREFIX}"

if [[ ! -f "${SNP_LIST}" ]]; then
  echo "ERROR: SNP list not found: ${SNP_LIST}" >&2
  exit 1
fi

if [[ ! -f "${VCF_IN}" ]]; then
  echo "ERROR: input VCF not found: ${VCF_IN}" >&2
  exit 1
fi

echo "Converting ${SNP_LIST} to BED format..."

awk '{
  n=split($0,a,"_");
  pos=a[n];
  chr=a[1];
  for(i=2;i<n;i++) chr=chr"_"a[i];
  print chr"\t"pos-1"\t"pos;
}' "${SNP_LIST}" > "${BED_OUT}"

echo "Filtering VCF to AIM regions..."

vcftools \
  --gzvcf "${VCF_IN}" \
  --bed "${BED_OUT}" \
  --recode \
  --recode-INFO-all \
  --out "${VCF_PREFIX}"

echo "Compressing filtered VCF..."
bgzip -f "${VCF_PREFIX}.recode.vcf"

echo "Indexing filtered VCF..."
tabix -p vcf "${VCF_PREFIX}.recode.vcf.gz"

echo "Counting AIM SNPs retained in filtered VCF..."
bcftools view -H "${VCF_PREFIX}.recode.vcf.gz" | wc -l

echo "Done."
echo "Filtered VCF: ${VCF_PREFIX}.recode.vcf.gz"
echo "Index:        ${VCF_PREFIX}.recode.vcf.gz.tbi"
