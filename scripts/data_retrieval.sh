#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-8}"

META="metadata/metadata.tsv"
SRA_DIR="data/sra"
FASTQ_DIR="data/raw_fastq"

mkdir -p "$SRA_DIR" "$FASTQ_DIR"

# extracting SRR accessions from metadata(skipping the header)
mapfile -t SRRS < <(cut -f2 "$META" | tail -n +2)

echo "Downloading ${#SRRS[@]} runs..."

for srr in "${SRRS[@]}"; do
  echo "==> ${srr}"

  # downloading .sra
  prefetch -O "$SRA_DIR" "$srr"

  # converting to FASTQ (using --split-files to output _1/_2 if paired end)
  SRA_FILE="${SRA_DIR}/${srr}/${srr}.sra"
  fasterq-dump --threads "$THREADS" --split-files -O "$FASTQ_DIR" "$SRA_FILE"

  # compressing FASTQs for this SRR
  gzip -f "$FASTQ_DIR/${srr}"*.fastq

  # removing .sra to save space
  rm -f "$SRA_FILE"
done

echo
echo "Download complete. FASTQs:"
ls -lh "$FASTQ_DIR" | head
