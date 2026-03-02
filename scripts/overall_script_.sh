
#01: getting reads
# writing a script and saving it separately
cat > scripts/data_retrieval.sh << 'EOF'
#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-8}"

META="metadata/metadata.tsv"
SRA_DIR="data/sra"
FASTQ_DIR="data/raw_fastq"

mkdir -p "$SRA_DIR" "$FASTQ_DIR"

# extracting SRR accessions from metadata (skipping the header)
mapfile -t SRRS < <(cut -f2 "$META" | tail -n +2)

echo "Downloading ${#SRRS[@]} runs..."

for srr in "${SRRS[@]}"; do
  echo "==> ${srr}"

  # downloading .sra
  prefetch -O "$SRA_DIR" "$srr"

  # converting to FASTQ
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
EOF
chmod +x scripts/data_retrieval.sh
#run script made
./scripts/data_retrieval.sh

#sanity checking:
ls data/raw_fastq | wc -l
#single end or paired end
ls data/raw_fastq | head
#confirming all SRRs are listed in metadata
cut -f2 metadata/metadata.tsv | tail -n +2 | while read -r srr; do
  ls data/raw_fastq/${srr}*.fastq.gz >/dev/null && echo "OK $srr" || echo "MISSING $srr"
done

# QC step prior to quantification

#fast qc was previously downloaded, for multiqc, I decided to use docker to pull an image

#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-8}"
IN_DIR="data/raw_fastq"
OUT_FASTQC="results/qc/fastqc"
OUT_MULTI="results/qc/multiqc"

mkdir -p "$OUT_FASTQC" "$OUT_MULTI"

# Runing FastQC locally (since previously downloaded)
fastqc -t "$THREADS" -o "$OUT_FASTQC" "$IN_DIR"/*.fastq.gz

# running MultiQC in docker
docker run --rm -t \
  -u "$(id -u)":"$(id -g)" \
  -v "$(pwd)":"$(pwd)" -w "$(pwd)" \
  multiqc/multiqc \
  multiqc "$OUT_FASTQC" -o "$OUT_MULTI"

echo "QC complete."
echo "FastQC outputs: $OUT_FASTQC"
echo "MultiQC report: $OUT_MULTI/multiqc_report.html"

'''
docker run --rm -t \
  -u "$(id -u)":"$(id -g)" \
  -v "$(pwd)":"$(pwd)" -w "$(pwd)" \
  multiqc/multiqc \
  multiqc results/qc/fastqc -o results/qc/multiqc
'''
# quantification

#getting the reference genome, transcriptome
cd ref

# S. cerevisiae (S288C / R64) reference from NCBI
wget -O transcripts.fa.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_rna.fna.gz"
wget -O genome.fa.gz      "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz"
wget -O genes.gtf.gz      "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz"

gunzip -f *.gz
cd ..
#making the decoy aware gentrome needed
cat ref/transcripts.fa ref/genome.fa > ref/gentrome.fa
grep "^>" ref/genome.fa | sed 's/^>//' | cut -d' ' -f1 > ref/decoys.txt

#building an index using docker image
mkdir -p ref/salmon_index
docker run --rm -t \
  -u "$(id -u)":"$(id -g)" \
  -v "$(pwd)":"$(pwd)" -w "$(pwd)" \
  combinelab/salmon \
  salmon index -t ref/gentrome.fa -d ref/decoys.txt -i ref/salmon_index -k 31

#quantifying each sample, catering to the single end reads, using default parameters found for --fldMean and --fldSD


#sanity check for output
ls results/quants | wc -l
head -n 5 results/quants/$(ls results/quants | head -n 1)/quant.sf
#recording one quant stat
for srr in $(cut -f2 metadata/metadata.tsv | tail -n +2); do
  echo "== $srr =="
  grep -E "Mapping rate|Percent mapped" -n results/quants/${srr}/logs/salmon_quant.log || true
done

#version check: docker run --rm multiqc/multiqc multiqc --version, docker run --rm combinelab/salmon salmon --version
