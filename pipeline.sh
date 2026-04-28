#!/bin/bash
source ~/.bashrc
conda activate rna_seq_env

REF_GENOME="[path_to_your_genome]/hg19.fa"
REF_ANNOTATION="[path_to_your_annotation]/hg19.gtf"

### Downloading and preparing data
echo "Downloading data..."
files=("SRR3194428" "SRR3194429" "SRR3194430" "SRR3194431" "SRR3191542" "SRR3191543" "SRR3191544" "SRR3191545")

for id in "${files[@]}"; do
    if [ ! -f "$id.fastq" ] && [ ! -f "${id}_1.fastq" ]; then
        fastq-dump $id --split-files
    fi
done

# Remove the _1 suffix for single-end reads to unify naming
for id in SRR3194428 SRR3194429 SRR3194430 SRR3194431; do
    if [ -f "${id}_1.fastq" ]; then
        mv "${id}_1.fastq" "${id}.fastq"
    fi
done

### FastQC Quality Control
echo "Running FastQC..."
fastqc *.fastq
mkdir -p qc_fastqc
mv *.html *.zip qc_fastqc/

### Trimming with Trimmomatic
echo "Running Trimmomatic..."

single_end=("SRR3194428" "SRR3194429" "SRR3194430" "SRR3194431")
paired_end=("SRR3191542" "SRR3191543" "SRR3191544" "SRR3191545")

for file in "${single_end[@]}"; do
    trimmomatic SE -threads 8 "$file.fastq" "$file.fq" SLIDINGWINDOW:4:30 TRAILING:30 MINLEN:50
done

for file in "${paired_end[@]}"; do
    r1="${file}_1"
    r2="${file}_2"
    trimmomatic PE -threads 8 "$r1.fastq" "$r2.fastq" "$r1.fq" "${r1}_unpaired.fq" "$r2.fq" "${r2}_unpaired.fq" SLIDINGWINDOW:4:30 TRAILING:30 MINLEN:50
done

# Move raw reads to a separate directory to keep the workspace clean
mkdir -p raw_fastq
mv *.fastq raw_fastq/ 2>/dev/null

### Indexing and Mapping (HISAT2)
echo "Running HISAT2..."

# Build index if it doesn't exist
if [ ! -f "hg19.1.ht2" ]; then
    hisat2-build -f -p 8 $REF_GENOME hg19
fi

for file in "${single_end[@]}"; do
    echo "Processing $file..."
    hisat2 -p 8 -x hg19 -U "$file.fq" -S "$file.single.sam"
done

for file in "${paired_end[@]}"; do
    r1="${file}_1"
    r2="${file}_2"
    echo "Processing $file..."
    hisat2 -q -p 8 -x hg19 -1 "$r1.fq" -2 "$r2.fq" -S "${file}.sam"
done

### SAM to BAM conversion and sorting
echo "Converting SAM to sorted BAM..."
sam_files=("SRR3194428.single" "SRR3194429.single" "SRR3194430.single" "SRR3194431.single" "SRR3191542" "SRR3191543" "SRR3191544" "SRR3191545")

for file in "${sam_files[@]}"; do
    samtools view -Sb -@ 8 "$file.sam" > "$file.bam"
    samtools sort -@ 8 "$file.bam" -o "${file}_sorted.bam"
    samtools index "${file}_sorted.bam"
done

### Read counting - FeatureCounts
echo "Counting reads with FeatureCounts..."
featureCounts -T 8 -a $REF_ANNOTATION -o counts_SE.txt SRR31944*_sorted.bam 
featureCounts -T 8 -p -a $REF_ANNOTATION -o counts_PE.txt SRR319154*_sorted.bam 

echo "Pipeline finished successfully!"
