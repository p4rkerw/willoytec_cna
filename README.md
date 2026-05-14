# willoytec_cna
wgs analysis for somatic cna and mca in kidney cell lines with crispr loy

bwa=/mnt/g/software/bwa-mem2-2.2.1_x64-linux/bwa-mem2
ref=/mnt/g/reference/chm13v2.0.fa.gz

$bwa mem -t 16 $ref 24182D-08-01_S0_L001_R1_001.fastq.gz 24182D-08-01_S0_L001_R2_001.fastq.gz | samtools sort -o 24182D-08-01_S0_L001.bam
$bwa mem -t 16 $ref 24182D-08-02_S0_L001_R1_001.fastq.gz 24182D-08-02_S0_L001_R2_001.fastq.gz | samtools sort -o 24182D-08-02_S0_L001.bam
$bwa mem -t 16 $ref 24182D-13-01_S109_L004_R1_001.fastq.gz 24182D-13-01_S109_L004_R2_001.fastq.gz | samtools sort -o 24182D-13-01_S109_L004.bam
$bwa mem -t 16 $ref 24182D-13-02_S110_L004_R1_001.fastq.gz 24182D-13-02_S110_L004_R2_001.fastq.gz | samtools sort -o 24182D-13-02_S110_L004.bam

set -euo pipefail

bwa=/mnt/g/software/bwa-mem2-2.2.1_x64-linux/bwa-mem2
ref=/mnt/g/reference/chm13v2.0.fa.gz
samtools=/mnt/g/software/samtools-1.17/samtools

threads=16
t_sort=8
t_markdup=8

samples=(
24182D-08-01_S0_L001
24182D-08-02_S0_L001
24182D-13-01_S109_L004
24182D-13-02_S110_L004
)

for sample in "${samples[@]}"; do
  echo "Processing $sample"

  r1=${sample}_R1_001.fastq.gz
  r2=${sample}_R2_001.fastq.gz

  [[ -f $r1 && -f $r2 ]] || { echo "Missing FASTQ for $sample"; exit 1; }

  $bwa mem -t $threads -R "@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA" $ref $r1 $r2 | \
    $samtools sort -n -l 1 -@ $t_sort -m 2G - | \
    $samtools fixmate -m - - | \
    $samtools sort -l 1 -@ $t_sort -m 2G - | \
    $samtools markdup -@ $t_markdup -s - ${sample}.dedup.bam 2> ${sample}.dup_metrics.txt

  $samtools index ${sample}.dedup.bam
  $samtools flagstat ${sample}.dedup.bam > ${sample}.flagstat.txt

done



cnvkit=cnvkit.py
threads=16

normal=24182D-08-01_S0_L001
tumors=(
24182D-13-01_S109_L004
24182D-13-02_S110_L004
24182D-08-02_S0_L001
)


for sample in "${samples[@]}"; do
$samtools idxstats ${sample}.dedup.bam | grep chrY
done

# Configure the sources where conda will find packages
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install cnvkit

# Activate environment
conda activate cnvkit

# Verify installation
cnvkit.py version

$cnvkit batch ${normal}.dedup.bam \
  --normal \
  --method wgs \
  --fasta $ref \
  --output-reference reference.cnn \
  --processes $threads


for sample in "${tumors[@]}"; do
  echo "CNV calling for $sample"

  $cnvkit batch ${sample}.dedup.bam \
    --reference reference.cnn \
    --method wgs \
    --processes $threads \
    --output-dir cnvkit_output

done

for sample in "${tumors[@]}"; do
  $cnvkit scatter cnvkit_output/${sample}.cnr \
    -s cnvkit_output/${sample}.cns \
    -o cnvkit_output/${sample}.scatter.pdf
done

for sample in "${tumors[@]}"; do
  echo "Chromosome summary for $sample"

  $cnvkit genemetrics cnvkit_output/${sample}.cnr \
    -s cnvkit_output/${sample}.cns \
    > cnvkit_output/${sample}.genemetrics.txt

done

for sample in "${tumors[@]}"; do
  echo "chrY signal for $sample"

  grep -w "chrY" cnvkit_output/${sample}.cns || true

done

for sample in "${tumors[@]}"; do
  $cnvkit call cnvkit_output/${sample}.cns \
    -o cnvkit_output/${sample}.call.cns
done
