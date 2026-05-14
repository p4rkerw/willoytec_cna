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

cnvkit_docker="etal/cnvkit:latest"
workdir=$(pwd)

normal=24182D-08-01_S0_L001
tumors=(
24182D-08-02_S0_L001
24182D-13-01_S109_L004
24182D-13-02_S110_L004
)

gunzip -c /mnt/g/reference/chm13v2.0.fa.gz > /mnt/g/reference/chm13v2.0.fa
ref_file=chm13v2.0.fa
ref_dir=/mnt/g/reference
mkdir -p cnvkit_output

# Build reference from normal
docker run --rm \
  -v "$workdir":/data \
  -v "$ref_dir":/ref \
  -w /data \
  $cnvkit_docker \
  cnvkit.py batch ${normal}.dedup.bam \
    --normal \
    --method wgs \
    --fasta /ref/$ref_file \
    --output-reference cnvkit_output/reference.cnn \
    --processes 8

for sample in "${tumors[@]}"; do
  echo "CNVkit processing $sample"

  docker run --rm \
    -v "$workdir":/data \
    -v "$ref_dir":/ref \
    -w /data \
    $cnvkit_docker \
    cnvkit.py batch ${sample}.dedup.bam \
      --reference cnvkit_output/reference.cnn \
      --method wgs \
      --output-dir cnvkit_output \
      --processes 8

done

for sample in "${tumors[@]}"; do
  docker run --rm \
    -v $workdir/cnvkit_output:/cnvkit_output \
    -v "$ref_dir":/ref \
    -w /data \
    -w $workdir \
    $cnvkit_docker \
    cnvkit.py scatter \
      /cnvkit_output/${sample}.dedup.cnr \
      -s /cnvkit_output/${sample}.dedup.cns \
      -o /cnvkit_output/${sample}.dedup.scatter.pdf
done

for sample in "${tumors[@]}"; do
  docker run --rm \
    -v $workdir/cnvkit_output:/cnvkit_output \
    -v "$ref_dir":/ref \
    -w /data \
    -w $workdir \
    $cnvkit_docker \
    cnvkit.py genemetrics \
      /cnvkit_output/${sample}.dedup.cnr \
      -s /cnvkit_output/${sample}.dedup.cns \
      > /cnvkit_output/${sample}.dedup.genemetrics.txt

  docker run --rm \
    -v $workdir/cnvkit_output:/cnvkit_output \
    -v "$ref_dir":/ref \
    -w /data \
    -w $workdir \
    $cnvkit_docker \
    cnvkit.py call \
      /cnvkit_output/${sample}.dedup.cns \
      -o /cnvkit_output/${sample}.dedup.call.cns
done

for sample in "${tumors[@]}"; do
  echo "==== $sample chrY ===="
  grep -w "chrY" cnvkit_output/${sample}.dedup.cns || true
done

for sample in "${tumors[@]}"; do
  echo "Summarizing chromosome-level CNAs for $sample"

  awk '
  BEGIN{OFS="\t"}

  # -------- PASS 1: collect values --------
  NR==FNR {
    if ($1=="chromosome") next

    chr=$1
    len=$3-$2
    val=$5

    # weighted mean accumulators
    sum[chr]+=len*val
    total[chr]+=len

    # store for median (segment-level, unweighted)
    vals[chr][n[chr]++]=val

    next
  }

  # -------- PASS 2: compute + print --------
  FNR==1 {
    print "chromosome","mean_log2","median_log2"
    next
  }

  END {
    for (c in sum) {

      # sort values for median
      asort(vals[c])

      nvals=n[c]
      if (nvals==0) {
        med="NA"
      } else if (nvals%2==1) {
        med=vals[c][(nvals+1)/2]
      } else {
        med=(vals[c][nvals/2]+vals[c][nvals/2+1])/2
      }

      mean=sum[c]/total[c]

      print c, mean, med
    }
  }
  ' cnvkit_output/${sample}.dedup.cns \
  > cnvkit_output/${sample}.chromosome_summary.txt

done

for sample in "${tumors[@]}"; do
  echo "==== Large CNAs (MEAN) in $sample ===="
  awk 'NR>1 && ($2 > 0.4 || $2 < -0.4)' \
    cnvkit_output/${sample}.chromosome_summary.txt
done

for sample in "${tumors[@]}"; do
  echo "==== Large CNAs (MEDIAN) in $sample ===="
  awk 'NR>1 && ($3 > 0.4 || $3 < -0.4)' \
    cnvkit_output/${sample}.chromosome_summary.txt
done

