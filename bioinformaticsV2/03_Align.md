
# Alignment and BAM Processing Pipeline for Paired-End Sequencing Data

**Objective**:  
This script aligns paired-end sequencing data to a reference genome using BWA-MEM, sorts the resulting BAM files, and generates alignment statistics. By efficiently processing trimmed reads (`R1` and `R2`), it prepares the data for downstream analyses such as variant calling and further quality assessment.

```bash
apptainer exec "\$BWA" bwa mem -t 30 -R "\$RG" "\$GENOME" "\$R1" "\$R2" | \
  apptainer exec "\$SAMTOOLS_IMAGE" samtools view -b -o "\$BAM_OUTPUT" -

apptainer exec "\$SAMTOOLS_IMAGE" samtools sort -@ 30 \$BAM_OUTPUT -o "\$SAVE"
```

By aligning the reads to a reference genome and sorting the BAM files, this pipeline ensures that sequencing data is properly organized and ready for further analysis.

### Key Parameters:
1. **Read Group (`RG`)**: The read group information (`ID`, `SM`, `PL`) is included in the alignment to track metadata like sample name and platform. This is essential for accurate downstream processing, especially in multi-sample projects.
2. **`--threads 30` (`-t 30`)**: Uses 30 threads to speed up the alignment and BAM file sorting processes, making the pipeline more efficient on high-performance clusters.
3. **BWA and SAMtools**: The script utilizes the BWA-MEM algorithm for alignment and SAMtools for BAM file processing (conversion, sorting, and generating statistics), ensuring high performance and scalability for large datasets.

### Output:
- **`sorted.bam`**: The sorted BAM file is the main output, which is necessary for downstream analysis tools such as variant callers.
- **`flagstat` statistics**: A summary of alignment statistics, providing insights into the quality and completeness of the alignment.

### Directory Structure:
- BAM files are stored in a `bam/` directory.
- Alignment statistics are stored in a `stats/` directory.

This pipeline ensures efficient alignment and processing of sequencing data, while maintaining clean and organized output for future analyses.

[← download script](./scripts/03_Align.sh)

```bash


#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/hcm_project/xutl
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/HCM_project/Xutl
script_name=Align
script_number=03
mem=240G
cpus_per_task=40



GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/GRCh38.primary.genome.fa.gz"

#create folders
logs_path=$scratch_path/output/logs
scripts_path=$scratch_path/scripts
temp_path=$scratch_path/temp

#mkdir -p $output_path
mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path


# arrays setup
total_arrays=$(( $(echo "$fasta_paths" | wc -l) - 1 ))



setup_sbatch="$scripts_path/${script_number}_${script_name}.sbatch"

echo "Creating scripts... $setup_sbatch"
cat >"$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=$script_name
#SBATCH --output=$logs_path/${script_number}_${script_name}_%A_%a.out
#SBATCH --array=0-$total_arrays%5
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --mem=$mem
#SBATCH --cpus-per-task=$cpus_per_task

GENOME=$GENOME

file_path=\$(find $output_path/*/trim -type f -name "*R1*fq.gz")


IFS=$'\r\n' GLOBIGNORE='*' command eval 'INPUTS=(\$(echo "\$file_path"))'

INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"

# R1 and R2 files
R1="\$INPUT"
R2="\${R1//_R1/_R2}"
R2="\${R2//val_1/val_2}"

echo "R1: \$R1"
echo "R2: \$R2"

sample=\$(basename "\$INPUT" | cut -d "_" -f 1)
full_name=\$(basename "\$R1" | awk -F '_R1' '{print \$1}')

# Create a directory for the sample in a temporary folder with job ID
temp_output="$temp_path/\${sample}_align_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"



BAM_OUTPUT="\$temp_output/\$full_name.bam"  # Updated this line
SAVE="\$temp_output/\$full_name""_sorted.bam"
STATS_OUTPUT="\$temp_output/\$full_name""_stats.txt"

echo "full_name: \$full_name"
echo "folder name: \$sample"
echo "sorted_output: \$SAVE"
echo "bam output: \$BAM_OUTPUT"
echo "stats output: \$STATS_OUTPUT"

SM="\$sample"
ID="\$full_name""_id"
RG="@RG\tID:\$ID\tSM:\$SM\tPL:ILLUMINA"    # AQUI DEVIA TER TIDO CUIDADO POR TER LANES. O SM é a SAMPLE C2040neg e o ID é o unique id.

echo "\$RG"

# IMAGES
BWA="/mnt/beegfs/apptainer/images/bwa_latest.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"

# BAM MEM ALIGNMENT
#srun apptainer exec "\$BWA" bwa mem -t 30 -R "\$RG" "\$GENOME" "\$R1" "\$R2" > "\$BAM_OUTPUT"

#srun apptainer exec "\$SAMTOOLS_IMAGE" samtools sort -@ 30 "\$BAM_OUTPUT" -o "\$SAVE"

apptainer exec "\$BWA" bwa mem -t 30 -R "\$RG" "\$GENOME" "\$R1" "\$R2" | \
  apptainer exec "\$SAMTOOLS_IMAGE" samtools view -b -o "\$BAM_OUTPUT" -

apptainer exec "\$SAMTOOLS_IMAGE" samtools sort -@ 30 \$BAM_OUTPUT -o "\$SAVE"

# STATISTICS
apptainer exec "\$SAMTOOLS_IMAGE" samtools flagstat -@ 30 "\$SAVE" > "\$STATS_OUTPUT"

apptainer exec "\$SAMTOOLS_IMAGE" samtools view -H "\$SAVE" | grep "^@RG"

bam_dir="$output_path/\$sample/bam/"
stats_dir="$output_path/\$sample/stats/"
mkdir -p \$bam_dir
mkdir -p \$stats_dir

mv "\$temp_output"/*sorted.bam "\$bam_dir/\${sample}.bam"
mv "\$temp_output"/*stats.txt "\$stats_dir"

rm -rf "\$temp_output"



EOL



```

[Previous -> Trimming](./02_TrimGalore.md) | 3 | [Next → Coverage](./04_Coverage.md) | [← Return to main list](../README.md)

