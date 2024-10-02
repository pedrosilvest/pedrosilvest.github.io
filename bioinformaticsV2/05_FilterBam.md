# BAM File Filtering, Sorting, Indexing, and Statistics Pipeline

**Objective**:  
This script processes BAM files by filtering low-quality reads, sorting the resulting BAM files, indexing them for efficient access, and generating statistics. This ensures that only high-quality data is retained for downstream analyses, enhancing the reliability of genomic studies.

```bash
srun apptainer exec "\$SAMTOOLS_IMAGE" samtools view -@ 30 -q 10 -F 256 -F 4 -b "\$INPUT" > "\$filtered_bam"

srun apptainer exec "\$SAMTOOLS_IMAGE" samtools sort -@ 30 "\$filtered_bam" -o "\$sorted_file"
```

By employing these steps, the pipeline ensures that only high-quality, properly organized BAM files are available for downstream analyses, making it essential for accurate genomic interpretation.

### Key Steps:
1. **Filtering**: 
   - The script filters out low-quality reads from the BAM file using the `samtools view` command with flags:
     - `-q 10`: Discards reads with a mapping quality lower than 10.
     - `-F 256`: Excludes secondary alignments.
     - `-F 4`: Excludes unmapped reads.
   - This step ensures that only reliable reads are kept for further processing.

2. **Quick Check**:
   - A quick check on the filtered BAM file is performed to ensure it is not corrupted, which is crucial for maintaining data integrity before proceeding to sorting.

3. **Sorting**:
   - The filtered BAM file is sorted using `samtools sort`, which organizes the reads in a specific order, typically by their genomic coordinates. This organization is vital for efficient data retrieval and downstream analysis.

4. **Indexing**:
   - The sorted BAM file is indexed using `samtools index`, creating an index file that allows for rapid access to specific regions of the BAM file. This is essential for many bioinformatics tools that require quick lookups.

5. **Statistics Generation**:
   - The script generates alignment statistics using `samtools flagstat`, which provides a summary of read alignments and is useful for assessing the quality of the sequencing data.

### Output:
- **Sorted BAM File**: The processed BAM file is saved in the specified output directory, ready for downstream analyses.
- **Index File**: An index file corresponding to the sorted BAM file is also saved, allowing for efficient access.
- **Statistics Report**: A statistics report summarizing the alignment quality is generated and stored, aiding in the evaluation of the sequencing experiment.

### Directory Structure:
- The sorted BAM files and their indices are moved to a `bam/` directory.
- Alignment statistics are moved to a `stats/` directory, ensuring organized output that is easy to navigate.

This pipeline effectively filters, organizes, and summarizes BAM files, providing a robust foundation for subsequent analyses in genomic research.

[← download script](./scripts/05_FilterBam.sh)

```bash

#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/hcm_project/xutl
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/HCM_project/Xutl
script_name=FilterBam
script_number=05
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

# Define the path to the Singularity container image with Samtools
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"

# Define an array of file paths for BAM files to process
file_path=\$(find $output_path/*/bam -type f -name "*.bam")


echo "all paths to do: \$file_path"

IFS=$'\r\n' GLOBIGNORE='*' command eval 'INPUTS=(\$(echo "\$file_path"))'

INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"

echo "Input file: \$INPUT"
echo "Dest dir: \$(dirname "\$INPUT")"
# Extract the sample name from the input BAM file
sample=\$(basename "\$INPUT" .bam)
echo "Sample: \$sample"

# Create a temporary output directory
temp_output="$temp_path/\${sample}_filter_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"

# Print debugging information
echo "Processing BAM file: \$INPUT"
echo "Output directory: \$temp_output"

# Filtering
echo "FILTERING"
filtered_bam="\$temp_output/\$(basename "\$INPUT")"

srun apptainer exec "\$SAMTOOLS_IMAGE" samtools view -@ 30 -q 10 -F 256 -F 4 -b "\$INPUT" > "\$filtered_bam"

# Quick check on the filtered BAM file
if srun apptainer exec "\$SAMTOOLS_IMAGE" samtools quickcheck "\$filtered_bam"; then
    echo "The BAM file \$filtered_bam passed the quick check and is not corrupted."
    
    # Sorting
    echo "SORTING"
    sorted_file="\$temp_output/\${sample}_sorted.bam"
    srun apptainer exec "\$SAMTOOLS_IMAGE" samtools sort -@ 30 "\$filtered_bam" -o "\$sorted_file"
    
    # Indexing
    echo "INDEXING"
    srun apptainer exec "\$SAMTOOLS_IMAGE" samtools index -@ 30 "\$sorted_file"
    
    # Stats
    echo "GENERATING STATS"
    STATS_OUTPUT="\$temp_output/\${sample}_sorted_stats.txt"
    srun apptainer exec "\$SAMTOOLS_IMAGE" samtools flagstat -@ 30 "\$sorted_file" > "\$STATS_OUTPUT"

fi


BAM_DIR=$output_path/\$sample/bam
STATS_DIR=$output_path/\$sample/stats

echo "MOVING \$sorted_file to \$BAM_DIR"
mv "\$sorted_file" "\$BAM_DIR"
mv "\$sorted_file.bai" "\$BAM_DIR"
mv "\$STATS_OUTPUT" "\$STATS_DIR"

rm -rf \$temp_output





EOL

```

[Previous -> Coverage](./04_Coverage.md) | 5 | [Next → Bam Filtering](./05_FilterBam.md) | [← Return to main list](../README.md)