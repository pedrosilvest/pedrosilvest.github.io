# MarkDuplicates Pipeline for BAM Files

**Objective**:  
This script is designed to mark duplicate reads in BAM files using GATK's MarkDuplicates tool. By identifying and marking duplicates, this process enhances the accuracy of downstream analyses, such as variant calling, by ensuring that duplicate reads do not skew the results.

```bash
srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=36 -XX:ConcGCThreads=36" MarkDuplicates \
  -I "\$input_file" \
  -O "\$output_file" \
  -M "\$metrics_file" \
  --TMP_DIR "\$temp_output"
```

This step ensures that duplicate reads are properly marked in the output BAM file, facilitating accurate genomic interpretation in subsequent analysis steps.

### Key Steps:

1. **Setup Paths**:
   - Paths are defined for scratch space, output directories, logging, and temporary files, ensuring a well-organized workflow.

2. **Creating Temporary Output Directory**:
   - For each sample, a temporary directory is created to store intermediate results, preventing clutter and ensuring smooth processing.

3. **Quick Check on BAM File**:
   - A preliminary check is performed on the input BAM file using `samtools quickcheck` to verify its integrity before processing.

4. **Marking Duplicates**:
   - The `gatk MarkDuplicates` command is executed to mark duplicate reads in the BAM file. Key parameters include:
     - `-I`: Specifies the input BAM file.
     - `-O`: Specifies the output BAM file, which will contain marked duplicates.
     - `-M`: Defines the metrics file to store duplicate statistics.
     - `--TMP_DIR`: Indicates the temporary directory for intermediate files.

5. **Indexing the Output BAM File**:
   - After marking duplicates, the output BAM file is indexed using `samtools index`, which allows for efficient access to specific regions in the BAM file.

6. **Moving Results**:
   - The output BAM file and metrics file are moved to their respective directories for organization and ease of access.

7. **Cleanup**:
   - Temporary files and directories are removed to free up space and maintain a clean working environment.

### Output:
- **Output BAM File**: The BAM file containing marked duplicates is saved in the specified output directory.
- **Metrics File**: A metrics file containing statistics about the duplicates is generated and stored for later review.
- **Organized Directory Structure**: The output and metrics files are stored in designated directories, keeping the project organized.

### Script Structure:
The script utilizes a job array to handle multiple BAM files in parallel, maximizing computational efficiency and resource utilization. The use of Singularity containers ensures that the required tools are consistently available in a controlled environment.

This pipeline plays a crucial role in preprocessing BAM files for subsequent analyses, contributing to the overall reliability and accuracy of genomic research outcomes.

[← download script](./scripts/06_MarkDuplicates.sh)

```bash

#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/hcm_project/xutl
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/HCM_project/Xutl
script_name=MarkDuplicates
script_number=06
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

file_path=\$(find $output_path/*/bam -type f -name "*_sorted.bam")

echo "all file paths: \$file_path"

IFS=$'\r\n' GLOBIGNORE='*' command eval 'INPUTS=(\$(echo "\$file_path"))'


INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"


sample=\$(basename "\$INPUT" _sorted.bam)


input_file=\$INPUT
echo "Input file: \$input_file"
echo "Sample: \$sample"



# Define the output directory
temp_output="$temp_path/\${sample}_markdup_\${SLURM_JOB_ID}"

# Create the output directory if it doesn't exist
mkdir -p "\$temp_output"


################## output_file name is .dup.bam ############### 

output_file=\$temp_output/\$(basename "\$input_file" .bam).dup.bam
metrics_file=\$temp_output/\$(basename "\$input_file" .bam)_metrics_dup.txt

echo "Output: \$output_file"
echo "Metrics: \$metrics_file"

######### IMAGE #########
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk_latest.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"
########


srun apptainer exec "\$SAMTOOLS_IMAGE" samtools quickcheck "\$input_file"

if [ \$? -eq 0 ]; then
    echo "The BAM file \$input_file is ok"
else
    echo "The BAM file \$input_file failed the quick check and may be corrupted."
fi


START_TIME=\$(date +%s)

srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=36 -XX:ConcGCThreads=36" MarkDuplicates \
  -I "\$input_file" \
  -O "\$output_file" \
  -M "\$metrics_file" \
  --TMP_DIR "\$temp_output"


END_TIME=\$(date +%s)
ELAPSED_TIME=\$((END_TIME - START_TIME))

echo "Job finished in \$((\$ELAPSED_TIME / 3600)) hours, \$(((\$ELAPSED_TIME % 3600) / 60)) minutes, and \$((\$ELAPSED_TIME % 60)) seconds."


output_file_DIR=\$(dirname "\$INPUT")

metrics_file_DIR=\$(dirname "\$INPUT" | sed 's/bam\$/markdup_metrics/')
mkdir -p \$metrics_file_DIR


# Indexing
echo "INDEXING"
srun apptainer exec "\$SAMTOOLS_IMAGE" samtools index -@ 30 "\$output_file"
    
if [ \$? -eq 0 ]; then
    mv \$output_file* \$output_file_DIR
    mv \$metrics_file \$metrics_file_DIR
    rm -rf \$temp_output
else
    echo "The BAM file \$input_file failed the quick check and may be corrupted."
fi




EOL

```

[5 ← Bam Filtering](./05_FilterBam.md) | 6 | [7 → BQSR](./07_Bqsr.md)

[← Return to main list](../README.md)