# Base Quality Score Recalibration (BQSR) Pipeline for BAM Files

**Objective**:  
This script performs Base Quality Score Recalibration (BQSR) on BAM files using the GATK (Genome Analysis Toolkit). BQSR is a critical step in preprocessing genomic data, helping to correct systematic biases in base quality scores that can arise from the sequencing process. By recalibrating these scores, the accuracy of downstream analyses, such as variant calling, is significantly improved.

```bash
srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=40 -XX:ConcGCThreads=40" BaseRecalibrator \
    -I "\$input_file" \
    -R "\$GENOME" \
    --known-sites "\$KNOWN" \
    -O "\$TABLE" \
    --tmp-dir "\$TMP"
```

This command utilizes the GATK BaseRecalibrator tool to generate a recalibration table based on known variants, which is then applied to the input BAM file to produce a corrected output BAM file.

### Key Steps:

1. **Setup**:
   - Define necessary paths for the genome reference, known variant sites, and temporary output directories. This organization is crucial for managing files and ensuring that all necessary inputs are available.

2. **Finding BAM Files**:
   - The script locates all sorted and duplicate-marked BAM files within the specified output path.

3. **Creating Temporary Directories**:
   - For each sample, a temporary output directory is created to store intermediate results during processing, minimizing clutter.

4. **File Integrity Check**:
   - The input BAM file is checked for integrity using `samtools quickcheck`. If the file passes the check, processing continues; otherwise, an error message is displayed.

5. **Base Quality Score Recalibration**:
   - The `gatk BaseRecalibrator` command is executed to generate a recalibration table, taking the input BAM file, reference genome, and known variants into account. Key parameters include:
     - `-I`: Input BAM file.
     - `-R`: Reference genome.
     - `--known-sites`: Specifies known variant sites to guide recalibration.
     - `-O`: Output file for the recalibration table.

6. **Applying BQSR**:
   - The `gatk ApplyBQSR` command is then run to apply the generated recalibration table to the input BAM file, producing a corrected output BAM file.

7. **Final Checks and Indexing**:
   - After generating the output BAM file, the script performs another integrity check using `samtools quickcheck`. If successful, it indexes the output BAM file using `samtools index` for efficient access.

8. **Organizing Outputs**:
   - The script organizes final output files and directories, moving the recalibration table and output BAM files to their respective directories. This step ensures that all results are easily accessible for subsequent analysis.

### Output:
- **Output BAM File**: The corrected BAM file is saved in the designated output directory, reflecting the recalibrated base quality scores.
- **Recalibration Table**: A table containing recalibration data is generated and stored for future reference and analysis.
- **Organized Directory Structure**: The output and recalibration files are stored in a well-structured directory to facilitate easy access.

### Script Structure:
The script is structured to handle multiple BAM files using job arrays, allowing for parallel processing and efficient resource utilization. The use of Singularity containers ensures that all necessary tools are available in a consistent environment, reducing the potential for discrepancies and errors.

This BQSR pipeline is an essential step in genomic data processing, contributing to the accuracy and reliability of subsequent analyses, such as variant calling and genomic interpretation.# Base Quality Score Recalibration (BQSR) Pipeline for BAM Files

**Objective**:  
This script performs Base Quality Score Recalibration (BQSR) on BAM files using the GATK (Genome Analysis Toolkit). BQSR is a critical step in preprocessing genomic data, helping to correct systematic biases in base quality scores that can arise from the sequencing process. By recalibrating these scores, the accuracy of downstream analyses, such as variant calling, is significantly improved.

```bash
srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=40 -XX:ConcGCThreads=40" BaseRecalibrator \
    -I "\$input_file" \
    -R "\$GENOME" \
    --known-sites "\$KNOWN" \
    -O "\$TABLE" \
    --tmp-dir "\$TMP"
```

This command utilizes the GATK BaseRecalibrator tool to generate a recalibration table based on known variants, which is then applied to the input BAM file to produce a corrected output BAM file.

### Key Steps:

1. **Setup**:
   - Define necessary paths for the genome reference, known variant sites, and temporary output directories. This organization is crucial for managing files and ensuring that all necessary inputs are available.

2. **Finding BAM Files**:
   - The script locates all sorted and duplicate-marked BAM files within the specified output path.

3. **Creating Temporary Directories**:
   - For each sample, a temporary output directory is created to store intermediate results during processing, minimizing clutter.

4. **File Integrity Check**:
   - The input BAM file is checked for integrity using `samtools quickcheck`. If the file passes the check, processing continues; otherwise, an error message is displayed.

5. **Base Quality Score Recalibration**:
   - The `gatk BaseRecalibrator` command is executed to generate a recalibration table, taking the input BAM file, reference genome, and known variants into account. Key parameters include:
     - `-I`: Input BAM file.
     - `-R`: Reference genome.
     - `--known-sites`: Specifies known variant sites to guide recalibration.
     - `-O`: Output file for the recalibration table.

6. **Applying BQSR**:
   - The `gatk ApplyBQSR` command is then run to apply the generated recalibration table to the input BAM file, producing a corrected output BAM file.

7. **Final Checks and Indexing**:
   - After generating the output BAM file, the script performs another integrity check using `samtools quickcheck`. If successful, it indexes the output BAM file using `samtools index` for efficient access.

8. **Organizing Outputs**:
   - The script organizes final output files and directories, moving the recalibration table and output BAM files to their respective directories. This step ensures that all results are easily accessible for subsequent analysis.

### Output:
- **Output BAM File**: The corrected BAM file is saved in the designated output directory, reflecting the recalibrated base quality scores.
- **Recalibration Table**: A table containing recalibration data is generated and stored for future reference and analysis.
- **Organized Directory Structure**: The output and recalibration files are stored in a well-structured directory to facilitate easy access.

### Script Structure:
The script is structured to handle multiple BAM files using job arrays, allowing for parallel processing and efficient resource utilization. The use of Singularity containers ensures that all necessary tools are available in a consistent environment, reducing the potential for discrepancies and errors.

This BQSR pipeline is an essential step in genomic data processing, contributing to the accuracy and reliability of subsequent analyses, such as variant calling and genomic interpretation.

[← download script](./scripts/07_Bqsr.sh)

```bash

#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/hcm_project/xutl
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/HCM_project/Xutl
script_name=Bqsr
script_number=07
mem=240G
cpus_per_task=40

# arrays setup
total_arrays=0



GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/GRCh38.primary.genome.fa.gz"

#create folders
logs_path=$scratch_path/output/logs
scripts_path=$scratch_path/scripts
temp_path=$scratch_path/temp

#mkdir -p $output_path
mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path





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

file_path=\$(find $output_path/*/bam -type f -name "*_sorted.dup.bam")

echo "all file paths: \$file_path"

IFS=$'\r\n' GLOBIGNORE='*' command eval 'INPUTS=(\$(echo "\$file_path"))'


INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"

sample=\$(basename "\$INPUT" _sorted.dup.bam)


input_file=\$INPUT

echo "Input file: \$input_file"
echo "Sample: \$sample"



temp_output="$temp_path/\${sample}_\${SLURM_JOB_ID}_bqsr"
mkdir -p "\$temp_output"


#GENOME AND KNOWN LOCATION
GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa"
KNOWN="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/known-sites/Homo_sapiens_assembly38.dbsnp138.vcf"

#BASE RECAL OUTPUT TABLE
TABLE=\$temp_output/\$(basename "\$input_file" .bam).recal_data.table

#TEMPORARY FILES
TMP=\$temp_output

#BQSR OUTPUT
OUTPUT=\$temp_output/\$(basename "\$input_file" .bam).bqsr.bam


echo "input: \$input_file"
echo "output table: \$TABLE"
echo "output: \$OUTPUT"



######### IMAGE #########
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk_latest.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"
########


#CHECK INPUT TO COPY
srun apptainer exec "\$SAMTOOLS_IMAGE" samtools quickcheck "\$input_file"

if [ \$? -eq 0 ]; then
    echo "The BAM file \$input_file is ok"
else
    echo "The BAM file \$input_file failed the quick check and may be corrupted."
fi


START_TIME=\$(date +%s)
echo "BASE RECAL: \$input_file"

srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=40 -XX:ConcGCThreads=40" BaseRecalibrator \
    -I "\$input_file" \
    -R "\$GENOME" \
    --known-sites "\$KNOWN" \
    -O "\$TABLE" \
    --tmp-dir "\$TMP"

echo "BQSR: \$input_file"
echo "bqsr table: \$TABLE"
srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=40 -XX:ConcGCThreads=40" ApplyBQSR \
    -I "\$input_file" \
    -R "\$GENOME" \
    --bqsr-recal-file "\$TABLE" \
    -O "\$OUTPUT" \
    --tmp-dir "\$TMP"


srun apptainer exec "\$SAMTOOLS_IMAGE" samtools quickcheck "\$OUTPUT"

srun apptainer exec "\$SAMTOOLS_IMAGE" samtools index -@ 40 "\$OUTPUT"


output_bam_DIR=\$(dirname "\$INPUT")

bqsr_table_DIR=\$(dirname "\$INPUT" | sed 's/bam\$/bqsr_table/')
mkdir -p \$bqsr_table_DIR


echo "final paths: bam - \$output_bam_DIR table - \$bqsr_table_DIR"

if [ \$? -eq 0 ]; then
    echo "EVERYTHING IS FINE"
    mv \$TABLE \$bqsr_table_DIR
    mv \$temp_output/\$sample* \$output_bam_DIR
else
    echo "The BAM file \$OUTPUT failed the quick check and may be corrupted."
fi

EOL

```

[6 ← MarkDuplicates](./06_MarkDuplicates.md) | 7 | [8 → Bam Filtering](./08_HaplotypeCaller.md)

[← Return to main list](../README.md)

