# Haplotype Caller Pipeline for GVCF Generation

**Objective**:  
This script is designed to call variants from BAM files using the GATK (Genome Analysis Toolkit) HaplotypeCaller. The output is a GVCF (Genomic Variant Call Format) file, which is essential for variant discovery and subsequent analysis, especially in population genomics and clinical genomics applications.

```bash
srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=40 -XX:ConcGCThreads=40" HaplotypeCaller \
    -I "\$input_file" \
    -R "\$GENOME" \
    -ERC GVCF \
    -O "\$OUTPUT" \
    --tmp-dir "\$temp_output"
```

This command runs the HaplotypeCaller within a Singularity container to generate a GVCF file from the input BAM file, utilizing a reference genome for accurate variant calling.

### Key Steps:

1. **Setup**:
   - Define the necessary paths for the output files, reference genome, and temporary directories. This organization ensures that all files are properly managed and easily accessible.

2. **Finding BAM Files**:
   - The script uses `find` to locate all processed BAM files that have been sorted, marked for duplicates, and recalibrated (BQSR).

3. **Creating Temporary Directories**:
   - For each sample, a temporary output directory is created to store intermediate results during processing, helping keep the workspace organized.

4. **File Integrity Check**:
   - The input BAM file undergoes an integrity check using `samtools quickcheck`. If the file is intact, the script proceeds; otherwise, an error message is logged.

5. **Read Group (RG) Information Check**:
   - The script checks the read group information in the BAM file's header to ensure that the sample name (SM) is correctly specified, which is crucial for accurate variant calling.

6. **Indexing**:
   - The script checks if the BAM file index exists. If it does not, it creates the index using `samtools index`, which is essential for random access to the BAM file during processing.

7. **Running HaplotypeCaller**:
   - The `gatk HaplotypeCaller` command is executed to perform variant calling, generating a GVCF file. Key parameters include:
     - `-I`: Input BAM file.
     - `-R`: Reference genome.
     - `-ERC GVCF`: Specifies that the output should be in GVCF format.
     - `-O`: Output file for the GVCF.
     - `--tmp-dir`: Temporary directory for intermediate files.

8. **Final Output Management**:
   - If the HaplotypeCaller completes successfully, the output GVCF file and its associated index file (.tbi) are moved to the specified output directory for further analysis. If there’s an error during the calling process, an error message is logged.

### Output:
- **GVCF File**: The GVCF file contains detailed information about variant calls, along with their associated quality scores and other relevant metrics.
- **Indexed GVCF**: An index file for the GVCF, allowing for efficient access and visualization of the variant data.

### Script Structure:
The script utilizes SLURM job arrays for parallel processing of multiple BAM files, maximizing computational resource utilization. The use of Singularity containers ensures a consistent environment for running GATK, reducing potential discrepancies due to varying software environments.

This pipeline is essential for generating high-quality variant calls from genomic data, facilitating subsequent analyses, including variant filtering, annotation, and integration into larger genomic studies.

[← download script](./scripts/08_HaplotypeCaller.sh)

```bash

#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/hcm_project/xutl
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/HCM_project/Xutl
script_name=HaplotypeCaller
script_number=08
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

file_path=\$(find $output_path/*/bam -type f -name "*_sorted.dup.bqsr.bam")

echo "all file paths: \$file_path"

IFS=\$'\r\n' GLOBIGNORE='*' command eval 'INPUTS=(\$(echo "\$file_path"))'


INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"

sample=\$(basename "\$INPUT" _sorted.dup.bqsr.bam)



# Construct the full path to the input BAM file
input_file=\$INPUT

# Create a temporary output directory specific to this job
temp_output="$temp_path/\${sample}_haplotype_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"

# Define the path to the reference genome
GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa"

# Define the path for the output GVCF file
OUTPUT="\$temp_output/\$(basename "\$input_file" _sorted.dup.bqsr.bam).g.vcf.gz"

# Specify Docker container images for GATK and Samtools
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk_latest.sif"
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk4.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"

# Print important variables for debugging
echo "Input file: \$input_file"
echo "Sample: \$sample"
echo "Output: \$OUTPUT"

# Check the integrity of the input BAM file
srun apptainer exec "\$SAMTOOLS_IMAGE" samtools quickcheck "\$input_file"

if [ \$? -eq 0 ]; then
    echo "The BAM file \$input_file is OK."
else
    echo "The BAM file \$input_file failed the quick check and may be corrupted."
fi

#check if @RG SM is correct now since i needed to correct before in ReagGroup_change >_<
srun apptainer exec "\$SAMTOOLS_IMAGE" samtools view -H \$input_file | grep '@RG'

# check if index of bam exists 
if [ ! -f "\$input_file.bai" ]; then
    srun apptainer exec "\$SAMTOOLS_IMAGE" samtools index -@ 40 "\$input_file"
else
    echo "\$input_file.bai already exists."
fi

# Run HaplotypeCaller within the GATK Docker container
srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=40 -XX:ConcGCThreads=40" HaplotypeCaller \
    -I "\$input_file" \
    -R "\$GENOME" \
    -ERC GVCF \
    -O "\$OUTPUT" \
    --tmp-dir "\$temp_output"

# Create the final output directory and move the results
if [ \$? -eq 0 ]; then
    echo "HaplotypeCaller completed successfully."

    # Create the final output directory and move the results
    output_dest_DIR=\$(dirname "\$INPUT" | sed 's/bam\$/VarCal/')
    mkdir -p "\$output_dest_DIR"
    mv "\$OUTPUT" "\$output_dest_DIR/"
    mv "\${OUTPUT}.tbi" "\$output_dest_DIR/"
else
    echo "Error: HaplotypeCaller failed."
fi






EOL

```

[7 ← Read Alignment](./07_Bqsr.md) | 8 | [9 → CombineGVCFs](./09_CombineGVCFs.md)

[← Return to main list](../README.md)