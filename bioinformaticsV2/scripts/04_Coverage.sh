#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project
script_name=Coverage
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
# Should be a list of paths!!
file_paths=(
"/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/CBE7negB"
"/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/CBE7ttdB"
)

# Get the total number of files (array length)
total_files=${#file_paths[@]}

# Subtract 1 from the total number of files
total_arrays=$((total_files - 1))


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



file_paths=($(printf '"%s"\n' "${file_paths[@]}"))

INPUT="\${file_paths[\$SLURM_ARRAY_TASK_ID]}"
echo "\$INPUT"

sample=\$(basename "\$INPUT")

bam_file=(\$INPUT/bam/\${sample}.bam)


SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"
COVERAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_coverage_tools_latest.sif"

#------
echo "Starting: Mosdepth coverage: \$bam_file"

temp_output="$temp_path/\${sample}_\${SLURM_JOB_ID}_mosdepth"
mkdir -p "\$temp_output"

gene_panel=\${temp_output}/genes.bed
GTF_GENCODE="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/gencode/gencode.v44.chr_patch_hapl_scaff.annotation.gtf.gz"

zgrep -E 'HAVANA.*\bgene\b' "\$GTF_GENCODE" | awk -v OFS='\t' '{ gsub(/[";]/, "", \$10); print \$1, \$4, \$5, \$10 }' > "\$gene_panel"

srun apptainer exec "\$COVERAGE" mosdepth -t 4 --by \$gene_panel \$temp_output/\${sample} \$bam_file

mv \$temp_output \$INPUT/stats

#-----
echo "Starting Samtools depth: \$bam_file"


temp_output="$temp_path/\${sample}_\${SLURM_JOB_ID}_samdepth"
mkdir -p "\$temp_output"


# Define output file paths
output_depth="\$temp_output/\${sample}_samdepth.txt"
output_coverage="\$temp_output/\${sample}_samcoverage.txt"

echo "samtools depth to get coverage"
srun apptainer exec "\$SAMTOOLS_IMAGE" samtools depth -a "\$bam_file" > "\$output_depth"

# -a to include 0 coverage positions
echo "calculating average coverage"
awk '{sum+=\$3} END {print "Average coverage: " sum/NR}' "\$output_depth" > "\$output_coverage"

# Correct the destination paths for 'mv' commands
mv "\$temp_output" "\$INPUT/stats"


EOL