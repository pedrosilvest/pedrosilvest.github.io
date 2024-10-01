# To DO:

## Coverage Analysis Pipeline for BAM Files

This pipeline is designed to compute the depth and coverage of sequencing data using Mosdepth and Samtools. The coverage analysis focuses on specific gene regions, helping to evaluate how well the sequencing data covers important genomic regions, such as gene panels of interest.

Once the sequencing reads have been aligned to the reference genome and stored in BAM files, the next possible step is to assess how well the data covers the genomic regions of interest. This information is critical for ensuring that the sequencing depth is sufficient for accurate variant detection or other analyses. The pipeline calculates both gene-specific and overall coverage using Mosdepth and Samtools, ensuring comprehensive coverage metrics for your dataset.

In the script below, BAM files are processed to compute coverage over gene panels and across the entire genome using two tools:

1. Mosdepth: Provides detailed coverage information for specific regions of interest, such as gene panels, at a faster speed with reduced memory requirements.
2. Samtools depth: Offers a straightforward calculation of depth at each position across the genome, including regions with zero coverage, to compute the average coverage.
Here are the key parameters used:

Mosdepth is run with the --by flag to restrict coverage calculations to a gene panel, extracted from a GTF file, ensuring focus on the relevant genomic regions.
Samtools depth computes the depth at each genomic position and is followed by an awk command to calculate the average coverage across all positions.
bash
Copy code
# Extract paths for BAM files
bam_file=($INPUT/bam/${sample}.bam)

# Run Mosdepth for gene-specific coverage
srun apptainer exec "$COVERAGE" mosdepth -t 4 --by $gene_panel $temp_output/${sample} $bam_file

# Run Samtools depth to calculate coverage at each position
srun apptainer exec "$SAMTOOLS_IMAGE" samtools depth -a "$bam_file" > "$output_depth"

# Calculate average coverage
awk '{sum+=$3} END {print "Average coverage: " sum/NR}' "$output_depth" > "$output_coverage"
In this process:

Gene Panel Extraction: A BED file is generated from a GTF file to focus coverage calculations on specific genes of interest.
Coverage Calculation: Mosdepth and Samtools depth provide both gene-level and genome-wide coverage statistics.
Output: Results are saved in a stats directory for each sample, including depth and average coverage metrics.
By performing coverage analysis at this stage, you can ensure that your data has sufficient depth across key regions, which is essential for the reliability of downstream variant discovery or other genomic analyses.






