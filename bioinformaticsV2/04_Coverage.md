# Coverage Analysis Pipeline for BAM Files

**Objective**:  
This script analyzes coverage for BAM files generated from sequencing data, using both Mosdepth and Samtools to assess the depth of coverage across genes. It generates coverage statistics, which are crucial for evaluating the quality and completeness of genomic data in downstream analyses such as variant calling.

```bash
srun apptainer exec "\$COVERAGE" mosdepth -t 4 --by \$gene_panel \$temp_output/\${sample} \$bam_file

srun apptainer exec "\$SAMTOOLS_IMAGE" samtools depth -a "\$bam_file" > "\$output_depth"
```

By calculating coverage metrics, this pipeline helps ensure that the sequencing data is adequately represented across the genome, which is vital for reliable variant detection and other analyses.

### Key Steps:
1. **Gene Panel Creation**: The script extracts gene information from a GTF file to create a BED file (`genes.bed`) that defines regions of interest for coverage calculation. This ensures that the analysis focuses on relevant genomic features.
2. **Mosdepth Coverage Calculation**: Mosdepth is executed to compute coverage metrics for the specified genes, providing detailed coverage information that helps assess how well the target regions are represented in the sequencing data.
3. **Samtools Depth Calculation**: Samtools is used to calculate coverage depth at every position in the BAM file, including positions with zero coverage, which is critical for a comprehensive understanding of the data.
4. **Average Coverage Calculation**: The script calculates the average coverage across all positions from the output of the Samtools depth analysis, giving a single value that summarizes the overall coverage quality.

### Output:
- **Coverage Statistics**: Results from Mosdepth and Samtools are stored in the `stats` directory within the input BAM file's path, making it easy to locate coverage information alongside the BAM files.
- **Average Coverage Report**: The average coverage is reported, providing a quick metric for evaluating the overall sequencing depth.

This coverage analysis pipeline ensures thorough evaluation of BAM files, enhancing the reliability of downstream genomic analyses.

[← download script](./scripts/04_Coverage.sh)

```bash

to edit

```

[3 ← Read Alignment](./03_Align.md) | 4 | [5 → Bam Filtering](./05_FilterBam.md)

[← Return to main list](../README.md)

