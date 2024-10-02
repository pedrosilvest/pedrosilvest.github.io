# Quality Control Pipeline for Paired-End Sequencing Data:

This script performs an initial quality control (QC) assessment of paired-end sequencing data using FastQC. It processes both R1 and R2 reads from sequencing files, ensuring the data is suitable for downstream analyses by generating detailed quality reports.

Quality control is a crucial first step in any sequencing data analysis pipeline. By identifying potential issues like low-quality reads, adapter contamination, or sequence biases early on, this pipeline helps avoid costly errors and ensures that only high-quality data is used in subsequent analysis steps. This initial assessment provides transparency and reliability for the rest of your workflow, making it essential to perform before moving forward.

In the below script, since I had both R1 and R2 files in the same path, I checked whether the current file was an R1 file. If it was, I extracted the corresponding R2 file by replacing "R1" with "R2" in the file name. I then ran FastQC on both files using 20 threads to ensure efficient processing.

[← download script](./scripts/01_QualityControl.sh)

```bash


```

[1 ← Quality Control](./01_QualityControl.md) | 2 | [3 → Read Alignment](./03_Align.md)

[← Return to main list](../README.md)