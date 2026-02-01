# My Personal Genomics

## vg2vcf

A robust, Unvalidated, Bash pipeline for Whole Genome Sequencing (WGS) analysis using the **HPRC Pangenome Graph** (v1.1).

This pipeline performs read mapping against a pangenome graph (VG Giraffe) and projects the data onto **two linear references simultaneously** (GRCh38 and CHM13). It outputs fully compliant BAMs, VCFs (SNPs/Indels), Structural Variants (SVs), and Copy Number Variants (CNVs).

## 🚀 Key Features
* **Graph-Based Mapping:** Uses `vg giraffe` to map against the complex HPRC graph.
* **Dual-Reference Output:** Automatically surjects reads to both **GRCh38** and **CHM13** (T2T) coordinate systems.
* **State-of-the-Art Callers:**
    * **Small Variants:** Google DeepVariant (WGS model).
    * **Structural Variants:** `vg call` (Graph-based SV detection).
    * **CNVs:** CNVkit (Moving average approach).
* **Resumable:** Tracks progress in a `.txt` file; automatically skips completed steps.
* **Dockerized:** All tools run in isolated containers (VG, DeepVariant, BCFTools, CNVkit, MultiQC).

    ```
