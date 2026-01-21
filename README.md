### Drosophila Variant Pipeline

## Introduction

The primary objective of this project is to develop a reproducible, automated pipeline to identify high-confidence genomic variations in Drosophila melanogaster by processing real-world sequencing data. 

By utilizing the Snakemake workflow management system, the project aims to transform raw, short-read sequencing data into a curated list of genetic variants, specifically Single Nucleotide Polymorphisms (SNPs) and Small Indels. 

The expected result is a statistically validated set of biological markers that differentiate a specific fly line from the standard reference genome, providing actionable insights for evolutionary and functional genomics research.

## Background

* The Reference Genome: A standardized map of an organism's DNA, such as the BDGP6.46 assembly for Drosophila, which provides a baseline but does not capture the full diversity of a living population.

* Genetic Variants: Discrete changes in the DNA sequence that drive the majority of biological differences between individuals.

* Common Variant Types:
    * SNPs (Single Nucleotide Polymorphisms): Single-letter swaps within the DNA sequence, such as changing an A to a G.

    * Indels: The insertion or deletion of small segments of DNA.

* Alignment (Mapping): The process of comparing millions of short DNA fragments (Reads) against the reference genome to identify mismatches that indicate potential variants.

* Transition to Transversion (Ti/Tv) Ratio: A critical quality metric used to evaluate the biological validity of discovered variants.

* Biological Basis: In natural evolution, transitions (purine-to-purine or pyrimidine-to-pyrimidine) occur roughly twice as often as transversions.

* Validation: A Ti/Tv ratio near 2.0 suggests the variants follow natural biological patterns, whereas a significantly lower ratio typically indicates technical noise or sequencing errors.

