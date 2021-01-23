#!/bin/sh

# Example:
# Identify allelic imbalance (AI) given a tumor 
# exome VCF file generated using the UnifiedGenotyper 
# of the GATK.  This involves the following 3 steps.

printf "STEP 1: PHASING 1KG HET SITES ...\n"
mkdir example_output
python ../scripts/simple_phaser.py \
  --ldmap ../ldmap/hg19.exome.ldmap \
  --vcf example_input/tumor_exome.vcf \
  -o example_output/tumor_exome

printf "\nSTEP 2: IDENTIFYING REGIONS OF AI ...\n"
../haplohseq \
  --vcf example_output/tumor_exome.hap.vcf \
  --phased example_output/tumor_exome.hap \
  --event_prevalence 0.1 \
  -d example_output \
  -p tumor_exome_haplohseq

printf "\nSTEP 3: PLOTTING HAPLOHSEQ GENOMIC AI PROFILE ...\n"
Rscript ../scripts/haplohseq_plot.R \
  --file example_output/tumor_exome_haplohseq.posterior.dat \
  --out example_output \
  --prefix tumor_exome_haplohseq
