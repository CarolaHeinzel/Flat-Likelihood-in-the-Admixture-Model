#!/bin/bash
base_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
for chr in {1..22} X Y; do
  file="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  wget "${base_url}/${file}"
done
