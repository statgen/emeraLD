#!/bin/sh

## example using small chr20 region from 1000G P3
## emeraLD is most useful for larger data sets with 10s of Ks of samples 

## m3vcf files save space and computation 
input="example/chr20.1KG.25K_m.m3vcf.gz"
## see m3vcftools (http://genome.sph.umich.edu/wiki/M3vcftools) to learn more

region="20:60479-238197"
prefix="example"

if [ -e bin/emeraLD ]; then
    printf  "\nexample command: bin/emeraLD -i $input --region $region --out $prefix\n\n"
    time bin/emeraLD -i $input --region $region --out $prefix
else
    printf "\n\tbin/emeraLD does not exist\n"
    printf "\n\tdid you compile emeraLD? try running 'make'\n\n"
fi


