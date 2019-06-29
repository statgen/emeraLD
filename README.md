## emeraLD
Tools for rapid on-the-fly LD calculation
#### About
- Exploits sparsity and haplotype structure to efficiently calculate LD 
- Uses tabix indexes to support rapid querying of genomic regions
- Supports VCF (phased or unphased) and M3VCF formats
- Provides easy-to-use R interface to avoid storing/precomputing LD. This can save space without compromising speed for GWAS analysis in R 
- In data sets with 10Ks of samples, emeraLD is 5-20x faster than PLINK-1.9 and 100s-Ks of times faster than existing tools for VCF files
#### Installing 
```bash
git clone https://github.com/statgen/emeraLD.git  
cd emeraLD  
make  
```
#### Usage 
- Example usage from command line  
```bash
# example usage for calculating LD in a region:
bin/emeraLD -i example/chr20.1KG.25K_m.m3vcf.gz --region 20:60479-438197 --stdout | bgzip -c > output.txt.gz

- Example usage from R interface (see emeraLD2R.r file) 
```R
## "$chr" functions as chromosome wildcard 
in_path <- "example/chr$chr.1KG.25K_m.vcf.gz"

## emeraLD2R creates an LD retrieval function 
getLD <- emeraLD2R(path = in_path)

## calling an LD retrieval function invokes emeraLD
ld_data <- getLD(region = "20:83061-92955")
## emeraLD processes genotypes; LD output is passed to R

## check LD output 
head(ld_data$Sigma[, 1:10], 10)
head(ld_data$info)
```
#### Citation
- To cite emeraLD, you can use [Quick et al. (2018) *Bioinformatics*](https://doi.org/10.1093/bioinformatics/bty547).

#### Feedback and bug reports
- Feel free to contact Corbin Quick (corbinq@gmail.com) with bug reports or feedback
