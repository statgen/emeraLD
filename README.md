## emeraLD
Tools for rapid on-the-fly LD calculation
#### About
- Exploits sparsity and haplotype structure to enable rapid LD calculation with massive data sets
- Uses tabix indexes to support rapid querying of genomic regions
- Supports VCF and M3VCF formats
- Provides easy-to-use R interface to avoid storing/precomputing LD. This can save space without compromising speed for GWAS analysis in R 
- In data sets with 10Ks of samples, emeraLD is 100s-Ks of times faster than standard tools (e.g. VCFtools and RAREMETALWORKER) and 20-80 times faster than tools for custom genotype formats (e.g., m3vcftools and LDstore)
#### Installing 
```bash
git clone https://github.com/statgen/emeraLD.git  
cd emeraLD  
make cloneLib  
make  
```
Or if libStatGen is already present,  
```bash
git clone https://github.com/statgen/emeraLD.git  
cd emeraLD  
make LIBSTATGEN_PATH=my/path/to/libStatGen  
```
#### Usage 
- Example usage from command line  
```bash
# example usage for custom region:
bin/emeraLD -i example/chr20.1KG.25K_m.m3vcf.gz --region 20:60479-438197 --stdout | bgzip -c > my_LD.txt.gz

# using --dstats flag to include D and D' statistics:
bin/emeraLD -i example/chr20.1KG.25K_m.m3vcf.gz --region 20:60479-438197 --stdout --dstats | head
```
- Example usage from R interface (see emeraLD2R.r file) 
```R
## use "$chr" when genotype files are separated by chromosome
## "$chr" is automatically replaced when a region is specified
in_path <- "example/chr$chr.1KG.25K_m.m3vcf.gz"

## emeraLD2R creates an LD retrieval function 
emeraLD <- emeraLD2R(path = in_path)

## calling an LD retrieval function invokes emeraLD
ld_data <- getLD(region = "20:83061-92955")
## emeraLD processes genotypes; LD output is passed to R

## check LD output 
head(ld_data$Sigma[, 1:10], 10)
head(ld_data$info)
```
#### Feedback and bug reports
- Feel free to contact Corbin Quick (corbinq@gmail.com) with bug reports or feedback
