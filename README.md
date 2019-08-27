## emeraLD
Tools for rapid on-the-fly LD calculation
#### About
- Exploits sparsity and haplotype structure to efficiently calculate LD 
- Uses tabix indexes to support rapid querying of genomic regions
- Supports VCF (phased or unphased) and M3VCF formats
- Supports integration with Python and R

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

```

#### Software References
Libraries and resources used or adapted in emeraLD:
- Data structures: [Boost C++ libraries](https://www.boost.org/)
- Python Integration: [pybind11, pybind dev team](https://github.com/pybind/pybind11)
- Tabix and HTSLIB: [tabixpp, ekg et al.](https://github.com/ekg/tabixpp) and [htslib, samtools team](https://github.com/samtools/htslib)

#### Contributors
Special thanks to Daniel Taliun and Ryan Welch

#### Citation
- To cite emeraLD, you can use [Quick et al. (2018) *Bioinformatics*](https://doi.org/10.1093/bioinformatics/bty547).

#### Feedback and bug reports
- Feel free to contact Corbin Quick (corbinq@gmail.com) with bug reports or feedback
