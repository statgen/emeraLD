#include <string>
#include "processGenotypes.hpp"
#include "calcLD.hpp"
#include "emeraldpy/core.hpp"

using namespace std;

double MIN_PRINT = 0.00001;
int MAX_DIST = 1000000;

foptions make_default_options() {
  foptions fopts;
  fopts.max_mac = 100000000;
  fopts.min_mac = 1;
  fopts.m_fold = 0;
  fopts.m_pad = 50;
  fopts.mmac = 1000;
  fopts.max_sample = 1000;
  fopts.one_vs_all = 0;
  fopts.force_phased = 0;
  fopts.phased = 1;
  fopts.region_mode = 0;
  fopts.region = "";
  fopts.region_e = "";
  return fopts;
}

// I'm assuming return value optimization gets rid of the copy on returning this vector,
// but would be interesting to try and verify that sometime.
vector<LDResult> ld_refsnp(string infile, string variant, string region) {
  static vector<LDResult> empty;

  // Setup necessary objects
  foptions fopts = make_default_options();
  gdata gdat;
  hdata hdat;
  idata idat;
  snpinfo sinfo;

  // Setup default fopts
  fopts.region_mode = 1;
  fopts.region = region;
  fopts.one_vs_all = 1;

  // Check region
  vector<string> region_v = getRegion(region);
  if (region_v.size() < 3) {
    cerr << "ERROR: region input format must follow --region chr:start-end\n\n";
    return empty;
  }

  // Parse the EPACTS formatted variant (chr:pos_ref/alt)
  targetinfo target = parseEpactsVariant(variant);

  // Check infile is accessible
  if (access( infile.c_str(), F_OK ) == -1 ) {
    cerr << "ERROR: input file '" << infile << "' does not exist...\n\n";
    return empty;
  }

  // Check for tabix index
  string infile_tbi = infile + ".tbi";
  if (access( infile_tbi.c_str(), F_OK ) == -1 ) {
    cerr << "ERROR: index file '" << infile_tbi << "' does not exist...\n";
    cerr << "\temeraLD requires bgzipped & tabixed input files";
    cerr << "\tuse 'tabix -p vcf " << infile << "' to generate index\n\n";
    return empty;
  }

  // Are we reading from a m3vcf?
  int m3vcf = 0;
  if( infile.find("m3vcf") != string::npos ) {
    m3vcf = 1;
  }

  setThresh(MIN_PRINT);
  setMaxSample(fopts.max_sample);
  setPhased(fopts.phased);
  setOptions(fopts);

  // Scan VCF and collect information on samples?
  idat.filter_mode = false;
  idat.process(infile);

  /**
   * Initial scan of m3vcf or VCF. Accumulates variant information in sinfo.
   * Determines whether the reference (target) variant exists in the VCF, and fills in
   * information about it from the VCF if necessary, such as:
   *    - REF and ALT alleles, if only rsID was specified
   * Accumulate haplotypes and their count
   */
  int n_haps;
  if (!m3vcf) {
    if (!fopts.phased) {
      cerr << "Assuming unphased data (reporting diploid genotype LD)...\n";
    }
    if (read_tabixed_vcf(infile, target, gdat, sinfo, idat, n_haps, fopts.phased) > 0 ) {
      cerr << "ERROR: check vcf file " << infile << "\n";
      return empty;
    }
  } else{
    if (!fopts.phased) {
      cerr << "ERROR: --nophase (genotype LD) currently only supported for VCF format\n";
      return empty;
    }
    if (read_tabixed_m3vcf(infile, target, gdat, sinfo, idat, hdat, n_haps) > 0) {
      cerr << "ERROR: check m3vcf file " << infile << "\n";
      return empty;
    }
  }

  setSize(n_haps);

  if (fopts.one_vs_all) {
    if (target.matches < 1) {
      cerr << "ERROR: variant " << target.epacts << " not found!\n";
      return empty;
    }
  }

  if (gdat.mac.size() < 1) {
    if (fopts.region_mode) {
      cerr << "ERROR: Found no SNPs in region " << fopts.region << " ... \n";
      cerr << "hint: make sure chr:pos formatting matches input file\n\n";
    } else {
      cerr << "ERROR: For no SNPs in input file!\n";
      cerr << "hint: check input file formatting\n\n";
    }
    return empty;
  }

  double r, d, dprime;

//  fprintf(stdout, "#CHR\tPOS1\tRSID1\tREF:ALT1\tPOS2\tRSID2\tREF:ALT2\t");
//  fprintf(stdout, "R\tRsq\tD\tDprime\n");

  vector<LDResult> results;
  for (int i = 0; i < sinfo.size(); i++) {
    bool not_target = sinfo.pos[i] != target.pos && sinfo.ref[i] != target.ref && sinfo.alt[i] != target.alt;
    if (abs(target.pos - sinfo.pos[i]) < MAX_DIST && not_target) {
      getCorr(r, d, dprime, i, target.index, gdat, hdat);

      LDResult result;
      result.chrom = sinfo.chr[i];
      result.pos1 = target.pos;
      result.rsid1 = target.rsid;
      result.ref1 = target.ref;
      result.alt1 = target.alt;
      result.epacts1 = sinfo.chr[i] + ":" + to_string(target.pos) + "_" + target.ref + "/" + target.alt;
      result.pos2 = sinfo.pos[i];
      result.rsid2 = sinfo.rsid[i];
      result.ref2 = sinfo.ref[i];
      result.alt2 = sinfo.alt[i];
      result.epacts2 = sinfo.chr[i] + ":" + to_string(sinfo.pos[i]) + "_" + sinfo.ref[i] + "/" + sinfo.alt[i];
      result.r = r;
      result.rsq = pow(r,2);
      result.d = d;
      result.dprime = dprime;
      results.push_back(result);

//      if (abs(r) > MIN_PRINT) {
//        fprintf(stdout,"%s\t", sinfo.chr[i].c_str());
//        fprintf(stdout,"%u\t", target.pos);
//        string outl = target.rsid + "\t" + target.ref + ":" + target.alt + "\t";
//        fprintf(stdout,"%s", outl.c_str());
//        fprintf(stdout,"%u\t", sinfo.pos[i]);
//        outl = sinfo.rsid[i] + "\t" + sinfo.ref[i] + ":" + sinfo.alt[i] + "\t";
//        fprintf(stdout,"%s", outl.c_str());
//        fprintf(stdout,"%.5f\t%.5f\t%.5f\t%.5f\n", r, pow(r,2), d, dprime );
//      }
    }
  }

  return results;
}

int main (int argc, char *argv[]) {
  vector<LDResult> results = ld_refsnp("tests/data/test_vcf.vcf.gz", "20:69094_G/A", "20:61098-80071");
  return 0;
}
