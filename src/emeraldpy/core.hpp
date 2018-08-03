#include <string>
#include "processGenotypes.hpp"
#include "calcLD.hpp"

using namespace std;

foptions make_default_options();

class LDResult {
public:
  string chrom;
  int pos1;
  string rsid1;
  string ref1;
  string alt1;
  string epacts1;
  int pos2;
  string rsid2;
  string ref2;
  string alt2;
  string epacts2;
  double r;
  double rsq;
  double d;
  double dprime;
  LDResult()
    : chrom(""),
      pos1(-1),
      rsid1(""),
      ref1(""),
      alt1(""),
      epacts1(""),
      pos2(-1),
      rsid2(""),
      ref2(""),
      alt2(""),
      epacts2(""),
      r(NAN),
      rsq(NAN),
      d(NAN),
      dprime(NAN)
  { }
};

vector<LDResult> ld_refsnp(string infile, string variant, string region);
