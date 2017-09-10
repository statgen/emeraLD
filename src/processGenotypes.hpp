#include "tabix_util/tabix.hpp"
#include <iostream>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <stdio.h>

using namespace std;

vector<int> getRegion (string str);

string asRegion (int chr, int pos, int end);

int read_tabixed_vcf(string &vcf_path, string &region, int &region_mode, int &r_start, int &r_end, int &one_vs_all, string &ref_rsid, int &ref_snp, string &ref_ref,  string &ref_alt, int &j_ref, int &matches, vector<vector<int>> &carriers, vector<vector<bool>> &genotypes, vector<int> &dirs, vector<int> &chrs, vector<int> &poss, vector<string> &rsids, vector<string> &refs, vector<string> &alts, int &k, int &n_haps, vector<int> &BLOCK_ID);

int read_tabixed_m3vcf(string &m3vcf_path, string &region, int &region_mode, int &r_start, int &r_end, int &one_vs_all, string &ref_rsid, int &ref_snp, string &ref_ref,  string &ref_alt, int &j_ref, int &matches, vector<vector<int>> &carriers, vector<vector<bool>> &genotypes, vector<int> &dirs, vector<int> &chrs, vector<int> &poss, vector<string> &rsids, vector<string> &refs, vector<string> &alts, int &k, int &n_haps, vector<vector<int>> &HAP_IID, vector<vector<int>> &HAP_CTS, vector<int> &N_HAP, vector<int> &BLOCK_ID, vector<vector<int>> &HAP_MAP, vector<int> &MACS);

