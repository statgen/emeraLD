#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <unordered_set>

#ifndef PROCGENO_H
#define PROCGENO_H

#include "tabix_util/tabix.hpp"

using namespace std;

class snpinfo
{
	public:
	vector<string> chr;
	vector<int> pos;
	vector<string> rsid;
	vector<string> ref;
	vector<string> alt;
	unsigned int size() ;
	void push (string, int, string, string, string) ;
};

class idata
{
	public:
	unordered_set<string> ids;
	vector<bool> kcols;
	bool filter_mode;
	bool keep_mode;
	bool keep(string);
	bool keep(int);
	bool process(string);
	void open(string, bool);
};

class targetinfo
{
	public:
	string chr;
	int pos;
	string rsid;
	string chrpos;
	string ref;
	string alt;
	int index;
	int matches;
targetinfo() : chr(""), pos(-1), chrpos(""), rsid(""), ref(""), alt(""), index(-1), matches(0) {}
};

class hdata
{
	public:
	vector<vector<int>> iids;
	vector<vector<int>> hcts;
	vector<int> nhap;
	vector<vector<int>> map;
	void push_block (vector<int>, vector<int>, int) ;
	void push_map (vector<int>) ;
};

struct gprob {
	double p0;
	double p1;
	double p2;
};

class gdata
{
	public:
//	vector<vector<int>> carriers;
	vector<vector<int>> subsample;
	vector<vector<bool>> genotypes;
	vector<vector<uint_fast8_t>> ugenos;
//	vector<vector<int>> hets;
	vector<vector<int>> hets_ss;
//	vector<vector<int>> homs;
	vector<vector<int>> homs_ss;
	vector<gprob> p_gts;
	vector<int> dir;
	vector<int> mac;
	vector<int> block;
	void push (vector<int>, vector<bool>, int, int, int);
	void push (vector<uint_fast8_t>,vector<int>,vector<int>,double,double,double,int,int,int);
};

class foptions
{
	public:
	string region;
	string region_e;
	int min_mac;
	int max_mac;
	int max_sample;
	int mmac;
	int no_snps;
	int no_indels;
	int region_mode;
	int one_vs_all;
	int phased;
	int force_phased;
	double qual;
};

void setOptions(foptions);

vector<string> getRegion (string str);

string asRegion (int chr, int pos, int end);
string asRegion (string chr, int pos, int end);
string asRegion (string chr, string pos, string end);

int read_tabixed_vcf(string &vcf_path, targetinfo &target, gdata &gdat, snpinfo &sinfo, idata &idat, int &n_haps, int &ph);

int read_tabixed_m3vcf(string &m3vcf_path, targetinfo &target, gdata &gdat, snpinfo &sinfo, idata &idat, hdata &hdat, int &n_haps);

#endif /* PROCGENO_H */
