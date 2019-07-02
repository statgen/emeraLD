#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <random>
#include <unordered_set>

#include "tabix_util/tabix.hpp"
#include "boost/dynamic_bitset.hpp"

#ifndef PROCGENO_H
#define PROCGENO_H

using namespace std;

typedef boost::dynamic_bitset<> haploVec;

//typedef vector<bool> haploVec;

class diploVec
{
	public:
	haploVec het;
	haploVec hom;
	int operator[] (const int& x) {
		if( het[x] ){
			return ( hom[x] ? 2 : 1 );
		}else{
			return 0;
		}
	}
//	diploVec() : {};
//	diploVec(int);
	void resize (int&);
//	void reserve (int&);
	void assign(int, int&);
	void push_back (int);
	unsigned int size();
	void flip();
};

class snpinfo
{
	public:
	vector<string> chr;
	vector<int> pos;
	vector<string> rsid;
	vector<string> ref;
	vector<string> alt;
	unsigned int size() ;
	void push (string&, int&, string&, string&, string&);
};

class idata
{
	public:
	unordered_set<string> ids;
	haploVec kcols;
	bool filter_mode;
	bool keep_mode;
	bool keep(string&);
	bool keep(int);
	bool process(string&);
	void open(string&, bool);
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
	string epacts;
	int index;
	int matches;
targetinfo() : chr(""), pos(-1), rsid(""), chrpos(""), ref(""), alt(""), epacts(""), index(-1), matches(0) {}
};

class hdata
{
	public:
	vector<vector<int>> iids;
	vector<vector<int>> hcts;
	vector<haploVec> dense;
	vector<int> nhap;
	vector<vector<int>> map;
//	vector<vector<int>> count_matrix;
//	int b1;
//	int b2;
//	void get_matrix(int bi, int bj);
	void push_block (vector<int>&, vector<int>&, int) ;
	void push_map (vector<int>&, haploVec&);
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
	vector<bool> has_ss;
	vector<haploVec> genotypes;
	vector<diploVec> ugenos;
//	vector<vector<int>> hets;
	vector<vector<int>> hets_ss;
//	vector<vector<int>> homs;
	vector<vector<int>> homs_ss;
	vector<gprob> p_gts;
	vector<int> dir;
	vector<int> mac;
	vector<int> block;
	void getCarriers(int);
	void push (vector<int>&, haploVec&, int, int, int);
	void push (diploVec&,vector<int>&,vector<int>&,double&,double&,double&,int,int,int);
	void print();
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
	int m_pad;
	int m_fold;
	double qual;
};

void setOptions(foptions);

void print_hv(haploVec&,int&);
void print_dv(diploVec&,int&);

vector<string> getRegion (string str);

string asRegion (int chr, int pos, int end);
string asRegion (string chr, int pos, int end);
string asRegion (string chr, string pos, string end);
targetinfo parseEpactsVariant(std::string& variant);

int read_tabixed_vcf(string &vcf_path, targetinfo &target, gdata &gdat, snpinfo &sinfo, idata &idat, int &n_haps, int &ph);

int read_tabixed_m3vcf(string &m3vcf_path, targetinfo &target, gdata &gdat, snpinfo &sinfo, idata &idat, hdata &hdat, int &n_haps);

#endif /* PROCGENO_H */
