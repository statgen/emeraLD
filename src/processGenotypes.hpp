#include "tabix_util/tabix.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <unordered_set>

using namespace std;

class snpinfo
{
	public:
		vector<int> chr; 
		vector<int> pos;
		vector<string> rsid; 
		vector<string> ref; 
		vector<string> alt; 
		unsigned int size() ;
		void push (int, int, string, string, string) ;
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
		int chr; 
		int pos; 
		string rsid;
		string chrpos;
		string ref;  
		string alt; 
		int index; 
		int matches;
		targetinfo() : chr(-1), pos(-1), chrpos(""), rsid(""), ref(""), alt(""), index(-1), matches(0) {}
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

class gdata 
{
	public:
		vector<vector<int>> carriers; 
		vector<vector<bool>> genotypes;
		vector<int> dir;
		vector<int> mac;
		vector<int> block;
		void push (vector<int>, vector<bool>, int, int, int) ;
};

vector<int> getRegion (string str);

string asRegion (int chr, int pos, int end);

int read_tabixed_vcf(string &vcf_path, string &region, int &region_mode, int &one_vs_all, targetinfo &target, gdata &gdat, snpinfo &sinfo, idata &idat, int &n_haps);

int read_tabixed_m3vcf(string &m3vcf_path, string &region, int &region_mode, int &one_vs_all, targetinfo &target, gdata &gdat, snpinfo &sinfo, idata &idat, hdata &hdat, int &n_haps);

