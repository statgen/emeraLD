#include "processGenotypes.hpp"
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <cmath>

using namespace std;

void setThresh( double &p);
void setMaxSample(int &m);
void setPhased(int&);

double sizeIntersection(vector<int> &iid1, double &mac1, vector<bool> &genotypes2, int &n, vector<int> &sid1 );

void corr (double &R, double &D, double &DPRIME, vector<int> &G1, vector<int> &G2, vector<bool> &geno1, vector<bool> &geno2, int &N, int &D1, int &D2, vector<int> &ss1, vector<int> &ss2 );

void corr_within (double &R, double &D, double &DPRIME, vector<int> &idx_i, vector<int> &idx_j, vector<int> &hac, int &m_i, int &m_j, int &n, int &D1, int &D2);

void corr_unph (double &R, double &D, double &DPRIME, int i, int j, gdata& GDATA);
