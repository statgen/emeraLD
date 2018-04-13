#include "calcLD.hpp"
#include<iostream>

using namespace std;

double THRESH;
int max_sample;
int N;
int phased;

void setThresh(double &p){
	THRESH = p;
}

void setMaxSample(int &m){
	max_sample = m;
}

void setPhased(int &p){
	phased = p;
}

void setSize(int& n){
	N = n;
}

bool filter(double& p, double& q){
	if( p > q ){
		return filter(q, p);
	}else{
		return ( p*(1-q)/((1-p)*q) <= THRESH );
	}
}

double sizeIntersection(int &mac1, haploVec &genotypes2, vector<int> &sid1 ){
	double out = 0;
	int sz = sid1.size();
	double ms = min(mac1, sz);
	for(int k = 0; k < ms; k++){
		if( genotypes2[ sid1[k] ] ){
			out++;
		}
	}
	return (double) mac1 * (out/ms);
}

void getDprime(double& d_prime, double& cov, double &pi, double &pj){
	if( cov < 0 ){
		if( pi * pj < (1-pi)*(1-pj)){
			d_prime = abs(cov)/(pi * pj);
			return;
		}else{
			d_prime = abs(cov)/((1-pi)*(1-pj));
			return;
		}
	}else{
		if( (1-pi) * pj < pi*(1-pj)){
			d_prime = cov/((1-pi) * pj);
			return;
		}else{
			d_prime = cov/(pi * (1-pj));
			return;
		}
	}
}

bool isMonorph(int& MAC1, int& MAC2){
	return ( min(MAC1, MAC2) <= 0 || max(MAC1, MAC2) >= N || N == 0 );
}

void corr (double &R, double &D, double &DPRIME, int &MAC1, int &MAC2, haploVec &geno1, haploVec &geno2, int &D1, int &D2, vector<int> &S1, vector<int> &S2 ) {
	double EG2 = (double) MAC1 / N;
	double EG1 = (double) MAC2 / N;

	if( filter(EG1, EG2) ){
		R = 0;
		D = 0;
		DPRIME = 0;
		return;
	}
	
	double E12;
	
	if( MAC1 <= MAC2 ){
		E12 = sizeIntersection(MAC1, geno2, S1)/N;
	}else{
		E12 = sizeIntersection(MAC2, geno1, S2)/N;
	}
	
	D = D1 * D2 * (E12 - EG1*EG2);	
	getDprime(DPRIME, D, EG1, EG2);
	
	double denom = sqrt(EG1*(1-EG1)*EG2*(1-EG2));
	if( denom > 0 ){
		R  = D / denom;	
		if( abs(R) > 1 ){
			R = (R > 0) ? 1 : -1;
		}
	}else{
		R = -99;
	}
	return;
}

void corr_within (double &R, double &D, double &DPRIME, const vector<int> &idx_i, const haploVec &pos_j, vector<int> &hac, int &m_i, int &m_j, int &D1, int &D2){
	
	if( isMonorph(m_i, m_j) ){
		R = -99;
		D = -99;
		DPRIME = -99;
		return;
	}
	
	double pi = m_i/(double)N;
	double pj = m_j/(double)N;
	
	if( filter(pi, pj) ){
		R = 0; D = 0; DPRIME = 0;
		return;
	}
	
	double pij = 0;
	for( const int& k : idx_i ){
		if( pos_j[k] ) pij += hac[ k ];
	}
	pij /= N;
	D = D1 * D2 * (pij - pi*pj );
	R = D/sqrt(pi * (1- pi) * pj * (1 - pj));
	getDprime(DPRIME, D, pi, pj);
	return;
}
/*
void corr_between (double &R, double &D, double &DPRIME, int i, int j, gdata& GDATA, hdata& HDATA) {
	double pi = GDATA.mac[i]/(double) N;
	double pj = GDATA.mac[j]/(double) N;
	if( filter(pi, pj) ){
		R = 0; D = 0; DPRIME = 0;
		return;
	}
	if(GDATA.block[i] != HDATA.b1){
		if( GDATA.block[i] == HDATA.b2 &&  GDATA.block[j] == HDATA.b1 ){
			swap(i,j);
		}else{
			HDATA.get_matrix(GDATA.block[i], GDATA.block[j]);
		}
	}
	double pij = 0;
	for(const int& a : HDATA.map[i] ){
		for(const int& b : HDATA.map[j] ){
			pij += HDATA.count_matrix[a][b];
		}
	}
	pij /= N;
	D = GDATA.dir[i] * GDATA.dir[j] * (pij - pi*pj );
	R = D/sqrt(pi * (1- pi) * pj * (1 - pj));
	getDprime(DPRIME, D, pi, pj);
	return;
}
*/
void corr_unph (double &R, double &D, double &DPRIME, int i, int j, gdata& GDATA) {
	double maf_1 = 0.5*GDATA.p_gts[i].p1 + GDATA.p_gts[i].p2;
	double maf_2 = 0.5*GDATA.p_gts[j].p1 + GDATA.p_gts[j].p2;
	maf_1 = (maf_1 <= 0.5 ? maf_1 : 1 - maf_1);
	maf_2 = (maf_2 <= 0.5 ? maf_2 : 1 - maf_2);
	if( filter(maf_1, maf_2) ){
		R = 0; D = 0; DPRIME = 0;
		return;
	}
	if( GDATA.p_gts[i].p1 + GDATA.p_gts[i].p2 > GDATA.p_gts[j].p1 + GDATA.p_gts[j].p2 ){
		swap(i, j);
	}
	
	GDATA.getCarriers(i);
	DPRIME = -99;
	
	double EGj_k1 = 0;
	double EGj_k2 = 0;
	for( const int&  k : GDATA.hets_ss[i] ) {
		EGj_k1 += (double) GDATA.ugenos[j][k];
	}
	if( GDATA.hets_ss[i].size() > 0 ){
		EGj_k1 /= (double) GDATA.hets_ss[i].size();
	}
	for( const int&  k : GDATA.homs_ss[i] ) {
		EGj_k2 += (double) GDATA.ugenos[j][k];
	}
	if( GDATA.homs_ss[i].size() > 0 ){
		EGj_k2 /= (double) GDATA.homs_ss[i].size();
	}
	
	double p_i = GDATA.p_gts[i].p2 + 0.5*GDATA.p_gts[i].p1;
	double p_j = GDATA.p_gts[j].p2 + 0.5*GDATA.p_gts[j].p1;

	double s_i = sqrt( 4*GDATA.p_gts[i].p2 + GDATA.p_gts[i].p1 - 4*p_i*p_i );
	double s_j = sqrt( 4*GDATA.p_gts[j].p2 + GDATA.p_gts[j].p1 - 4*p_j*p_j );
	
	if( s_i <= 0 || s_j <= 0 ){
		R = -99;
		D = -99;
		return;
	}
	D = GDATA.p_gts[i].p1*EGj_k1 + 2*GDATA.p_gts[i].p2*EGj_k2 - 4*p_i*p_j;
	D *= GDATA.dir[i]*GDATA.dir[j];
	R = D / ( s_i * s_j ); 
}

void getCorr(double &r, double &d, double &dprime, int i, int j, gdata& gdat, hdata& hdat){
	if( gdat.block[i] == gdat.block[j] ){
		bool do_swap = (hdat.map[i].size() > hdat.map[j].size());
		if( do_swap ) swap(i,j);
		corr_within (r, d, dprime, hdat.map[i], hdat.dense[j], hdat.hcts[gdat.block[i]], gdat.mac[i], gdat.mac[j], gdat.dir[i], gdat.dir[j]);
		if( do_swap ) swap(i,j);
	}else{
		if( phased ){
			if(max_sample > 1) gdat.getCarriers( gdat.mac[i] <= gdat.mac[j] ? i : j );
			corr(r, d, dprime, gdat.mac[i], gdat.mac[j], gdat.genotypes[i], gdat.genotypes[j], gdat.dir[i], gdat.dir[j],  gdat.subsample[i], gdat.subsample[j] );
		}else{
			corr_unph(r, d, dprime, i, j, gdat);
		}
	}
}

