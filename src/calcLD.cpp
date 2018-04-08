#include "calcLD.hpp"
#include<iostream>

using namespace std;

double THRESH;
int max_sample;
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

double sizeIntersection(vector<int> &iid1, double &mac1, vector<bool> &genotypes2, int &n, vector<int> &sid1 ){
	double out = 0;
	if( mac1 > max_sample && max_sample > 0 ){
		int mmac = mac1;
		for(int k = 0; k < max_sample; k++){
			if( genotypes2[ sid1[k] ] ){
				out++;
			}
		}
		return (out/max_sample) * (mac1/n);
	}else{
		for(int k = 0; k < mac1; k++){
			if( genotypes2[ iid1[k] ] ){
				out++;
			}
		}
		return out/n;
	}
}

void corr (double &R, double &D, double &DPRIME, vector<int> &G1, vector<int> &G2, vector<bool> &geno1, vector<bool> &geno2, int &N, int &D1, int &D2, vector<int> &S1, vector<int> &S2 ) {
	
	double SS1 = G1.size(); double SS2 = G2.size();
	if( N == 0 || SS1 == 0 || SS2 == 0 || SS1 == N || SS2 == N ){
		R = 0;
		D = 0;
		DPRIME = 0;
		return;
	}
	
	double EG1 = SS1 / N;
	double EG2 = SS2 / N;
	double E12;
	
	if( max( EG1 * (1-EG2) / ((1-EG1)*EG2), EG2 * (1-EG1) / ((1-EG2)*EG1) ) < THRESH ){
		R = 0;
		D = 0;
		DPRIME = 0;
		return;
	}
	
	if( SS1 < SS2 ){
		E12 = sizeIntersection(G1, SS1, geno2, N, S1);
	}else{
		E12 = sizeIntersection(G2, SS2, geno1, N, S2);
	}
	
	double SD1 = sqrt(EG1)*sqrt(1 - EG1);
	double SD2 = sqrt(EG2)*sqrt(1 - EG2);
	double C12 = E12 - EG1*EG2;
	
	R = 0;
	
	D = D1 * D2 * C12;
	
	if( D < 0 ){
		if( EG1 * EG2 < (1-EG1)*(1-EG2)){
			DPRIME = abs(D)/(EG1 * EG2);
		}else{
			DPRIME = abs(D)/((1-EG1)*(1-EG2));
		}
	}else{
		if( (1-EG1) * EG2 < EG1*(1-EG2)){
			DPRIME = D/((1-EG1) * EG2);
		}else{
			DPRIME = D/(EG1 * (1-EG2));
		}
	}
	
	if( SD1 * SD2 > 0 ){
		R  = D / ( SD1 * SD2 );
	}
	
	if( abs(R) > 1 ){
		R = (R > 0) ? 1 : -1;
	}
}

void corr_within (double &R, double &D, double &DPRIME, vector<int> &idx_i, vector<int> &idx_j, vector<int> &hac, int &m_i, int &m_j, int &n, int &D1, int &D2){
	
	if( m_i == 0 || m_i == n || m_j == 0 || m_j == n ){
		R = 0;
		D = 0;
		DPRIME = 0;
		return;
	}
	
	double pi = m_i/(double)n;
	double pj = m_j/(double)n;
	double pij = 0;
	int x_i = 0; int x_j = 0;
	while(x_i < idx_i.size() && x_j < idx_j.size() ){
		if( idx_i[x_i] > idx_j[x_j] ){
			x_j++;
		}else if( idx_i[x_i] < idx_j[x_j] ){
			x_i++;
		}else if( idx_i[x_i] == idx_j[x_j] ){
			pij += hac[ idx_i[x_i] ];
			x_i++;
			x_j++;
		}
	}
	pij = pij / n;
	if( pij > pi ){
		pij = pi;
	}
	if( pij > pj ){
		pij = pj;
	}
	D = D1 * D2 * (pij - pi*pj );
	R = D/sqrt(pi * (1- pi) * pj * (1 - pj));
	
	if( D < 0 ){
		if( pi * pj < (1-pi)*(1-pj)){
			DPRIME = abs(D)/(pi * pj);
		}else{
			DPRIME = abs(D)/((1-pi)*(1-pj));
		}
	}else{
		if( (1-pi) * pj < pi*(1-pj)){
			DPRIME = D/((1-pi) * pj);
		}else{
			DPRIME = D/(pi * (1-pj));
		}
	}
}

void corr_unph (double &R, double &D, double &DPRIME, int i, int j, gdata& GDATA) {
	if( GDATA.p_gts[i].p1 + GDATA.p_gts[i].p2 > GDATA.p_gts[j].p1 + GDATA.p_gts[j].p2 ){
		int i_orig = i;
		i = j;
		j = i_orig;
	}
	
	DPRIME = -99;
	
	double EGj_1 = 0;
	double EGj_2 = 0;
	for(vector<int>::iterator k = GDATA.hets_ss[i].begin(); k != GDATA.hets_ss[i].end(); ++k) {
		EGj_1 += (double) GDATA.ugenos[j][*k];
	}
	if( GDATA.hets_ss[i].size() > 0 ){
		EGj_1 /= (double) GDATA.hets_ss[i].size();
	}
	for(vector<int>::iterator k = GDATA.homs_ss[i].begin(); k != GDATA.homs_ss[i].end(); ++k) {
		EGj_2 += (double) GDATA.ugenos[j][*k];
	}
	if( GDATA.homs_ss[i].size() > 0 ){
		EGj_2 /= (double) GDATA.homs_ss[i].size();
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
	D = GDATA.p_gts[i].p1*EGj_1 + 2*GDATA.p_gts[i].p2*EGj_2 - 4*p_i*p_j;
	D *= GDATA.dir[i]*GDATA.dir[j];
	R = D / ( s_i * s_j ); 
}

