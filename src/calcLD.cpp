#include "calcLD.hpp"
#include<iostream>

using namespace std;

double sizeIntersection(vector<int> &iid1, double &mac1, vector<bool> &genotypes2, int &n, int &maxmac, vector<double> &runif ){
    double out = 0;
    if( mac1 > maxmac && maxmac > 0 ){
        // subsampling occurs only if max( MAC_1, MAC_2 ) > maxmac
        int mmac = mac1;
        for(int k = 0; k < maxmac; k++){
            if( genotypes2[ iid1[ (int) (runif[k]*(mmac-1)+0.5) ] ] ){
                out++;
            }
        }
        return (out/maxmac) * (mac1/n);
    }else{
        for(int k = 0; k < mac1; k++){
            if( genotypes2[ iid1[k] ] ){
                out++;
            }
        }
        return out/n;
    }
}

void corr (double &R, double &D, double &DPRIME, vector<int> &G1, vector<int> &G2, vector<bool> &geno1, vector<bool> &geno2, int &N, int &D1, int &D2, int mmax, vector<double> &uvec ) {
    
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
    
    if( SS1 < SS2 ){
        E12 = sizeIntersection(G1, SS1, geno2, N, mmax, uvec);
    }else{
        E12 = sizeIntersection(G2, SS2, geno1, N, mmax, uvec);
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
	return;
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
	return;
}
