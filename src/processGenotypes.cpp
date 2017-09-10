#include "processGenotypes.hpp"

using namespace std;

vector<int> getRegion (string str){
	vector<int> v ; 
	size_t prev_pos = 0, pos;
	while ((pos = str.find_first_of(":-,", prev_pos)) != std::string::npos)
	{
		if (pos > prev_pos){
			v.push_back(stoi(str.substr(prev_pos, pos-prev_pos)));
		}
		prev_pos= pos+1;
    }
    if (prev_pos< str.length()){
        v.push_back(stoi(str.substr(prev_pos, std::string::npos)));
	}
	return v;
}

string asRegion (int chr, int pos, int end){
	string out = to_string(chr) + ":" + to_string(pos) + "-" + to_string(end);
	return out;
}

int read_tabixed_vcf(string &vcf_path, string &region, int &region_mode, int &r_start, int &r_end, int &one_vs_all, string &ref_rsid, int &ref_snp, string &ref_ref,  string &ref_alt, int &j_ref, int &matches, vector<vector<int>> &carriers, vector<vector<bool>> &genotypes, vector<int> &dirs, vector<int> &chrs, vector<int> &poss, vector<string> &rsids, vector<string> &refs, vector<string> &alts, int &k, int &n_haps, vector<int> &BLOCK_ID){
	
	Tabix tfile(vcf_path);
	
	tfile.setRegion(region);

	string line;

	while ( tfile.getNextLine(line) ) {

		if( line[0] != '#' && line.length() > 0 ){
			
			string rsid, ref, alt, qual, filter, info, format;
			
			int chr, pos, n, n1, n0;
			n = n1 = n0 = 0;
			vector<int> id_0, id_1;
			vector<bool> genov;
			
			istringstream iss(line);
			
			iss >> chr >> pos >> rsid >> ref >> alt >> qual >> filter >> info >> format;
			
			if( region_mode < 0 || (pos >= r_start && pos <= r_end) ){
			
				if( one_vs_all > 0){
					if( rsid==ref_rsid || pos==ref_snp ){
						matches++;
						ref_ref = ref;
						ref_alt = alt;
						if(ref_rsid == ""){
							ref_rsid=rsid;
						}
						if( ref_snp < 0 ){
							ref_snp = pos;
						}
						j_ref = k;
						if( matches > 1 ){
							cerr << "\nERROR: found duplicate matches for target SNP " << chr << ":" << pos << " (rsid " << rsid <<")\n";
							return 1;
						}
					}
				}
				chrs.push_back(chr); poss.push_back(pos); rsids.push_back(rsid); refs.push_back(ref); alts.push_back(alt);
				
				string geno;
				while(iss >> geno){
					for (int idx = 0; idx < 3; idx += 2){
						if( geno[idx] == '0' ){
							id_0.push_back(n);
							genov.push_back(false);
							n0++;
							n++;
						}else if( geno[idx] == '1' ){
							id_1.push_back(n);
							genov.push_back(true);
							n1++;
							n++;
						}
					}
				}
				if( n1 < n0 ){
					carriers.push_back( id_1 );
					dirs.push_back(1);
				}else{
					carriers.push_back( id_0 );
					dirs.push_back(-1);
					genov.flip();
				}
				genotypes.push_back( genov );
				n_haps = n;
				BLOCK_ID.push_back(k);
				k++;
				
			}
		}
	}
	return 0;
}

int read_tabixed_m3vcf(string &m3vcf_path, string &region, int &region_mode, int &r_start, int &r_end, int &one_vs_all, string &ref_rsid, int &ref_snp, string &ref_ref,  string &ref_alt, int &j_ref, int &matches, vector<vector<int>> &carriers, vector<vector<bool>> &genotypes, vector<int> &dirs, vector<int> &chrs, vector<int> &poss, vector<string> &rsids, vector<string> &refs, vector<string> &alts, int &k, int &n_haps, vector<vector<int>> &HAP_IID, vector<vector<int>> &HAP_CTS, vector<int> &N_HAP, vector<int> &BLOCK_ID, vector<vector<int>> &HAP_MAP, vector<int> &MACS){
	
	Tabix tfile(m3vcf_path);
	
	tfile.setRegion(region);

	string line;

	vector<int> h_iid;
	vector<int> h_cts;
	int m_block = 0;
	int begun = 0;
	int N_B = 0;

	while( tfile.getNextLine(line) ){
		if( line[0] != '#' && line.length() > 0 ){
			
			string rsid, ref, alt, qual, filter, info, format;
			
			int chr, pos, n1, n0;
			n1 = n0 = 0;
			vector<int> id_0, id_1;
			vector<bool> genov;
			vector<int> h_map;
			
			if( line.find("BLOCK") != string::npos ){
				begun++;
				if( m_block > 0 ){
					N_B++;
					m_block = 0;
					HAP_IID.push_back(h_iid);
					HAP_CTS.push_back(h_cts);
					N_HAP.push_back(h_cts.size());
				}
				h_iid.clear();
				h_cts.clear();
				istringstream iss(line);
				string tmp;
				for(int i = 0; i < 9; i++){
					iss >> tmp;
				}
				string geno;
				int genotype;
				int max_h = 0;
				n_haps = 0;
				while(iss >> geno){
					genotype = stoi(geno);
					if(genotype > max_h ){
						max_h = genotype;
					}
					h_iid.push_back(genotype);
					n_haps++;
				}
				for (int i = 0; i < max_h; i++){
					h_cts.push_back(0);
				}
				for (int i = 0; i < h_iid.size(); i++){
					h_cts[h_iid[i]]++;
				}
			}else{
				istringstream iss(line);
				
				iss >> chr >> pos >> rsid >> ref >> alt >> qual >> filter >> info >> format;
				//	  if( ! keyExists( observed, rsid ) ){
				//	    observed.insert({rsid, true});
				if( begun > 0 && (region_mode < 0 || (pos >= r_start && pos <= r_end) ) ){
					if( one_vs_all > 0 ){
						if( rsid==ref_rsid || pos==ref_snp ){
							matches++;
							ref_ref = ref;
							ref_alt = alt;
							if(ref_rsid == ""){
								ref_rsid=rsid;
							}
							if( ref_snp < 0 ){
								ref_snp = pos;
							}
							j_ref = k;
							if( matches > 1 ){
								cerr << "\nERROR: found duplicate matches for target SNP " << chr << ":" << pos << " (rsid " << rsid << ")\n";
								return 1;
							}
						}
					}
					vector<int> positives;
					for(string::size_type i = 0; i < format.size(); ++i) {
						if(format[i]=='1'){
							positives.push_back(i);
						}
					}
					HAP_MAP.push_back(positives);
					BLOCK_ID.push_back(N_B);
					chrs.push_back(chr); poss.push_back(pos); rsids.push_back(rsid); refs.push_back(ref); alts.push_back(alt);
					for( int i = 0; i < h_iid.size(); i++ ){
						
						if(  binary_search( positives.begin(), positives.end(), h_iid[i] ) ){
							id_0.push_back(i);
							genov.push_back(false);
							n0++;
						}else{
							id_1.push_back(i);
							genov.push_back(true);
							n1++;
						}
					}
					MACS.push_back(n1);
					if( n1 < n0 ){
						carriers.push_back( id_1 );
						dirs.push_back(1);
					}else{
						carriers.push_back( id_0 );
						dirs.push_back(-1);
						genov.flip();
					}
					genotypes.push_back( genov );
					n_haps = n1 + n0;
					k++;
					m_block++;
				}
			}
		}
		//line.clear();
	}
	if( m_block > 0 ){
		HAP_IID.push_back(h_iid);
		HAP_CTS.push_back(h_cts);
		N_HAP.push_back(h_cts.size());
	}
	return 0;
}


