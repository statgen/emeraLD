#include "processGenotypes.hpp"

using namespace std;

void snpinfo::push(int ch, int po, string rs, string rf, string al)  
{
	chr.push_back(ch);
	pos.push_back(po);
	rsid.push_back(rs);
	ref.push_back(rf);
	alt.push_back(al);
}

unsigned int snpinfo::size()  {
	return pos.size();
}

void gdata::push (vector<int> ca, vector<bool> ge, int di, int ma, int bl) 
{
	carriers.push_back(ca);
	genotypes.push_back(ge);
	dir.push_back(di);
	mac.push_back(ma);
	block.push_back(bl);
}

void hdata::push_block (vector<int> ii, vector<int> hc, int nh) 
{
	iids.push_back(ii);
	hcts.push_back(hc);
	nhap.push_back(nh);
}

void hdata::push_map (vector<int> ma) 
{
	map.push_back(ma);
}


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


int read_tabixed_vcf(string &vcf_path, string &region, int &region_mode, int &r_start, int &r_end, int &one_vs_all, targetinfo &target, gdata &gdat, snpinfo &sinfo, int &n_haps){
	Tabix tfile(vcf_path);
	
	tfile.setRegion(region);

	string line;
	
	int k = 0;

	while ( tfile.getNextLine(line) ) {

		if( line[0] != '#' && line.length() > 0 ){
			
			string rsid, ref, alt, qual, filter, info, format;
			
			int chr, n, n1, n0;
			int pos;
			n = n1 = n0 = 0;
			vector<int> id_0, id_1;
			vector<bool> genov;
			
			istringstream iss(line);
			
			iss >> chr >> pos >> rsid >> ref >> alt >> qual >> filter >> info >> format;
			
			if( region_mode < 0 || (pos >= r_start && pos <= r_end) ){
			
				if( one_vs_all > 0){
					if( rsid == target.rsid || pos == target.pos ){
						target.matches++;
						target.ref = ref;
						target.alt = alt;
						if(target.rsid == ""){
							target.rsid = rsid;
						}
						if( target.pos < 0 ){
							target.pos = pos;
						}
						target.index = k;
						if( target.matches > 1 ){
							cerr << "\nERROR: found duplicate target.matches for target SNP " << chr << ":" << pos << " (rsid " << rsid <<")\n";
							return 1;
						}
					}
				}
				
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
				sinfo.push(chr, pos, rsid, ref, alt);
				if( n1 < n0 ){
					gdat.push( id_1, genov, 1, n1, k );
				}else{
					genov.flip();
					gdat.push( id_0, genov, -1, n0, k );
				}
				n_haps = n;
				k++;
			}
		}
	}
	return 0;
}

int read_tabixed_m3vcf(string &m3vcf_path, string &region, int &region_mode, int &r_start, int &r_end, int &one_vs_all, targetinfo &target, gdata &gdat, snpinfo &sinfo, hdata &hdat, int &n_haps){
	Tabix tfile(m3vcf_path);
	
	tfile.setRegion(region);

	string line;

	vector<int> h_iid;
	vector<int> h_cts;
	int m_block = 0;
	int begun = 0;
	int N_B = 0;
	int k = 0;

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
					hdat.push_block(h_iid, h_cts, h_cts.size());
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
					if( one_vs_all > 0){
						if( rsid == target.rsid || pos == target.pos ){
							target.matches++;
							target.ref = ref;
							target.alt = alt;
							if(target.rsid == ""){
								target.rsid = rsid;
							}
							if( target.pos < 0 ){
								target.pos = pos;
							}
							target.index = k;
							if( target.matches > 1 ){
								cerr << "\nERROR: found duplicate target.matches for target SNP " << chr << ":" << pos << " (rsid " << rsid <<")\n";
								return 1;
							}
						}
					}
					vector<int> positives;
					vector<int> negatives;
					for(string::size_type i = 0; i < format.size(); ++i) {
						if(format[i]=='1'){
							positives.push_back(i);
						}else{
							negatives.push_back(i);
						}
					}
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
					sinfo.push(chr, pos, rsid, ref, alt);
					if( n1 < n0 ){
						gdat.push( id_1, genov, 1, n1, N_B );
						hdat.push_map(positives);
					}else{
						genov.flip();
						gdat.push( id_0, genov, -1, n0, N_B );
						hdat.push_map(negatives);
					}
					n_haps = h_iid.size();
					k++;
					m_block++;
				}
			}
		}
		//line.clear();
	}
	if( m_block > 0 ){
		hdat.push_block(h_iid, h_cts, h_cts.size());
	}
	return 0;
}


