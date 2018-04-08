#include "processGenotypes.hpp"

using namespace std;

foptions fopts;

int m_pad = 100;
int m_fold = 5;

int check_phase = 1;

uniform_real_distribution<double> unif(0,1);
random_device rd;
default_random_engine re(rd());

void setOptions(foptions in)
{
	fopts = in;
}

void snpinfo::push(string ch, int po, string rs, string rf, string al)
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

vector<int> getSubset(vector<int> &x, int &mac){
	if( mac > fopts.mmac + m_pad && fopts.mmac > 1 ){
		vector<int> out;
		if( mac > m_fold*fopts.mmac ){
			for(int k = 0; k < fopts.mmac; k++){
				out.push_back(x[ (int) ( unif(re)*(mac-1)+0.5) ]);
			}
			sort(out.begin(), out.end());
		}else{
			out = x;
			shuffle ( out.begin(), out.end(), re );
			out.resize(fopts.mmac);
		}
		sort(out.begin(), out.end());
		return out;
	}else{
		return x;
	}
}

vector<int> sampleVec(vector<int> &x, int &n){
	if( x.size() > n + m_pad ){
		vector<int> out;
		if( x.size() > m_fold*n ){
			for(int k = 0; k < n; k++){
				out.push_back(x[ (int) ( unif(re)*(n-1)+0.5) ]);
			}
		}else{
			out = x;
			shuffle ( out.begin(), out.end(), re );
			out.resize(n);			
		}
		sort(out.begin(), out.end());
		return out;
	}else{
		return x;
	}
}

void gdata::push (vector<int> ca, vector<bool> ge, int di, int ma, int bl)
{
	carriers.push_back(ca);
	subsample.push_back(getSubset(ca,ma));
	genotypes.push_back(ge);
	dir.push_back(di);
	mac.push_back(ma);
	block.push_back(bl);
}


void gdata::push (vector<uint_fast8_t>g,vector<int>ht,vector<int>hm,double p0, double p1, double p2, int d, int m, int bl)
{
	double pt = p0 + p1 + p2;
	p0 /= pt;
	p1 /= pt;
	p2 /= pt;
	for(int k = 0; k < hm.size(); k++){
		g[hm[k]] = 2;
	}
	if( ht.size() + hm.size() + m_pad > fopts.mmac && fopts.mmac > 1 ){
		int n2 = (int)( fopts.mmac*(2*p2)/(2*p2+p1) +0.5f);
		int n1 = fopts.mmac - n2;
		hets_ss.push_back(sampleVec(ht, n1));
		homs_ss.push_back(sampleVec(hm, n2));
	}else{
		hets_ss.push_back(ht);
		homs_ss.push_back(hm);
	}
	ugenos.push_back(g);
	hets.push_back(ht);
	homs.push_back(hm);
	gprob gp;
	gp.p0 = p0;
	gp.p1 = p1;
	gp.p2 = p2;
	p_gts.push_back(gp);
	dir.push_back(d);
	mac.push_back(m);
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

bool idata::process(string input){
	if( filter_mode ){
		Tabix tfile(input);
		string header;
		tfile.getHeader(header);
		header = header.substr(0, header.find_last_of("\n"));
		header = header.substr(header.find_last_of("\n")+1);
		string id;
		istringstream iss(header);
		for(int i = 0; i < 9; i++){
			iss >> id;
		}
		//		cout << "Processing VCF IIDs ... \n";
		while(iss >> id){
			//cout << "\t" << id << "\t" << keep(id) << "\n";
			kcols.push_back(keep(id));
		}
	}
	return true;
}

bool idata::keep( string id ){
	if( filter_mode ){
		//cout << id << "\t" << id.substr(0, id.find("_HAP_")) << "\n";
		if( ids.count(id.substr(0, id.find("_HAP_") ) ) > 0 ){
			if( keep_mode ){
				return true;
			}else{
				return false;
			}
		}else{
			if( keep_mode ){
				return false;
			}else{
				return true;
			}
		}
	}else{
		return true;
	}
}

bool idata::keep( int i ){
	if( filter_mode ){
		if( i < kcols.size() ){
			return kcols[i];
		}else{
			return false;
		}
	}else{
		return true;
	}
}

void idata::open(string idpath, bool kmode){
	if(idpath == ""){
		filter_mode = false;
	}else{
		filter_mode = true;
		keep_mode = kmode;
		ifstream idfile(idpath);
		//	cout << "Processing ID file ... \n";
		string id;
		while( idfile >> id ){
			//			cout << "\t" << id << "\n";
			ids.insert(id.substr(0, id.find("_HAP_")));
		}
		idfile.close();
	}
}

vector<string> getRegion (string str){
	vector<string> v ;
	size_t prev_pos = 0, pos;
	string ss;
	while ((pos = str.find_first_of(":-,", prev_pos)) != std::string::npos)
	{
		if (pos > prev_pos){
			ss = str.substr(prev_pos, pos-prev_pos);
			//ss = ss.substr(ss.find("chr")+1,ss.size());
			//v.push_back(stoi(ss));
			v.push_back(ss);
		}
		prev_pos= pos+1;
	}
	if (prev_pos< str.length()){
		ss = str.substr(prev_pos, std::string::npos);
		//ss = ss.substr(ss.find("chr")+1,ss.size());
		//v.push_back(stoi(ss));
		v.push_back(ss);
	}
	return v;
}

string asRegion (int chr, int pos, int end){
	string out = to_string(chr) + ":" + to_string(pos) + "-" + to_string(end);
	return out;
}

string asRegion (string chr, int pos, int end){
	string out = chr + ":" + to_string(pos) + "-" + to_string(end);
	return out;
}

string asRegion (string chr, string pos, string end){
	string out = chr + ":" + pos + "-" + end;
	return out;
}

int read_tabixed_vcf(string &vcf_path, targetinfo &target, gdata &gdat, snpinfo &sinfo, idata &idat, int &n_haps, int& ph){
	
	Tabix tfile(vcf_path);
	
	/*string chr_pfx;
	tfile.getHeader(chr_pfx);
	tfile.getNextLine(chr_pfx);
	chr_pfx = chr_pfx.substr(0,chr_pfx.find_first_of("\t"));
	chr_pfx = chr_pfx.substr(0,chr_pfx.find("chr")+3);*/
	
	int r_start, r_end;
	r_start = r_end = -1;
	string r_chr;
	
	string region_e;
	if( fopts.region_mode ){
		vector<string> region_v = getRegion(fopts.region);
		r_chr = region_v[0];
		r_start = stoi(region_v[1]);
		r_end = stoi(region_v[2]);
	}
	
	//	region = chr_pfx + region;
	tfile.setRegion( fopts.region );
	
	string line;
	
	int k = 0;
	
	if(sinfo.size() > 0){
		k = sinfo.size();
	}
	
	while ( tfile.getNextLine(line) ) {
		
		if( line[0] != '#' && line.length() > 0 ){
			
			string rsid, ref, alt, qual, filter, info, format;
			
			string chr;
			int n, n1, n0;
			int pos;
			n = n1 = n0 = 0;
			vector<int> id_0, id_1;
			vector<bool> genov;
			
			double p0, p1, p2;
			p0 = p1 = p2 = 0;
			vector<uint_fast8_t> ug;
			vector<int> ht;
			vector<int> hm0;
			vector<int> hm2;
			
			istringstream iss(line);
			
			iss >> chr >> pos >> rsid >> ref >> alt >> qual >> filter >> info >> format;
			
			if( !fopts.region_mode || (pos >= r_start && pos <= r_end && chr == r_chr) ){
				
				if( fopts.one_vs_all ){
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
				int ii = 0;
				int ni = 0;
				int mx = idat.ids.size() > 1 ? idat.ids.size() : 1e7;
				while(iss >> geno && ni < mx ){
					if( idat.keep(ii) ){
						ni++;
						if( check_phase ){
							if( !fopts.force_phased ){
								if( geno[1] != '|' ){
									cerr << "\nNOTE: genotype data appear to be unphased\n";
									cerr << "      reporting genotype LD rather than haplotype LD\n";
									cerr << "      use \"--phased\" option to override this behaviour\n";
									ph = 1;
									fopts.phased = 1;
								}
							}
							check_phase = 0;
						}
						if( fopts.phased ){
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
						}else{
							if( geno[0]=='0' && geno[2]=='0' ){
								hm0.push_back(n);
								ug.push_back(0);
								p0++;
								n++;
								n0 += 2;
							}else if( ( geno[0]=='0' && geno[2]=='1' )  || ( geno[0]=='1' && geno[2]=='0' ) ){
								ht.push_back(n);
								ug.push_back(1);
								p1++;
								n++;
								n0 += 1;
								n1 += 1;
							}else if( geno[0]=='1' && geno[2]=='1' ){
								hm2.push_back(n);
								ug.push_back(0);
								p2++;
								n++;
								n1 += 2;
							}
						}
					}
					ii++;
				}
				if( min(n1, n0) >= fopts.min_mac && max(n1, n0) <= fopts.max_mac ){
					sinfo.push(chr, pos, rsid, ref, alt);
					if( fopts.phased ){
						if( n1 < n0 ){
							gdat.push( id_1, genov, 1, n1, k );
						}else{
							genov.flip();
							gdat.push( id_0, genov, -1, n0, k );
						}
					}else{
						if( hm2.size() < hm0.size() ){
							gdat.push(ug, ht, hm2, p0, p1, p2,  1, 2*p2 + p1, k);
						}else{
							gdat.push(ug, ht, hm0, p2, p1, p0, -1, 2*p0 + p1, k);
						}
					}
					n_haps = n;
					k++;
				}
			}
		}
	}
	return 0;
}

int read_tabixed_m3vcf(string &m3vcf_path, targetinfo &target, gdata &gdat, snpinfo &sinfo, idata &idat, hdata &hdat, int &n_haps){
	Tabix tfile(m3vcf_path);
	
	/*	string chr_pfx;
	tfile.getHeader(chr_pfx);
	tfile.getNextLine(chr_pfx);
	chr_pfx = chr_pfx.substr(0,chr_pfx.find_first_of("\t"));
	chr_pfx = chr_pfx.substr(0,chr_pfx.find("chr")+3);*/
	
	int r_start, r_end;
	r_start = r_end = -1;
	int pad = 100000;
	string r_chr;
	
	string region_e;
	if( fopts.region_mode ){
		vector<string> region_v = getRegion(fopts.region);
		r_chr = region_v[0];
		r_start = stoi(region_v[1]);
		r_end = stoi(region_v[2]);
		if(r_start > pad + 1){
			region_e = asRegion(r_chr, r_start - pad, r_end + pad);
		}else{
			region_e = asRegion(r_chr, 0, r_end + pad);
		}
	}
	
	//	region_e = chr_pfx + region_e;
	tfile.setRegion(region_e);
	
	string line;
	
	vector<int> h_iid;
	vector<int> h_cts;
	int m_block = 0;
	int begun = 0;
	int N_B = 0;
	int k = 0;
	
	if(sinfo.size() > 0){
		k = sinfo.size();
		N_B = gdat.block[k-1];
	}
	
	while( tfile.getNextLine(line) ){
		if( line[0] != '#' && line.length() > 0 ){
			
			string rsid, ref, alt, qual, filter, info, format;
			
			string chr;
			int pos, n1, n0;
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
				int ii = 0;
				int ni = 0;
				int mx = idat.ids.size() > 1 ? idat.ids.size() : 1e7;
				while(iss >> geno && ni < mx ){
					if( idat.keep(ii) ){
						ni++;
						genotype = stoi(geno);
						if(genotype > max_h ){
							max_h = genotype;
						}
						h_iid.push_back(genotype);
						n_haps++;
					}
					ii++;
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
				if( begun > 0 && ( !fopts.region_mode || (pos >= r_start && pos <= r_end && chr == r_chr) ) ){
					if( fopts.one_vs_all ){
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
					if( min(n1, n0) >= fopts.min_mac && max(n1, n0) <= fopts.max_mac ){
						sinfo.push(chr, pos, rsid, ref, alt);
						if( n1 < n0 ){
							gdat.push( id_1, genov, 1, n1, N_B );
							hdat.push_map(negatives);
						}else{
							genov.flip();
							gdat.push( id_0, genov, -1, n0, N_B );
							hdat.push_map(positives);
						}
						n_haps = h_iid.size();
						k++;
						m_block++;
					}
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


