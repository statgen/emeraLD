#include "processGenotypes.hpp"

using namespace std;

foptions fopts;

int load_sparse = 0;
int do_sample = 1;
int check_phase = 1;

int N_SS = 0;

uniform_real_distribution<double> unif(0,1);
random_device rd;
default_random_engine re(rd());

void setOptions(foptions in)
{
	fopts = in;
	if( !fopts.phased ){
		check_phase = 0;
	}
}

void diploVec::resize(int& n){
        het.resize(n);
        hom.resize(n);
}

void diploVec::push_back(int i){
	switch(i){
		case 0: 
			het.push_back(false);
			hom.push_back(false);
			break;
		case 1:
			het.push_back(true);
			hom.push_back(false);
			break;
		case 2:
			het.push_back(true);
			hom.push_back(true);
			break;	
	}
	return;
}

void diploVec::assign(int val, int& index){
	switch(val){
		case 0:
			//het[index] = false;
			//hom[index] = false;
			break;
		case 1:
			het[index] = true;
			//hom[index] = false;
			break;
		case 2:
			het[index] = true;
			hom[index] = true;
			break;
	}
	return;
}

unsigned int diploVec::size(){
	return het.size();
}

void diploVec::flip(){
	hom.flip();
	return;
}

vector<int> getSparse(haploVec v){
	vector<int> out;
	int n = v.find_first();
	while( n < v.size() ){
		out.push_back(n);
		n = v.find_next(n);
	}
	return out;
}

void snpinfo::push(string& ch, int& po, string& rs, string& rf, string& al)
{
	chr.emplace_back(ch);
	pos.emplace_back(po);
	rsid.emplace_back(rs);
	ref.emplace_back(rf);
	alt.emplace_back(al);
}

unsigned int snpinfo::size()  {
	return pos.size();
}

vector<int> getSubset(vector<int>& x, int& mac){
	if( mac > fopts.mmac + fopts.m_pad && fopts.mmac > 1 ){
		vector<int> out;
		if( mac > fopts.m_fold*fopts.mmac ){
			out.reserve(fopts.mmac);
			for(int k = 0; k < fopts.mmac; k++){
				out.emplace_back(x[ (int) ( unif(re)*(mac-1)+0.5) ]);
			}
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

vector<int> sampleVec(vector<int>& x, int& n){
	if( x.size() > n + fopts.m_pad ){
		vector<int> out;
		if( x.size() > fopts.m_fold*n ){
			out.resize(n);
			for(int k = 0; k < n; k++){
				out[k] = x[ (int) ( unif(re)*(n-1)+0.5) ];
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

void gdata::push (vector<int>& ca, haploVec& ge, int di, int ma, int bl)
{
//	carriers.emplace_back(ca);
	if( load_sparse && do_sample ){
		subsample.emplace_back(getSubset(ca,ma));
	}else{
		subsample.emplace_back(ca);
	}
	has_ss.push_back(load_sparse);
	genotypes.emplace_back(ge);
	dir.emplace_back(di);
	mac.emplace_back(ma);
	block.emplace_back(bl);
}


void gdata::push (diploVec& g,vector<int>& ht,vector<int>& hm,double& p0, double& p1, double& p2, int d, int m, int bl)
{
	double pt = p0 + p1 + p2;
	p0 /= pt;
	p1 /= pt;
	p2 /= pt;
	if( load_sparse ){ 
		if( ht.size() + hm.size() + fopts.m_pad > fopts.mmac && fopts.mmac > 1 ){
			int n2 = (int)( fopts.mmac*(2*p2)/(2*p2+p1) +0.5f);
			int n1 = fopts.mmac - n2;
			hets_ss.push_back(sampleVec(ht, n1));
			homs_ss.push_back(sampleVec(hm, n2));
		}else{
			hets_ss.push_back(ht);
			homs_ss.push_back(hm);
		}
	}else{
		hets_ss.push_back(ht);
		homs_ss.push_back(hm);
	}
	has_ss.push_back(load_sparse);
	ugenos.push_back(g);
//	hets.push_back(ht);
//	homs.push_back(hm);
	gprob gp;
	gp.p0 = p0;
	gp.p1 = p1;
	gp.p2 = p2;
	p_gts.push_back(gp);
	dir.push_back(d);
	mac.push_back(m);
	block.push_back(bl);
}

void gdata::getCarriers(int i){
	if( !has_ss[i] ){
		if( fopts.phased ){
			int k = genotypes[i].find_first();
			int n = 0;
			bool kp_all = !(mac[i] > fopts.mmac + fopts.m_pad && fopts.mmac > 1);
			double kp_prob = fopts.mmac/ (double) mac[i];
			subsample[i].reserve(min(mac[i], fopts.mmac) );
			int nt = min(mac[i], fopts.mmac + fopts.m_pad);
			while( n < nt && k < N_SS ){
				if( kp_all ){
					subsample[i].push_back(k);
				}else if( unif(re) <= kp_prob ){
					subsample[i].push_back(k);
				}
				n++;
				k = genotypes[i].find_next(k);
			}
		}else{
			int n12 = floor(N_SS*(1-p_gts[i].p0) + 0.01);
			bool kp_all = !(n12 > fopts.mmac + fopts.m_pad && fopts.mmac > 1);
			double kp_prob2 = 2*p_gts[i].p2/(2*p_gts[i].p2 + p_gts[i].p1);
			double kp_prob1 = (1 - kp_prob2);
			double NR = (kp_prob2*p_gts[i].p2 + kp_prob1*p_gts[i].p1)*N_SS;
			kp_prob2 *= fopts.mmac/NR;
			kp_prob1 *= fopts.mmac/NR;
			int k = ugenos[i].het.find_first();
			int n = 0;
			hets_ss[i].reserve( (int) NR*p_gts[i].p1 );
			homs_ss[i].reserve( (int) NR*p_gts[i].p2 );
			int nt = min(n12, fopts.mmac + fopts.m_pad);
			while( n < nt && k < N_SS ){
				n++;
				if( ugenos[i].hom[k] ){
					if( kp_all ){
						homs_ss[i].push_back(k);
					}else if( unif(re) <= kp_prob2 ){
						homs_ss[i].push_back(k);
					}
				}else{
					if( kp_all ){
						hets_ss[i].push_back(k);
					}else if( unif(re) <= kp_prob1 ){
						hets_ss[i].push_back(k);
					}
				}
				k = ugenos[i].het.find_next(k);
			}
			hets_ss[i].shrink_to_fit();
			homs_ss[i].shrink_to_fit();
		}
		has_ss[i] = true;
	}
	return;
}

void hdata::push_block (vector<int>& ii, vector<int>& hc, int nh)
{
	iids.push_back(ii);
	hcts.push_back(hc);
	nhap.push_back(nh);
}

void hdata::push_map (vector<int>& ma, haploVec& de)
{
	map.push_back(ma);
	dense.push_back(de);
}
/*
void hdata::get_matrix(int bi, int bj){
	if( b1==bi && b2==bj ) return;
	count_matrix.clear();
	count_matrix.reserve(nhap[bi]);
	int u_i = hcts[bi].size();
	int u_j = hcts[bj].size();
	for(int i=0; i < u_i; i++){
		count_matrix.emplace_back(vector<int>(u_j,0));
	}
	for(int i=0; i < N_SS; i++){
		count_matrix[iids[bi][i]][iids[bj][i]]++;
	}
	b1 = bi;
	b2 = bj;
}
*/
bool idata::process(string& input){
	Tabix tfile(input);
	string header;
	tfile.getHeader(header);
	header = header.substr(0, header.length() - 1);
	header = header.substr(header.find_last_of("\n")+1);
	string id;
	istringstream iss(header);
	for(int i = 0; i < 9; i++){
		iss >> id;
	}
	//		//cout << "Processing VCF IIDs ... \n";
	while(iss >> id){
		////cout << "\t" << id << "\t" << keep(id) << "\n";
		bool kp = keep(id);
		if(kp){
			N_SS++;
		}
		if( filter_mode ){
			kcols.push_back(kp);
		}
	}
	if( fopts.mmac <= 1 || fopts.mmac >= N_SS*0.05 ){
		load_sparse = 1;
	}
	if( fopts.mmac <= 1 ){
		do_sample = 0;
	}
	return true;
}

bool idata::keep( string& id ){
	if( filter_mode ){
		////cout << id << "\t" << id.substr(0, id.find("_HAP_")) << "\n";
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

void idata::open(string& idpath, bool kmode){
	if(idpath == ""){
		filter_mode = false;
	}else{
		filter_mode = true;
		keep_mode = kmode;
		ifstream idfile(idpath);
		//cout << "Processing ID file ... \n";
		string id;
		while( idfile >> id ){
			//cout << "\t" << id << "\n";
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
			vector<int> id_1;
			haploVec genov;
			
			double p0, p1, p2;
			p0 = p1 = p2 = 0;
			diploVec ug;
			vector<int> ht;
			vector<int> hm0;
			vector<int> hm2;
			
			if( fopts.phased & !check_phase ){
				genov.resize(N_SS);
				if(load_sparse){
					id_1.reserve(N_SS);
				}
			}else{
				ug.resize(N_SS);
				if(load_sparse){
					ht.reserve(N_SS);
					hm0.reserve(N_SS);
					hm2.reserve(N_SS);
				}
			}
			
			istringstream iss(line);
			
			iss >> chr >> pos >> rsid >> ref >> alt >> qual >> filter >> info >> format;
			
			if( !fopts.region_mode || (pos >= r_start && pos <= r_end && chr == r_chr) ){
				
				if( fopts.one_vs_all ){
					if (target.epacts != "") {
						// User specified a specific chr:pos_ref/alt and we want to match only against it.
						if (pos == target.pos && ref == target.ref && alt == target.alt) {
							target.matches++;
							target.index = k;
							target.rsid = rsid;
						}
					}
					else {
						// We have to match against rsid or just plain position in this case.
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
						}
					}

					if( target.matches > 1 ){
						cerr << "\nERROR: found duplicate target.matches for target SNP " << chr << ":" << pos << " (rsid " << rsid <<")\n";
						return 1;
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
									ug.resize(N_SS);
									if(load_sparse){
										ht.reserve(N_SS);
										hm0.reserve(N_SS);
										hm2.reserve(N_SS);
									}
								}else{
									N_SS *= 2;
									genov.resize(N_SS);
									if(load_sparse){
										id_1.reserve(N_SS);
									}
								}
							}else{
								N_SS *= 2;
								genov.resize(N_SS);
								if(load_sparse){
									id_1.reserve(N_SS);
								}
							}
							check_phase = 0;
						}
						if( fopts.phased ){
							for (int idx = 0; idx < 3; idx += 2){
								if( geno[idx] == '0' ){
									n0++;
									n++;
								}else if( geno[idx] == '1' ){
									if(load_sparse) id_1.emplace_back(n);
									genov[n] = true;
									n1++;
									n++;
								}
							}
						}else{
							if( geno[0]=='0' && geno[2]=='0' ){
								if(load_sparse) hm0.emplace_back(n);
								p0++;
								n++;
								n0 += 2;
							}else if( ( geno[0]=='0' && geno[2]=='1' )  || ( geno[0]=='1' && geno[2]=='0' ) ){
								if(load_sparse) ht.emplace_back(n);
								ug.assign(1,n);
								p1++;
								n++;
								n0 += 1;
								n1 += 1;
							}else if( geno[0]=='1' && geno[2]=='1' ){
								if(load_sparse) hm2.emplace_back(n);
								ug.assign(2,n);
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
							id_1.shrink_to_fit();
							gdat.push( id_1, genov, 1, n1, k );
						}else{
							genov.flip();
							if( load_sparse ){
								id_1.clear();
								id_1 = getSparse(genov);
							}
							gdat.push( id_1, genov, -1, n0, k );
						}
					}else{
						ht.shrink_to_fit();
						if( p2 < p0 ){
							hm2.shrink_to_fit();
							gdat.push(ug, ht, hm2, p0, p1, p2,  1, 2*p2 + p1, k);
						}else{
							hm0.shrink_to_fit();
							ug.flip();
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
			vector<int> id_1;
			haploVec genov;
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
				h_iid.reserve(N_SS);
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
				h_cts.resize(max_h+1);
				fill(h_cts.begin(), h_cts.end(), 0);
				for (const int& i : h_iid){
					h_cts[i]++;
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
					haploVec dense;
					vector<int> sparse1;
					int u_h = format.size();
					dense.resize(u_h);
					for(int i = 0; i < format.size(); ++i) {
						if(format[i]=='1'){
							dense[i] = true;
							sparse1.push_back(i);
						}
					}
					genov.resize(N_SS);
					if(load_sparse){
						id_1.reserve(N_SS);
					}
					for( int i = 0; i < h_iid.size(); i++ ){
						if(  !dense[h_iid[i]] ){
							n0++;
						}else{
							if(load_sparse) id_1.emplace_back(i);
							genov[i] = true;
							n1++;
						}
					}
					if( min(n1, n0) >= fopts.min_mac && max(n1, n0) <= fopts.max_mac ){
						sinfo.push(chr, pos, rsid, ref, alt);
						if( n1 < n0 ){
							id_1.shrink_to_fit();
							gdat.push( id_1, genov, 1, n1, N_B );
							hdat.push_map(sparse1, dense);
						}else{
							genov.flip();
							dense.flip();
							if(load_sparse){
								id_1.clear();
								id_1 = getSparse(genov);
							}
							sparse1.clear();
							sparse1 = getSparse(dense);
							gdat.push( id_1, genov, -1, n0, N_B );
							hdat.push_map(sparse1, dense);
						}
						n_haps = h_iid.size();
						k++;
						m_block++;
					}
				}
			}
		}
		line.clear();
	}
	if( m_block > 0 && h_cts.size() > 0){
		hdat.push_block(h_iid, h_cts, h_cts.size());
	}
	return 0;
}

