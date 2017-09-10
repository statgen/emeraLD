#include "processGenotypes.hpp"
#include "calcLD.hpp"
#include <stdlib.h>
#include <getopt.h>
#include <vector>
#include <random>

using namespace std;

void print_usage() {
    
    cerr << "\ncommand line options\n";
    cerr << "\tinput\n";
    cerr << "\t\t--in STR : vcf or m3vcf input file (bgzipped & tabixed)\n";
//    cerr << "\t\t--stdin : read from stdin rather than file\n";
    cerr << "\toutput\n";
    cerr << "\t\t--out STR : output file prefix \n";
    cerr << "\t\t--stdout : print output to stdout rather than file\n";
	cerr << "\t\t--dstats : print D and D' statistics in addition to R\n";
    cerr << "\t\t--matrix : print LD matrix rather than longtable\n";
    cerr << "\t\t--extra : print RSID, REF, and ALT in addition to CHR, POS\n";
    cerr << "\toptions\n";
    cerr << "\t\t--help : print this message and exit\n";
	cerr << "\t\t--region STR : calculate LD for SNPs in region (chr:start-end) \n";
    cerr << "\t\t--snp STR : only print pairwise LD for specified SNP (chr:pos) \n";
    cerr << "\t\t--rsid STR : only print pairwise LD for specified SNP \n";
    cerr << "\t\t--window INT : only calculate LD between SNPs within specified bp window (default: 1Mbp)\n";
    cerr << "\t\t--threshold DOUBLE : only print LD if abs(LD) > threshold (default: 1e-5)\n\n";
    //  cerr << "\t\t--nmax INT : LD precision parameter (default: 5000)\n\n";
}

int main (int argc, char *argv[]){
    cout.precision(5);
    
    int k = 0;
    
    vector<int> chrs, poss, dirs;
    vector<string> rsids, refs, alts;
    
    cerr << "emeraLD v0.1 (c) 2017 corbin quick (corbinq@gmail.com)\n";
    
    string infile = ""; 
	string outfile = "";
    
    double min_print = 0.00001;
    int max_dist = 1000000;
    int max_sample = 5000;
    
    string ref_rsid = "";
	string ref_snp_s = "";
	int ref_chr = -1;
    int ref_snp = -1;
    int j_ref = -1;
    
    int pstdin = -1;
    int pstdout = -1;
    
    int extra = -1;
	int extrastats = -1;
	
	//int use_tabix = -1;
	
	int region_mode = -1;
	string region = "";
	string region_e = "";
	int r_chr = -1;
	int r_start = -1;
	int r_start_e = -1;
	int r_end = -1;
	int r_end_e = -1;
	
    int help, matrix_out, one_vs_all;
    help = matrix_out = one_vs_all = -1;
    
    int opt = 0;
    
    static struct option long_options[] = {
    {"help",      no_argument,  NULL,  'h' },
	{"in",      required_argument,  NULL,  'i' },
	{"stdin", no_argument, NULL, 'x'},
	{"out", required_argument, NULL,  'o' },
	{"stdout", no_argument, NULL, 'z'},
	{"window",    required_argument, NULL,  'w' },
	{"region",    required_argument, NULL,  'r' },
	{"threshold",    required_argument, NULL,  't' },
	{"matrix",   no_argument, NULL,  'm' },
	{"dstats", no_argument, NULL, 'p'},
	{"snp",    required_argument, NULL,  's' },
	{"rsid",    required_argument, NULL,  'd' },
	{"nmax",    required_argument, NULL,  'n' },
	{"extra",    no_argument, NULL,  'e' }
	};
	int long_index =0;
	while ((opt = getopt_long(argc, argv,"hi:xo:zw:t:ms:",
	long_options, &long_index )) != -1) {
		switch (opt) {
			case 'h' : help = 1;
			break;
			case 'i' : infile = optarg;
			break;
			case 'x' : infile = "STDIN"; pstdin = 1;
			break;
			case 'o' : outfile = optarg;
			break;
			case 'z' : outfile = "STDOUT"; pstdout = 1;
			break;
			case 'w' : max_dist = atoi(optarg);
			break;
			case 't' : min_print = atof(optarg);
			break;
			case 'm' : matrix_out = 1;
			break;
			case 's' : one_vs_all = 1; ref_snp_s = optarg;
			break;
			case 'p' : extrastats = 1; 
			break;
			case 'd' : one_vs_all = 1; ref_rsid = optarg;
			break;
			case 'r' : region_mode = 1; region = optarg;
			break;
			case 'n' : max_sample = atoi(optarg);
			break;
			case 'e' : extra = 1;
			break;
			case '?': 
			fprintf(stderr, "\ninput error: '%c'.\n", optopt);
			help = 1;
			break;
		}
	}

	int n_haps;

	int matches = 0;

	int pillow = 100000;
	
	vector<vector<int>> carriers;
	vector<vector<bool>> genotypes;
	string ref_ref, ref_alt;
	
	vector<vector<int>> HAP_IID; // maps individuals to haplotypes 
	vector<vector<int>> HAP_CTS; // haplotype counts per block 
	
	vector<int> N_HAP; // number of distinct haplotypes in block, length(HAP_CTS)
	vector<int> M_HAP; // number of SNPs within block 
	
	vector<int> BLOCK_ID; // SNP to haplotype block
	
	vector<vector<int>> HAP_MAP; // SNP to haplotype IDs
	
	vector<int> MACS;
	
	if(ref_snp_s != ""){
		vector<int> region_v = getRegion(ref_snp_s);
		if( region_v.size() < 2 ){
			cerr << "\n\tERROR: SNP format must follow --snp chr:pos\n\n";
			return 1;
		}
		ref_chr = region_v[0];
		ref_snp = region_v[1];
	}
	
	if( region_mode > 0 ){
		vector<int> region_v = getRegion(region);
		if( region_v.size() < 3 ){
			cerr << "\n\tERROR: region input format must follow --region chr:start-end\n\n";
			return 1;
		}
		r_chr = region_v[0];
		r_start = region_v[1];
		r_end = region_v[2];
		
		if(r_start > pillow + 1){
			r_start_e = r_start - pillow;
		}else{
			r_start_e = 0;
		}
		r_end_e = r_end + pillow; 
		
		region_e = asRegion(r_chr, r_start_e, r_end_e);
		
	}
	
	if( one_vs_all > 0 && region_mode < 0 && ref_snp > 0 ){
		region_mode = 1;
		vector<int> region_v = getRegion(ref_snp_s);
		r_chr = region_v[0];
		if( r_start > max_dist +1 ){
			r_start = region_v[1] - max_dist;
		}else{
			r_start = 0;
		}
		r_end = region_v[1] + max_dist;
		
		if(r_start > pillow + 1){
			r_start_e = r_start - pillow;
		}else{
			r_start_e = 0;
		}
		r_end_e = r_end + pillow; 
		
		region_e = asRegion(r_chr, r_start_e, r_end_e);
		
	} 
	
	if( help > 0 ){
		print_usage();
		return 0;
	}

	if( infile.length() < 1 && pstdin < 0){
		cerr << "\n\t-i input file not specified...\n";
	//	cerr << "\tspecify --stdin or -i STDIN to read from stdin\n";
		cerr << "\n\tuse --help to see more options\n\n";
		//    print_usage();
		return 1;
	}else if(infile.length() > 1){
		if( access( infile.c_str(), F_OK ) == -1 ){
			cerr << "\nERROR: input file '" << infile << "' does not exist...\n\n";
			return 1;
		}
		string infile_tbi = infile + ".tbi";
		if( access( infile_tbi.c_str(), F_OK ) == -1 ){
			cerr << "\nERROR: index file '" << infile_tbi << "' does not exist...\n";
			cerr << "\n\temeraLD requires bgzipped & tabixed input files";
			cerr << "\n\tuse 'tabix -p vcf " << infile << "' to generate index\n\n";
			return 1;
		}
	}

	if( infile.length() > 1 && pstdin > 0 && infile != "STDIN" && infile != "-"){
		cerr << "\nERROR: choose either --infile [" << infile << "] or --stdin, but not both \n";
		cerr << "\n\tuse --help to see more options\n\n";
		//    print_usage();
		return 1;
	}

	if( outfile.length() > 0 ){
		if(outfile == "STDOUT" || outfile == "-"){
			outfile = "STDOUT";
			pstdout = 1;
		}else{
			if(pstdout > 0){
				cerr << "\nERROR: choose either --out [" << outfile << "] or --stdout, but not both\n";
				cerr << "\n\tuse --help to see more options\n\n";
				return 1;
			}
			outfile = outfile + ".LD.txt";
		}
	}else{
		if( pstdout < 0 ){
			cerr << "\t-o output file not specified...\n";
			cerr << "\tuse --stdout or -o STDOUT to print output to stdout\n";
			cerr << "\n\tuse --help to see more options\n\n";
			return 1;
		}else{
			outfile = "STDOUT";
			pstdout = 1;
		}
	}

	int m3vcf = -1;
	if( infile.find("m3vcf") != string::npos ){
		m3vcf = 1;
		cerr << "\nreading from m3vcf file...\n";
	}else{
		cerr << "\nreading from vcf file...\n";
	}
	
	//use_tabix = 1;
	//if( region_mode > 0 && pstdin < 0 ){
	//	use_tabix = 1;
	//}

	if( m3vcf < 0 ){
		if( read_tabixed_vcf(infile, region_e, region_mode, r_start, r_end, one_vs_all, ref_rsid, ref_snp, ref_ref, ref_alt, j_ref, matches, carriers, genotypes, dirs, chrs, poss, rsids, refs, alts, k, n_haps, BLOCK_ID ) > 0 ){
			cerr << "\nERROR: check vcf file " << infile << "\n";
			return 1;
		}
	}else{
		if( read_tabixed_m3vcf(infile, region_e, region_mode, r_start, r_end, one_vs_all, ref_rsid, ref_snp, ref_ref, ref_alt, j_ref, matches, carriers, genotypes, dirs, chrs, poss, rsids, refs, alts, k, n_haps, HAP_IID, HAP_CTS, N_HAP, BLOCK_ID, HAP_MAP, MACS ) > 0 ){
			cerr << "\nERROR: check m3vcf file " << infile << "\n";
			return 1;
		}
	}
	
	if( one_vs_all > 0 ){
		
		if( matches < 1 ){
			cerr << "\nERROR: SNP " << chrs[0] << ":" << ref_snp << " (rsid " << ref_rsid <<") not found!\n";
			return 1;
		}
	}
	
	cerr << "\nprocessed genotype data for " << n_haps << " haplotypes...\n";
	cerr << "\ncalculating LD for " << genotypes.size() << " SNPs...\n\n";
	
	//ifclose(inStream);
	
	uniform_real_distribution<double> unif(0,1);
	default_random_engine re;
	vector<double> vunif (0.5*n_haps);
	for (int i = 0; i < 0.5*n_haps; i++) {
		vunif[i] = unif(re);
	}

	FILE *outf;
	if( pstdout > 0 || outfile == "STDOUT"){
		outf = stdout;
	}else{
		outf = fopen( outfile.c_str(),  "w" );
	}
	
	double r, d, dprime;

	if( one_vs_all > 0 ){

		if(extra > 0){
			fprintf (outf, "#CHR\tPOS1\tRSID1\tREF:ALT1\tPOS2\tRSID2\tREF:ALT2\t");
		}else{
			fprintf (outf, "#CHR\tPOS1\tPOS2\t");
		}
		if(extrastats > 0 ){
			fprintf(outf, "R\tRsq\tD\tDprime\n");
		}else{
			fprintf(outf, "R\tRsq\n");
		}
		
		for (int i = 0; i < k; i++) {
			if( abs(ref_snp - poss[i]) < max_dist && poss[i] != ref_snp ){
				if( BLOCK_ID[i] == BLOCK_ID[j_ref] ){
					corr_within (r, d, dprime, HAP_MAP[i], HAP_MAP[j_ref], HAP_CTS[BLOCK_ID[i]], MACS[i], MACS[j_ref], n_haps );
				}else{
					corr(r, d, dprime, carriers[i], carriers[j_ref], genotypes[i], genotypes[j_ref], n_haps, dirs[i], dirs[j_ref], max_sample, vunif );
				}
				if(extra > 0){
					fprintf (outf,"%u\t", chrs[i]);
					fprintf (outf,"%u\t", ref_snp);
					string outl = ref_rsid + "\t" + ref_ref + ":" + ref_alt + "\t";
					fprintf (outf,"%s",outl.c_str());
					fprintf (outf,"%u\t", poss[i]);
					outl = rsids[i] + "\t" + refs[i] + ":" + alts[i] + "\t";
					fprintf (outf,"%s",outl.c_str());
				}else{
					fprintf (outf,"%u\t%u\t%u\t",chrs[i], ref_snp, poss[i]);
				}
				if(extrastats > 0 ){
					fprintf (outf,"%.5f\t%.5f\t%.5f\t%.5f\n", r, pow(r,2), d, dprime );
				}else{
					fprintf (outf,"%.5f\t%.5f\n", r, pow(r,2) );
				}
			}
		}
		
	}else if( matrix_out > 0){
		
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < k; j++) {
				if( i == j ){
					r = 1;
				}else if( abs(poss[i]-poss[j]) < max_dist ){
					if( BLOCK_ID[i] == BLOCK_ID[j] ){
						corr_within (r, d, dprime, HAP_MAP[i], HAP_MAP[j], HAP_CTS[BLOCK_ID[i]], MACS[i], MACS[j], n_haps );
					}else{
						corr(r, d, dprime, carriers[i], carriers[j], genotypes[i], genotypes[j], n_haps, dirs[i], dirs[j], max_sample, vunif );
					}
				}
				if( j == k-1 ){
					fprintf (outf, "%.5f\n", r);
				}else{
					fprintf (outf, "%.5f\t", r);
				}
			}
		}
		
	}else{
		
		if(extra > 0){
			fprintf (outf, "#CHR\tPOS1\tRSID1\tREF:ALT1\tPOS2\tRSID2\tREF:ALT2\t");
		}else{
			fprintf (outf, "#CHR\tPOS1\tPOS2\t");
		}
		if(extrastats > 0 ){
			fprintf(outf, "R\tRsq\tD\tDprime\n");
		}else{
			fprintf(outf, "R\tRsq\n");
		}
		
		for (int i = 0; i < k - 1; i++) {
			for (int j = i + 1; j < k; j++) {
				if( abs((poss[i])-(poss[j])) < max_dist && rsids[i]!=rsids[j] ){
					if( BLOCK_ID[i] == BLOCK_ID[j] && carriers[i].size() > 25 && carriers[j].size() > 25 ){
						corr_within (r, d, dprime, HAP_MAP[i], HAP_MAP[j], HAP_CTS[BLOCK_ID[i]], MACS[i], MACS[j], n_haps );
					}else{
						corr(r, d, dprime, carriers[i], carriers[j], genotypes[i], genotypes[j], n_haps, dirs[i], dirs[j], max_sample, vunif );
					}
					if( abs(r) > min_print ){
						if(extra > 0){
							fprintf (outf, "%u\t", chrs[i]);
							fprintf (outf, "%u\t", poss[i]);
							string outl = rsids[i] + "\t" + refs[i] + ":" + alts[i] + "\t";
							fprintf (outf,"%s", outl.c_str());
							fprintf (outf, "%u\t", poss[j]);
							outl = rsids[j] + "\t" + refs[j] + ":" + alts[j] + "\t";
							fprintf (outf,"%s", outl.c_str());
						}else{
							fprintf (outf, "%u\t%u\t%u\t", chrs[i], poss[i], poss[j]);
						}
						if(extrastats > 0 ){
							fprintf (outf, "%.5f\t%.5f\t%.5f\t%.5f\n", r, pow(r,2), d, dprime );
						}else{
							fprintf (outf, "%.5f\t%.5f\n", r, pow(r,2) );
						}
					}
				}
			}
		}
		
	}

	fclose(outf);
	
	cerr << "done!! thanks for using emeraLD\n\n";

	return 0;
}

