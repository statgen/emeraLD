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
	
    cerr << "emeraLD v0.1 (c) 2017 corbin quick (corbinq@gmail.com)\n";
    
    string infile = ""; 
	string outfile = "";
    
    double min_print = 0.00001;
    int max_dist = 1000000;
    int max_sample = 1000;
    
	snpinfo sinfo;
	targetinfo target;
	gdata gdat;
	hdata hdat;
    
    int pstdin = -1;
    int pstdout = -1;
    
    int extra = -1;
	int extrastats = -1;
	
	//int use_tabix = -1;
	
	int region_mode = -1;
	string region = "";
	string region_e = "";
	
    int help, matrix_out, one_vs_all;
    help = matrix_out = one_vs_all = -1;
    
    int opt = 0;
	
	int wmin = 0;
    
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
	{"extra",    no_argument, NULL,  'e' },
	{"wmin",    required_argument, NULL,  'q' }
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
			case 's' : one_vs_all = 1; target.chrpos = optarg;
			break;
			case 'p' : extrastats = 1; 
			break;
			case 'd' : one_vs_all = 1; target.rsid = optarg;
			break;
			case 'r' : region_mode = 1; region = optarg;
			break;
			case 'n' : max_sample = atoi(optarg);
			break;
			case 'e' : extra = 1;
			break;
			case 'q' : wmin = atoi(optarg);
			break;
			case '?': 
			fprintf(stderr, "\ninput error: '%c'.\n", optopt);
			help = 1;
			break;
		}
	}

	int n_haps;
	
	if(target.chrpos != ""){
		vector<int> region_v = getRegion(target.chrpos);
		if( region_v.size() < 2 ){
			cerr << "\n\tERROR: SNP format must follow --snp chr:pos\n\n";
			return 1;
		}
		target.chr = region_v[0];
		target.pos = region_v[1];
	}
	
	if( region_mode > 0 ){
		vector<int> region_v = getRegion(region);
		if( region_v.size() < 3 ){
			cerr << "\n\tERROR: region input format must follow --region chr:start-end\n\n";
			return 1;
		}
	}
	
	if( one_vs_all > 0 && region_mode < 0 && target.pos > 0 ){
		region_mode = 1;
		vector<int> region_v = getRegion(target.chrpos);
		int r_start = region_v[1];
		if( r_start > max_dist +1 ){
			r_start = region_v[1] - max_dist;
		}else{
			r_start = 0;
		}
		region = asRegion(region_v[0], r_start, region_v[1] + max_dist);
	}
	
	if( matrix_out > 0 && one_vs_all > 0 ){
		cerr << "\n\t--matrix cannot be used with --rsid or --snp\n";
		cerr << "\n\tuse --help to see more options\n\n";
		//    print_usage();
		return 1;
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
		if( read_tabixed_vcf(infile, region, region_mode, one_vs_all, target, gdat, sinfo, n_haps) > 0 ){
			cerr << "\nERROR: check vcf file " << infile << "\n";
			return 1;
		}
	}else{
		if( read_tabixed_m3vcf(infile, region, region_mode, one_vs_all, target, gdat, sinfo, hdat, n_haps) > 0 ){
			cerr << "\nERROR: check m3vcf file " << infile << "\n";
			return 1;
		}
	}
	
	if( one_vs_all > 0 ){
		
		if( target.matches < 1 ){
			cerr << "\nERROR: SNP " << sinfo.chr[0] << ":" << target.pos << " (rsid " << target.rsid <<") not found!\n";
			return 1;
		}
	}
	
	cerr << "\nprocessed genotype data for " << n_haps << " haplotypes...\n";
	cerr << "\ncalculating LD for " << gdat.genotypes.size() << " SNPs...\n\n";
	
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
		
		for (int i = 0; i < sinfo.size(); i++) {
			if( abs(target.pos - sinfo.pos[i]) < max_dist && sinfo.pos[i] != target.pos ){
				if( gdat.block[i] == gdat.block[target.index] && gdat.carriers[i].size() > wmin && gdat.carriers[target.index].size() > wmin){
					corr_within (r, d, dprime, hdat.map[i], hdat.map[target.index], hdat.hcts[gdat.block[i]], gdat.mac[i], gdat.mac[target.index], n_haps , gdat.dir[i], gdat.dir[target.index]);
				}else{
					corr(r, d, dprime, gdat.carriers[i], gdat.carriers[target.index], gdat.genotypes[i], gdat.genotypes[target.index], n_haps, gdat.dir[i], gdat.dir[target.index], max_sample, vunif );
				}
				if(extra > 0){
					fprintf (outf,"%u\t", sinfo.chr[i]);
					fprintf (outf,"%u\t", target.pos);
					string outl = target.rsid + "\t" + target.ref + ":" + target.alt + "\t";
					fprintf (outf,"%s",outl.c_str());
					fprintf (outf,"%u\t", sinfo.pos[i]);
					outl = sinfo.rsid[i] + "\t" + sinfo.ref[i] + ":" + sinfo.alt[i] + "\t";
					fprintf (outf,"%s",outl.c_str());
				}else{
					fprintf (outf,"%u\t%u\t%u\t",sinfo.chr[i], target.pos, sinfo.pos[i]);
				}
				if(extrastats > 0 ){
					fprintf (outf,"%.5f\t%.5f\t%.5f\t%.5f\n", r, pow(r,2), d, dprime );
				}else{
					fprintf (outf,"%.5f\t%.5f\n", r, pow(r,2) );
				}
			}
		}
		
	}else if( matrix_out > 0){
		int last_i = 0;
		int last_j = 0;
		for (int i = 0; i < sinfo.size(); i++) {
			if( sinfo.pos[i] > last_i ){
				if(extra > 0){
					fprintf (outf, "%u\t", sinfo.chr[i]);
					fprintf (outf, "%u\t", sinfo.pos[i]);
					string outl = sinfo.rsid[i] + "\t" + sinfo.ref[i] + ":" + sinfo.alt[i] + "\t";
					fprintf (outf,"%s", outl.c_str());
				}
				last_i = sinfo.pos[i];
				last_j = 0;
				for (int j = 0; j < sinfo.size(); j++) {
					if( sinfo.pos[j] > last_j ){
						last_j = sinfo.pos[j];
						if( i == j ){
							r = 1;
						}else if( abs(sinfo.pos[i]-sinfo.pos[j]) < max_dist ){
							if( gdat.block[i] == gdat.block[j] && gdat.carriers[i].size() > wmin && gdat.carriers[j].size() > wmin ){
								corr_within (r, d, dprime, hdat.map[i], hdat.map[j], hdat.hcts[gdat.block[i]], gdat.mac[i], gdat.mac[j], n_haps, gdat.dir[i], gdat.dir[j] );
							}else{
								corr(r, d, dprime, gdat.carriers[i], gdat.carriers[j], gdat.genotypes[i], gdat.genotypes[j], n_haps, gdat.dir[i], gdat.dir[j], max_sample, vunif );
							}
						}
						if( j == sinfo.size()-1 ){
							fprintf (outf, "%.5f\n", r);
						}else{
							fprintf (outf, "%.5f\t", r);
						}
					}
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
		
		for (int i = 0; i < sinfo.size() - 1; i++) {
			for (int j = i + 1; j < sinfo.size(); j++) {
				if( abs((sinfo.pos[i])-(sinfo.pos[j])) < max_dist && sinfo.rsid[i]!=sinfo.rsid[j] ){
					if( gdat.block[i] == gdat.block[j] && gdat.carriers[i].size() > wmin && gdat.carriers[j].size() > wmin ){
						corr_within (r, d, dprime, hdat.map[i], hdat.map[j], hdat.hcts[gdat.block[i]], gdat.mac[i], gdat.mac[j], n_haps, gdat.dir[i], gdat.dir[j] );
					}else{
						corr(r, d, dprime, gdat.carriers[i], gdat.carriers[j], gdat.genotypes[i], gdat.genotypes[j], n_haps, gdat.dir[i], gdat.dir[j], max_sample, vunif );
					}
					if( abs(r) > min_print ){
						if(extra > 0){
							fprintf (outf, "%u\t", sinfo.chr[i]);
							fprintf (outf, "%u\t", sinfo.pos[i]);
							string outl = sinfo.rsid[i] + "\t" + sinfo.ref[i] + ":" + sinfo.alt[i] + "\t";
							fprintf (outf,"%s", outl.c_str());
							fprintf (outf, "%u\t", sinfo.pos[j]);
							outl = sinfo.rsid[j] + "\t" + sinfo.ref[j] + ":" + sinfo.alt[j] + "\t";
							fprintf (outf,"%s", outl.c_str());
						}else{
							fprintf (outf, "%u\t%u\t%u\t", sinfo.chr[i], sinfo.pos[i], sinfo.pos[j]);
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

