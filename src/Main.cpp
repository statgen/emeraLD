#include "version.hpp"
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
	cerr << "\t\t--no-phase : calculate genotype rather than haplotype LD (do not assume known phase)\n";
	cerr << "\t\t--phased : calculate phased haplotype LD (assume known genotype phase)\n";
	cerr << "\t\t--nmax : number of minor-allele carriers sampled (default: 1000)\n";
	cerr << "\t\t--region STR : calculate LD for SNPs in region (chr:start-end) \n";
	cerr << "\t\t--snp STR : only print pairwise LD for specified SNP (chr:pos) \n";
	cerr << "\t\t--epacts STR : only print pairwise LD for specified SNP (chr:pos_ref/alt)\n";
	cerr << "\t\t--rsid STR : only print pairwise LD for specified SNP \n";
	cerr << "\t\t--window INT : only calculate LD between SNPs within specified bp window (default: 1Mbp)\n";
	cerr << "\t\t--threshold DOUBLE : only print LD if abs(LD) > threshold (default: 1e-5)\n";
	cerr << "\tvariant filtering\n";
	cerr << "\t\t--mac INT : minimum minor allele count \n";
	cerr << "\t\t--max-mac INT : maximum minor allele count \n";
	cerr << "\tindividual filtering\n";
	cerr << "\t\t--include STR : one-column file of individual IDs to include\n";
	cerr << "\t\t--exclude STR : one-column file of individual IDs to exclude\n";
	//  cerr << "\t\t--nmax INT : LD precision parameter (default: 5000)\n\n";
	cerr << "\n";
}

int main (int argc, char *argv[]){
	//cout.precision(5);
	
	cerr << "emeraLD v" << VERSION << " (c) 2018 corbin quick (corbinq@gmail.com)\n";
	
	string infile = "";
	string outfile = "";
	
	double min_print = 0.00001;
	int max_dist = 1000000;
	
	snpinfo sinfo;
	targetinfo target;
	gdata gdat;
	hdata hdat;
	
	int pstdin = 0;
	int pstdout = 0;
	
	int extra = 0;
	int extrastats = 0;
	
	int input_error = 0;
	
	foptions fopts;
	
	fopts.max_mac = 100000000;
	fopts.min_mac = 1;
	
	fopts.m_fold = 0;
	fopts.m_pad = 50;
	
	fopts.mmac = 1000;
	fopts.max_sample = 1000;
	fopts.one_vs_all = 0;
	
	string keepfile = "";
	string exclfile = "";
	
	string unit = "haplotypes";
	
	idata idat;
	
	fopts.force_phased = 0;
	fopts.phased = 1;
	fopts.region_mode = 0;
	fopts.region = "";
	fopts.region_e = "";
	
	int help = 0, matrix_out = 0;
	
	int opt = 0;
	
	static struct option long_options[] = {
		{"help",      no_argument,  &help,  1},
		{"in",        required_argument,  NULL,  'i' },
		{"stdin",     no_argument, &pstdin, 1},
		{"out",       required_argument, NULL,  'o' },
		{"stdout",    no_argument, &pstdout, 1},
		{"window",    required_argument, NULL,  'w' },
		{"region",    required_argument, NULL,  'r' },
		{"threshold", required_argument, NULL,  't' },
		{"matrix",    no_argument, &matrix_out,  1 },
		{"dstats",    no_argument, &extrastats, 1},
		{"phased",    no_argument, &fopts.force_phased, 1},
		{"no-phase",  no_argument, &fopts.phased, 0},
		{"snp",       required_argument, NULL,  's' },
		{"rsid",      required_argument, NULL,  'd' },
		{"epacts",      required_argument, NULL,  'e' },
		{"nmax",      required_argument, NULL,  'n' },
		{"include",   required_argument,  NULL,  'k' },
		{"exclude",   required_argument,  NULL,  'v' },
		{"mac",       required_argument,  NULL,  'f' },
		{"max-mac",   required_argument,  NULL,  'a' },
		{"extra",     no_argument, &extra,  1 },
		{"m-pad",     required_argument, NULL, 'b'},
		{"m-fold",    required_argument, NULL, 'c'},
		{NULL,        0,                 NULL,  0 }
	};

	while ((opt = getopt_long(argc, argv, "hi:o:w:t:s:r:n:k:v:f:a:b:c:e:", long_options, NULL)) != -1) {
		switch (opt) {
			case 'i' : infile = optarg;
				break;
			case 'o' : outfile = optarg;
				break;
			case 'w' : max_dist = atoi(optarg);
				break;
			case 't' : min_print = atof(optarg);
				break;
			case 'e' :
				fopts.one_vs_all = 1;
				target.epacts = optarg;
				break;
			case 's' : fopts.one_vs_all = 1; target.chrpos = optarg;
				break;
			case 'd' : fopts.one_vs_all = 1; target.rsid = optarg;
				break;
			case 'r' : fopts.region_mode = 1; fopts.region = optarg;
				break;
			case 'n' : fopts.max_sample = atoi(optarg); fopts.mmac = atoi(optarg);
				break;
			case 'k' : keepfile = optarg;
				break;
			case 'v' : exclfile = optarg;
				break;
			case 'f' : fopts.min_mac = atoi(optarg);
				break;
			case 'a' : fopts.max_mac = atoi(optarg);
				break;
			case 'b' : fopts.m_pad = atoi(optarg);
				break;
			case 'c' : fopts.m_fold = atoi(optarg);
				break;
			case '?' : input_error++;
		}
	}

        if( help ){
		print_usage();
        	return 0;
        }
        
	if( input_error ){
		cerr << "\nERROR: unrecognized command line argument\n\n";
		cerr << "\nuse \"--help\" to see valid arguments and options\n\n";
		return 1;
	}
	
	if( pstdin ){
		infile = "STDIN";
	}
	if( pstdout ){
		outfile = "STDOUT";
	}
	
	int n_haps;

	setThresh(min_print);
	setMaxSample(fopts.max_sample);
	setPhased(fopts.phased);
	
	if(target.chrpos != ""){
		vector<string> region_v = getRegion(target.chrpos);
		if( region_v.size() < 2 ){
			cerr << "\n\tERROR: SNP format must follow --snp chr:pos\n\n";
			return 1;
		}
		target.chr = region_v[0];
		target.pos = stoi(region_v[1]);
	}
	else if (target.epacts != "") {
		target = parseEpactsVariant(target.epacts);
	}
	
	if( keepfile != "" ){
		if( exclfile != "" ){
			cerr << "\n\tERROR: you can specify --keep-ids or --excl-ids, but not both\n\n";
			return 1;
		}else{
			idat.open(keepfile, true);
		}
	}else{
		idat.open(exclfile, false);
	}
	
	if( infile.length() <= 1 && !pstdin ){
		cerr << "\n\t-i input file not specified...\n";
		//	cerr << "\tspecify --stdin or -i STDIN to read from stdin\n";
		cerr << "\n\tuse --help to see more options\n\n";
		//	print_usage();
		return 1;
	}else if( infile.length() > 1 ){
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
	
	if( infile.length() > 1 && pstdin && infile != "STDIN" && infile != "-"){
		cerr << "\nERROR: choose either --infile [" << infile << "] or --stdin, but not both \n";
		cerr << "\n\tuse --help to see more options\n\n";
		//    print_usage();
		return 1;
	}
	
	setOptions(fopts);
	idat.process(infile);
	
	if( fopts.region_mode ){
		vector<string> region_v = getRegion(fopts.region);
		if( region_v.size() < 3 ){
			cerr << "\n\tERROR: region input format must follow --region chr:start-end\n\n";
			return 1;
		}
	}
	
	if( fopts.one_vs_all && !fopts.region_mode && target.pos > 0 ){
		fopts.region_mode = 1;
		vector<string> region_v = getRegion(target.chrpos);
		int r_start = stoi(region_v[1]);
		if( r_start > max_dist +1 ){
			r_start = stoi(region_v[1]) - max_dist;
		}else{
			r_start = 0;
		}
		fopts.region = asRegion(region_v[0], r_start, stoi(region_v[1]) + max_dist);
	}
	
	if( matrix_out && fopts.one_vs_all ){
		cerr << "\n\t--matrix cannot be used with --rsid or --snp\n";
		cerr << "\n\tuse --help to see more options\n\n";
		//    print_usage();
		return 1;
	}
	
	if( outfile.length() > 0 ){
		if(outfile == "STDOUT" || outfile == "-"){
			outfile = "STDOUT";
			pstdout = 1;
		}else{
			if(pstdout){
				cerr << "\nERROR: choose either --out [" << outfile << "] or --stdout, but not both\n";
				cerr << "\n\tuse --help to see more options\n\n";
				return 1;
			}
			outfile = outfile + ".LD.txt";
		}
	}else{
		if( !pstdout ){
			cerr << "\t-o output file not specified...\n";
			cerr << "\tuse --stdout or -o STDOUT to print output to stdout\n";
			cerr << "\n\tuse --help to see more options\n\n";
			return 1;
		}else{
			outfile = "STDOUT";
			pstdout = 1;
		}
	}
	
	int m3vcf = 0;
	if( infile.find("m3vcf") != string::npos ){
		m3vcf = 1;
		cerr << "\nreading from m3vcf file...\n";
	}else{
		cerr << "\nreading from vcf file...\n";
	}
	
	setOptions(fopts);
	
	if( !m3vcf ){
		if( !fopts.phased ){
			cerr << "\nassuming unphased data (reporting diploid genotype LD)...\n";
		}
		if( read_tabixed_vcf(infile, target, gdat, sinfo, idat, n_haps, fopts.phased) > 0 ){
			cerr << "\nERROR: check vcf file " << infile << "\n";
			return 1;
		}
	}else{
		if( !fopts.phased ){
			cerr << "\nERROR: --nophase (genotype LD) currently only supported for VCF format\n";
			return 1;
		}
		if( read_tabixed_m3vcf(infile, target, gdat, sinfo, idat, hdat, n_haps) > 0 ){
			cerr << "\nERROR: check m3vcf file " << infile << "\n";
			return 1;
		}
	}
	
	if( !fopts.phased ){
		unit = "individuals";
		if( extrastats ){
			cerr << "\nERROR: --nophase (genotype LD) cannot be used with --dstats \n";
			return 1;
		}
	}
	
	if( fopts.one_vs_all ){
		
		if( target.matches < 1 ){
			cerr << "\nERROR: SNP " << sinfo.chr[0] << ":" << target.pos << " (rsid " << target.rsid <<") not found!\n";
			return 1;
		}
	}
	
	if( gdat.mac.size() < 1 ){
		if( fopts.region_mode ){
			cerr << "\nERROR: Found no SNPs in region " << fopts.region << " ... \n";
			cerr << "\nhint: make sure chr:pos formatting matches input file\n\n";
		}else{
			cerr << "\nERROR: For no SNPs in input file!\n";
			cerr << "\nhint: check input file formatting\n\n";
		}
		return 1;
	}
	
	setPhased(fopts.phased);
	setSize(n_haps);
	
	cerr << "\nprocessed genotype data for " << n_haps << " " << unit << "...\n";
	cerr << "\ncalculating LD for " << gdat.mac.size() << " SNPs...\n\n";
	
	// ifclose(inStream);
	
	FILE *outf;
	if( pstdout > 0 || outfile == "STDOUT"){
		outf = stdout;
	}else{
		outf = fopen( outfile.c_str(),  "w" );
	}
	
	double r, d, dprime;
	
	if( fopts.one_vs_all ){
		
		if( extra ){
			fprintf (outf, "#CHR\tPOS1\tRSID1\tREF:ALT1\tPOS2\tRSID2\tREF:ALT2\t");
		}else{
			fprintf (outf, "#CHR\tPOS1\tPOS2\t");
		}
		if( extrastats ){
			fprintf(outf, "R\tRsq\tD\tDprime\n");
		}else{
			fprintf(outf, "R\tRsq\n");
		}
		
		for (int i = 0; i < sinfo.size(); i++) {
			bool not_target = sinfo.pos[i] != target.pos && sinfo.ref[i] != target.ref && sinfo.alt[i] != target.alt;
			if( abs(target.pos - sinfo.pos[i]) < max_dist && not_target ){
				getCorr(r, d, dprime, i, target.index, gdat, hdat);
				if(  abs(r) > min_print  ){
					if( extra ){
						//fprintf (outf,"%u\t", sinfo.chr[i]);
						fprintf (outf,"%s\t", sinfo.chr[i].c_str());
						fprintf (outf,"%u\t", target.pos);
						string outl = target.rsid + "\t" + target.ref + ":" + target.alt + "\t";
						fprintf (outf,"%s",outl.c_str());
						fprintf (outf,"%u\t", sinfo.pos[i]);
						outl = sinfo.rsid[i] + "\t" + sinfo.ref[i] + ":" + sinfo.alt[i] + "\t";
						fprintf (outf,"%s",outl.c_str());
					}else{
						//fprintf (outf,"%u\t%u\t%u\t",sinfo.chr[i], target.pos, sinfo.pos[i]);
						fprintf (outf,"%s\t%u\t%u\t",sinfo.chr[i].c_str(), target.pos, sinfo.pos[i]);
					}
					if( extrastats ){
						fprintf (outf,"%.5f\t%.5f\t%.5f\t%.5f\n", r, pow(r,2), d, dprime );
					}else{
						fprintf (outf,"%.5f\t%.5f\n", r, pow(r,2) );
					}
				}
			}
		}
		
	}else if( matrix_out ){
		int last_i = 0;
		int last_j = 0;
		for (int i = 0; i < sinfo.size(); i++) {
			if( sinfo.pos[i] > last_i ){
				if( extra ){
					//fprintf (outf, "%u\t", sinfo.chr[i]);
					fprintf (outf, "%s\t", sinfo.chr[i].c_str());
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
							getCorr(r, d, dprime, i, j, gdat, hdat);
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
		
		if( extra ){
			fprintf (outf, "#CHR\tPOS1\tRSID1\tREF:ALT1\tPOS2\tRSID2\tREF:ALT2\t");
		}else{
			fprintf (outf, "#CHR\tPOS1\tPOS2\t");
		}
		if( extrastats ){
			fprintf(outf, "R\tRsq\tD\tDprime\n");
		}else{
			fprintf(outf, "R\tRsq\n");
		}
		
		for (int i = 0; i < sinfo.size() - 1; i++) {
			for (int j = i + 1; j < sinfo.size(); j++) {
				if( abs((sinfo.pos[i])-(sinfo.pos[j])) < max_dist && (sinfo.pos[i] < sinfo.pos[j] || sinfo.rsid[i]!=sinfo.rsid[j] ) ){
					getCorr(r, d, dprime, i, j, gdat, hdat);
					if( abs(r) > min_print ){
						if( extra ){
							//fprintf (outf, "%u\t", sinfo.chr[i]);
							fprintf (outf, "%s\t", sinfo.chr[i].c_str());
							fprintf (outf, "%u\t", sinfo.pos[i]);
							string outl = sinfo.rsid[i] + "\t" + sinfo.ref[i] + ":" + sinfo.alt[i] + "\t";
							fprintf (outf,"%s", outl.c_str());
							fprintf (outf, "%u\t", sinfo.pos[j]);
							outl = sinfo.rsid[j] + "\t" + sinfo.ref[j] + ":" + sinfo.alt[j] + "\t";
							fprintf (outf,"%s", outl.c_str());
						}else{
							//fprintf (outf, "%u\t%u\t%u\t", sinfo.chr[i], sinfo.pos[i], sinfo.pos[j]);
							fprintf (outf, "%s\t%u\t%u\t", sinfo.chr[i].c_str(), sinfo.pos[i], sinfo.pos[j]);
						}
						if( extrastats ){
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

