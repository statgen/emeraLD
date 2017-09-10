#include "tabix.hpp"

Tabix::Tabix(void) { }

Tabix::Tabix(string& file) {
    filename = file;
    const char* cfilename = file.c_str();
    struct stat stat_tbi,stat_vcf;
    char *fnidx = (char*) calloc(strlen(cfilename) + 5, 1);
    strcat(strcpy(fnidx, cfilename), ".tbi");
    if ( bgzf_check_bgzf(cfilename)!=1 )
    {
        cerr << "\n\t was bgzip used to compress " << file << "?\n";
        free(fnidx);
        exit(1);
    }
    // Common source of errors: new VCF is used with an old index
    stat(fnidx, &stat_tbi);
    stat(cfilename, &stat_vcf);
    if ( stat_vcf.st_mtime > stat_tbi.st_mtime )
    {
        cerr << "\n\t" << file << ".tbi is older than " << file << "\n";
		cerr << "\n\tuse 'tabix -p vcf -f " << file << "' to overwrite or reindex\n";
        free(fnidx);
        exit(1);
    }
    free(fnidx);

    if ((t = ti_open(cfilename, 0)) == 0) {
        cerr << "\n\tcannot open " << file << "\n";
        exit(1);
    }

    if (ti_lazy_index_load(t) < 0) {
        cerr << "\n\tcannot load " << file << ".tbi\n";
        exit(1);
    }

    idxconf = ti_get_conf(t->idx);

    // set up the iterator, defaults to the beginning
    iter = ti_query(t, 0, 0, 0);

}

Tabix::~Tabix(void) {
    ti_iter_destroy(iter);
    ti_close(t);
}


void Tabix::getHeader(string& header) {
    header.clear();
    ti_iter_destroy(iter);
    iter = ti_query(t, 0, 0, 0);
    const char* s;
    int len;
    while ((s = ti_read(t, iter, &len)) != 0) {
        if ((int)(*s) != idxconf->meta_char) {
            firstline = string(s); // stash this line
            break;
        } else {
            header += string(s);
            header += "\n";
        }
    }
}

bool Tabix::setRegion(string& region) {
    if (ti_parse_region(t->idx, region.c_str(), &tid, &beg, &end) == 0) {
        firstline.clear();
        ti_iter_destroy(iter);
        iter = ti_queryi(t, tid, beg, end);
        return true;
    } else return false;
}

bool Tabix::getNextLine(string& line) {
    const char* s;
    int len;
    if (!firstline.empty()) {
        line = firstline; // recovers line read if header is parsed
        firstline.clear();
        return true;
    }
    if ((s = ti_read(t, iter, &len)) != 0) {
        line = string(s);
        return true;
    } else return false;
}
