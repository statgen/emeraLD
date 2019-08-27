#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include "../htslib/bgzf.h"
#include "../htslib/tabix.h"
#include <iostream>

#ifndef TABIXPP_H
#define TABIXPP_H

using namespace std;

class Tabix {

    tabix_t *t;
    ti_iter_t iter;
    const ti_conf_t *idxconf;
    int tid, beg, end;
    string firstline;

public:

    string filename;

    Tabix(void);
    Tabix(string& file);
    ~Tabix(void);

    void getHeader(string& header);
    bool setRegion(string& region);
    bool getNextLine(string& line);

};

#endif /* TABIXPP_H */
