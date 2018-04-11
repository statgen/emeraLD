#include "processGenotypes.hpp"
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <cmath>

using namespace std;

void setThresh(double&);
void setMaxSample(int&);
void setPhased(int&);
void setSize(int&);

void getCorr(double&, double&, double&, int&, int&, gdata&, hdata&);
