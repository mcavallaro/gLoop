#include <cstdlib>
#include <gsl/gsl_rng.h>

using namespace std;

typedef struct Rates {
    double alpha;
    double beta;
    double decay;
    double l_on;
    double l_off;
    int isind;
    double loop;
    double delta;
//    double thinning;
} rates;

typedef struct Systems {
	int pol2;
	int mRNA;
    int current;
    int is_open;
	vector <double> part_sum;
} System;

int update(System*, int ,  int, Rates);
