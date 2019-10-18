#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <vector>
#include <cstdio>
#include <iomanip>
#include "header.h"

int update(Systems *system, int i, int j, Rates br){
	if (j == 0) { // pol2 is recruited
		system[i].pol2++;
	}
	else if (j == 1) {  // transcription
        system[i].mRNA++ ;
		system[i].current++ ;
        system[i].pol2--;
	}
	else if (j == 2){ //mRNA decay
		system[i].mRNA--;
	}
	else if (j == 3){
		system[i].is_open--;
	}
	else if (j == 4){
		system[i].is_open++;
	}
	else if (j == 5){ // transcription (without loosing polymerase )
		system[i].mRNA++;
		system[i].current++;
	}
	else if (j == 6){ // loose polymerase from a compartment
		system[i].pol2--;
	}
//     else if (j==7){
// //        do nothing
//        fprintf(stderr, "move rejected\n");
//     }
	else {
        fprintf(stderr,"error %d %d %f\n", i, j, system[i].part_sum.back() );
        exit(777);
    }

	//UPDATE THE VECTOR PART_SUM
	system[i].part_sum[0] = br.alpha; // pol2 is recruited
	system[i].part_sum[1] = system[i].part_sum[0] + (1 - br.loop) * br.beta * system[i].pol2 * system[i].is_open; // transcription
	system[i].part_sum[2] = system[i].part_sum[1] + br.decay * system[i].mRNA; //mRNA decay
	system[i].part_sum[3] = system[i].part_sum[2] + br.l_off * system[i].is_open;
	system[i].part_sum[4] = system[i].part_sum[3] + br.l_on * (1 - system[i].is_open);
	system[i].part_sum[5] = system[i].part_sum[4] + br.loop *  br.beta * system[i].pol2 * system[i].is_open; // transcription (without loosing polymerase )
	system[i].part_sum[6] = system[i].part_sum[5] + br.delta * system[i].pol2; // loose polymerase from compartment    
    // system[i].part_sum[7] = system[i].part_sum[6]; // + (1 - br.thinning) * (1 + br.loop) * br.beta * system[i].pol2 * system[i].is_open; // do nothing
	return EXIT_SUCCESS;
}
