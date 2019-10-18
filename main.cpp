#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <cstdlib>
#include <vector>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm>
#include <ctime>
#include "header.h"

void usage(char *argv[]){
  printf("**************************************************************\n"
         "**                                                          **\n"
         "**                  Simulate mRNA population                **\n"
         "**                                                          **\n"
         "**************************************************************\n"
         "\n\n"
         " This is Free Software - You can use and distribute it under \n"
         " the terms of the GNU General Public License, version 3 or later\n\n"
         " (c) Massimo Cavallaro (m.cavallaro@warwick.ac.uk)\n\n");
  printf("Usage: %s [# simulation time] [DIM ENSEMBLE] [alpha] [beta] [decay] [l_on] [l_off] [loop(recycling rate)]\n\n" , argv[0]);
}

int main(int argc, char *argv[]) {

if (argc != 9){
    usage(argv);
    exit(1);
}

ofstream swi;
ofstream tran;
ofstream mRNA;  
ofstream ss;
ofstream PolIIss;
ofstream PolII;

char file_name_switch_events[150];
char file_name_transcription_events[150];
char file_name_mRNA[150];
char file_name_ss[150];
char file_name_PolIIss[150];
char file_name_PolII[150];

snprintf(file_name_switch_events, 150, "%s%s_%s_%s_%s_%s_%s_%s.dat", folder, "/trace/switch_", argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]);
snprintf(file_name_transcription_events, 150, "%s%s_%s_%s_%s_%s_%s_%s.dat", folder, "/trace/transcription", argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]);
snprintf(file_name_mRNA, 150, "%s%s_%s_%s_%s_%s_%s_%s.dat", folder, "/trace/mRNA_", argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]);
snprintf(file_name_ss, 150, "%s%s_%s_%s_%s_%s_%s_%s.dat", folder, "/trace/ss_mRNA_", argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]);
snprintf(file_name_PolIIss, 150, "%s%s_%s_%s_%s_%s_%s_%s.dat", folder, "/trace/Pol2ss_", argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]);
snprintf(file_name_PolII, 150, "%s%s_%s_%s_%s_%s_%s_%s.dat", folder, "/trace/Pol2_", argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]);

fprintf(stderr, "%s\n", file_name_ss);

mRNA.open(file_name_mRNA);
swi.open(file_name_switch_events);
tran.open(file_name_transcription_events);
PolII.open(file_name_PolII);

int num_of_rates = 7;
double T;
double t;
double R;
int DIM_ENSEMBLE;
int i;
int j;
gsl_rng *r;

r = gsl_rng_alloc(gsl_rng_mt19937);

gsl_rng_set (r,  1144); //time(NULL)

int status = 0;

T = atof(argv[1]);
DIM_ENSEMBLE  = atoi(argv[2]);

Rates br;

br.alpha = atof(argv[3]);
br.beta = atof(argv[4]);
br.decay = atof(argv[5]);
br.l_on = atof(argv[6]);
br.l_off = atof(argv[7]);
br.loop = atof(argv[8]);
br.delta = 1;
//br.thinning = 1;

vector<Systems> system(DIM_ENSEMBLE);

for(i=0; i<DIM_ENSEMBLE; i++){
    system[i].pol2=1;
    system[i].mRNA=0;
    system[i].part_sum.resize(num_of_rates);
    system[i].part_sum.assign(num_of_rates, 0);
    system[i].is_open = gsl_ran_bernoulli(r, 0.5);
    system[i].current = 0;
}

for(i=0; i<DIM_ENSEMBLE; i++){
    system[i].part_sum[0] = br.alpha; // pol2 is recruited
    system[i].part_sum[1] = system[i].part_sum[0] + (1 - br.loop) * br.beta * system[i].pol2 * system[i].is_open; // transcription
    system[i].part_sum[2] = system[i].part_sum[1] + br.decay * system[i].mRNA; //mRNA decay
    system[i].part_sum[3] = system[i].part_sum[2] + br.l_off * system[i].is_open;
    system[i].part_sum[4] = system[i].part_sum[3] + br.l_on * (1 - system[i].is_open);
    system[i].part_sum[5] = system[i].part_sum[4] + br.loop * br.beta * system[i].pol2 * system[i].is_open; // transcription (without loosing polymerase )
    system[i].part_sum[6] = system[i].part_sum[5] + br.delta * system[i].pol2;
}

/* * *
 * START THE SIMULATION
 */

for (i=0; i<DIM_ENSEMBLE; i++){
    t = 0;
    while (t < T) {

  //      br.thinning = 1; //int(sin(t/100.)+1);
        R = gsl_rng_uniform(r) * system[i].part_sum.back();
        j = upper_bound(system[i].part_sum.begin(), system[i].part_sum.end(), R) - system[i].part_sum.begin();
        status = update(&system[0], i, j, br);

        if (status!=EXIT_SUCCESS) {
            fprintf(stderr,"ERROR %f %d", t, j);
        }

        if(i==0){
            if(j==0){
                PolII << t << " " << system[0].pol2 << endl;
            }
            if (j==1){ // mRNA is created
                PolII << t << " " << system[0].pol2 << endl;
                mRNA << t << " " << system[0].mRNA << endl;
                tran << t << " " << system[0].current << endl;
            }
            else if(j==2){ // mRNA is destroied
                mRNA << t << " " << system[0].mRNA << endl;
            }
            else if(j==3){
                swi << t << " " << system[0].is_open << endl;
            }
            else if(j==4){
                swi << t << " " << system[0].is_open << endl;
            }
            else if (j==5){ // mRNA is created
                mRNA << t << " " << system[0].mRNA << endl;
                tran << t << " " << system[0].current << endl;
            }
            else if (j==6){
                PolII << t << " " << system[0].pol2 << endl;
            }
        }
        t = t + gsl_ran_exponential(r, 1./system[i].part_sum.back());

    }
}

double media_mRNA = 0;
double media_Pol2 = 0;
for (i=0; i<DIM_ENSEMBLE; i++){
    media_mRNA = media_mRNA + system[i].mRNA;
    media_Pol2 = media_Pol2 + system[i].pol2;
}
media_mRNA = media_mRNA / (double)DIM_ENSEMBLE;
media_Pol2 = media_Pol2 / (double)DIM_ENSEMBLE;

double mpol2 =  br.alpha / (br.delta + br.beta * (1 - br.loop));

printf("%s %f %f %s %f %f %f \n",
    "pol2",
    mpol2,
    /* br.alpha / (br.delta + br.beta),*/
    media_Pol2,
    "mRNA",
    br.beta * br.l_on / (br.l_on + br.l_off) / br.decay * 
    (br.alpha + br.loop * mpol2 * br.beta) / (br.delta + br.beta) ,    
    br.beta * mpol2 / br.decay * br.l_on / (br.l_on + br.l_off),
    media_mRNA);



ss.open(file_name_ss);
ss << system[0].mRNA;
for (i= 1; i<DIM_ENSEMBLE; i++){
    ss << ' ' << system[i].mRNA;
}
ss << endl;


PolIIss.open(file_name_PolIIss);
PolIIss << system[0].pol2;
for (i= 1; i<DIM_ENSEMBLE; i++){
    PolIIss << ' ' << system[i].pol2;
}
PolIIss << endl;

mRNA.close();
swi.close();
tran.close();
ss.close();
PolIIss.close();
PolII.close();

gsl_rng_free(r);

return EXIT_SUCCESS;
}
