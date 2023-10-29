#include "functions.h"
#include "classes.h"
#include "random.h"

using namespace std;

void initRand(Random &rnd);
double error(double av1, double av2, int n);

int main(int argc, char* argv[]){

    ifstream input("input.in");

    int thickness = 0;
    double r0 = 0.;
    double h;
    
    //parameters needed in input.in file
    input >> thickness;                   // thickness = 0 for slab model, 1 for exponential model, 2 for sech2 model
    input >> r0;                          // scale length of galactic disk
    input >> h;

    input.close();

    double alpha = 1/r0; 
    double hslab = 1.2 * pow(r0, 0.633);               
    double constant = 0.;                         //constant in front of the integral
    double delta = 0.;                            //parameter of metropolis algorithm

    Random rnd;          //random number generator
    initRand(rnd);

    distro dis(r0, thickness);  

    if(thickness == 0){
        constant = 1. /(alpha*hslab);
        delta = 6.;
    }
    if(thickness == 1){
        constant = 1.;
        delta = 6.;
    } 
    if(thickness == 2){
        constant = 1.;
        delta = 7.;
    }

    Metropolis metro(dis, rnd, 0.5, delta);   

    f bes();
    int nblocks = 100;                      //number of blocks
    int nthrows = 100000;                   //number of throws in each block
    BlockAverage av(bes, 1, nblocks);       //Class for block average

    ofstream data;                 //Output file
    data.open("data.out");
    int wd=18;

    double tilder = 0.;
    double tildermax = 15.;
    double increment = h;
    int nsteps = tildermax / increment;

    for(int i=1; i <= nsteps; i++){

        for(int iblk=1; iblk <= av.nblk; iblk++) {
    
            av.Reset(iblk);   //Reset block averages
            for(int istep=1; istep <= nthrows; istep++){ 

                metro.SamplingStep();
                //cout << "Acceptance rate: " << (double) metro.accepted/metro.attempted << endl;            //Debug
                av.Measure(tilder, metro.xold);
                av.Accumulate();    //Update block averages
            }
            av.Averages(iblk);   //Results for current block
        }  
 
        data << setw(wd) << tilder <<  setw(wd) << constant * tilder * av.glob_av[0]/(double)av.nblk << setw(wd) << constant * tilder * av.errore << endl;
        cout << "r = " << tilder << " done" << endl;

        tilder = (double) increment * i;
    }

    data.close();

    return 0;
}