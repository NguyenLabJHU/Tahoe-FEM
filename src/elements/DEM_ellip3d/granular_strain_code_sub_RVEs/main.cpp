// This is used to calculate the granular strain measures based on
// the experimental data from Professor Alshibli's group.
// Here we first calculate the granular strain measures for all the 
// domain and we calculate spatial dvdx and eulerian strain by rate form
// directly in this code. After that we will use matlab code to plot 
// eulerian strain by rate form, Lagrangian/Eulerian strain by rate-form 
// deformation gradient and Hencky strain.
// September 30, 2014 By Boning Zhang

#include "realtypes.h"
#include "vec.h"
#include "assembly.h"
#include <vector>
#include <iostream>


int main(){

    REAL xmin, xmax, ymin, ymax, zmin, zmax;
    xmin = 2495.56;  xmax = 13400.5;	// the coordinates of the computational domain
    ymin = 2273.52;  ymax = 13749.2;
    zmin = -20844.1;  zmax = -276.193;
    REAL RVE_xmin, RVE_xmax, RVE_ymin, RVE_ymax, RVE_zmin, RVE_zmax;	// the coordinates of each RVE
    REAL RVE_inter = 1453.4;	// the space interval of two neighbor RVE
    RVE_xmin = 7054; RVE_xmax = 8508;
    RVE_ymin = ymin; RVE_ymax = ymin;
    RVE_zmin = zmin; RVE_zmax = zmin;


    char stepsstr[6];
    char stepsfp[50];
    // loop over all RVEs in the computational domain
    while(RVE_ymin<ymax){
	RVE_ymax = RVE_ymax+RVE_inter;
	while(RVE_zmin<zmax){
	    RVE_zmax = RVE_zmax+RVE_inter;

	    // construct the name of output file for each RVE
	    sprintf(stepsstr, "%05d", int(RVE_xmin)); 
	    strcpy(stepsfp, "results"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    sprintf(stepsstr, "%05d", int(RVE_xmax)); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    sprintf(stepsstr, "%05d", int(RVE_ymin)); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    sprintf(stepsstr, "%05d", int(RVE_ymax)); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    sprintf(stepsstr, "%05d", int(RVE_zmin)); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    sprintf(stepsstr, "%05d", int(RVE_zmax)); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);

	    assembly sample;
   	    sample.calculateGranularStrain(8, 24720, "F35_400.txt", stepsfp, RVE_xmin, RVE_xmax, RVE_ymin, RVE_ymax, RVE_zmin, RVE_zmax);	// F35_400.txt is the input file

std::cout << "RVE_xmin: " << RVE_xmin << "  RVE_xmax: " << RVE_xmax << std::endl;
std::cout << "RVE_ymin: " << RVE_ymin << "  RVE_ymax: " << RVE_ymax << std::endl;
std::cout << "RVE_zmin: " << RVE_zmin << "  RVE_zmax: " << RVE_zmax << std::endl << std::endl << std::endl; 

	    RVE_zmin = RVE_zmin+RVE_inter;	// go to next RVE
	}

	RVE_ymin = RVE_ymin+RVE_inter;
	RVE_zmin = zmin; RVE_zmax = zmin;
    }


//    assembly sample;
//    sample.calculateGranularStrain(8, 24720, "F35_400.txt", "results_sub_domain", RVE_xmin, RVE_xmax, ymin, ymax, zmin, zmax);

    return 0;

} // end main
