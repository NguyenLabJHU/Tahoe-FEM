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
    xmin = 1732.8;  xmax = 13867.02;	// the coordinates of the computational domain
    ymin = 1894.1;  ymax =13746.81 ;
    zmin = -20720.4;  zmax = -2621.5;
   REAL RVE_xmin, RVE_xmax, RVE_ymin, RVE_ymax, RVE_zmin, RVE_zmax;	// the coordinates of each RVE
   REAL RVE_size = 6000;	// the size of the RVE
   REAL RVE_inter = 0.2*RVE_size;	// the space interval of two neighbor RVE
   RVE_xmin = 4456; RVE_xmax = 11140;
   RVE_ymin = 4456; RVE_ymax = 11140;
    RVE_zmin = -20141; RVE_zmax = -15596;


    char stepsstr[6];
    char stepsfp[50];
    // loop over all RVEs in the computational domain
    while(RVE_ymin <= ymax-RVE_size){
	RVE_ymax = RVE_ymax+RVE_size;
	while(RVE_zmin <= zmax-RVE_size){
	    RVE_zmax = RVE_zmax+RVE_size;

	    // construct the name of output file for each RVE
	    sprintf(stepsstr, "%05d", int(RVE_xmin)); 
	    strcpy(stepsfp, "results"); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    sprintf(stepsstr, "%05d", int(RVE_xmax)); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    sprintf(stepsstr, "%05d", int(RVE_ymin)); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    sprintf(stepsstr, "%05d", int(RVE_ymax)); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    sprintf(stepsstr, "%05d", int(RVE_zmin)); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);
	    sprintf(stepsstr, "%05d", int(RVE_zmax)); strcat(stepsfp, "_"); strcat(stepsfp, stepsstr);

	    assembly sample;
   	    sample.calculateGranularStrain(8, 24735, "F35_15.txt", stepsfp, RVE_xmin, RVE_xmax, RVE_ymin, RVE_ymax, RVE_zmin, RVE_zmax);	// F35_400.txt is the input file

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
