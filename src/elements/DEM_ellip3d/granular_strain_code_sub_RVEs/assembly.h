#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "realtypes.h"
#include "vec.h"
#include "particle.h"
#include "matrix.h"
#include "cell.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cassert>
#include <utility>
#include <cmath>
#include <stdlib.h>
#include <math.h>


class assembly{

public:
    assembly(){
	TotalNum = 0;
	Volume = 0;
    }

    ~assembly(){
    	std::vector<particle*>::iterator pt;
    	std::vector<cell*>::iterator ct;
    	// it is important to release memory pointed to by pointers in the container,
    	// otherwise memory leaking occurs
    	for(pt = ParticleVec.begin(); pt != ParticleVec.end(); ++pt)
      	    delete (*pt);
    	for(ct = cellVec.begin(); ct != cellVec.end(); ++ct)
           delete (*ct);
    	for(ct = cellVec_init.begin(); ct != cellVec_init.end(); ++ct)
          delete (*ct);
    
    	// in case of consecutive simulations
    	ParticleVec.clear();
    	cellVec.clear();
    	cellVec_init.clear();
    }

    void calculateGranularStrain(int num_step, int total_num, const char* iniptclfile, const char* progressfile,
				       REAL RVE_xmin, REAL RVE_xmax, REAL RVE_ymin, REAL RVE_ymax, REAL RVE_zmin, REAL RVE_zmax);

    void createInputForQhull() const;	// create input file for Qhull, int means the step
    void callQhull() const;	// based on input_for_Qhull that is from createInputForQhull() to tessellate and output "tess_info"
    void readTesse_finite(const char* str);	// create cell information into std::vector<cell> cellVec from Qhull out put files
    void setNumberingOrder();	// ensure that the numbering of each tetrahedron in std::vector<cell> cellVec is couter-clockwise
    matrix getAverage_dvdx() const;	// calculate average spatial velocity gradient tensor, June 24, 2013
    REAL getCellVolume(const cell&) const;	// calculate the volume of cell
    matrix getdvdx_curr(const cell&) const;	// calculate spatial velocity gradient tensor for each tet, June 24, 2013
    matrix getBigB(const cell&) const;
    void updateParticle(int step, REAL time);	// update the particle position and velocity based on the experimental data
    void readSample(const char*);



private:
    // particles property
    int  TotalNum;      // total number of particles
    int  RVE_number;	// number of particles in RVE
    int  nstep;		// the number of steps
    std::vector<particle*> ParticleVec; // a vector of pointers, each pointing to a particle
    std::vector<particle*> ParticleVec_part;	// a vector of pointers, each pointing to a particle in RVE, Juen 17, 2013      
    // container property
    REAL Volume;      // volume of the specimen
    REAL initVolume;  // initial volume of the specimen
  
    REAL strainThreshold;	// for granular strain calculation, used to control tessellate.

    std::vector<cell*> cellVec;		// used to store cell info, for finite granular strain calculation, March 26, 2013
    std::vector<cell*> cellVec_init;	// the initial cell vector, for lagrangian finite strain

}; 


#endif
