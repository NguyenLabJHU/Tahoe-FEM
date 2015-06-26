 
#include "assembly.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cassert>
#include <utility>
#include <cmath>
#include <stdlib.h>
#include <math.h>

using std::cout;
using std::setw;
using std::endl;
using std::flush;
using std::vector;


// This is used to calculate the granular strain measures based on
// the experimental data from Professor Alshibli's group.
// Here we first calculate the granular strain measures for all the 
// domain and we calculate spatial dvdx and eulerian strain by rate form
// directly in this code. After that we will use matlab code to plot 
// eulerian strain by rate form, Lagrangian/Eulerian strain by rate-form 
// deformation gradient and Hencky strain.
// September 30, 2014 By Boning Zhang

void assembly::calculateGranularStrain(int num_step, int total_num, const char* iniptclfile, const char* progressfile,
				       REAL RVE_xmin, REAL RVE_xmax, REAL RVE_ymin, REAL RVE_ymax, REAL RVE_zmin, REAL RVE_zmax){

    nstep = num_step;
    TotalNum = total_num;

    int OPREC = 9;
    int OWID = 20;
    std::ofstream progressinf;
    progressinf.open(progressfile);
    if(!progressinf) { cout << "stream error!" << endl; exit(-1);}
    progressinf.setf(std::ios::scientific, std::ios::floatfield);
    progressinf.precision(OPREC);
    progressinf << setw(OWID) << "spatial" 	// for spatial velocity gradient tensor, deformation rate tensor
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "spatial" 
		<< setw(OWID) << "from rate" 	// for strain based on spatial deformation rate tensor, July 8, 2013
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate" 
		<< setw(OWID) << "from rate"
		<< setw(OWID) << "particle" << endl
	        << setw(OWID) << "dvdx_11"	// for spatial velocity gradient tensor, deformation rate tensor
	        << setw(OWID) << "dvdx_12"
	        << setw(OWID) << "dvdx_13"
	        << setw(OWID) << "dvdx_21"
	        << setw(OWID) << "dvdx_22"
	        << setw(OWID) << "dvdx_23"
	        << setw(OWID) << "dvdx_31"
	        << setw(OWID) << "dvdx_32"
	        << setw(OWID) << "dvdx_33"
	        << setw(OWID) << "epsilon_11"	// for strain based on spatial deformation rate tensor
	        << setw(OWID) << "epsilon_12"
	        << setw(OWID) << "epsilon_13"
	        << setw(OWID) << "epsilon_21"
	        << setw(OWID) << "epsilon_22"
	        << setw(OWID) << "epsilon_23"
	        << setw(OWID) << "epsilon_31"
	        << setw(OWID) << "epsilon_32"
	        << setw(OWID) << "epsilon_33"
		<< setw(OWID) << "number" << std::endl;


    matrix spatial_dvdx(3,3);	// average spatial velocity gradient tensor
    matrix curr_rate(3,3);	// current spatial deformation rate tensor, July 8, 2013
    matrix curr_dvdx(3,3);	// current spatial velocity gradient tensor
    matrix prev_strain_rate(3,3);
    matrix curr_strain_rate(3,3);	// current strain based on deformation rate tensor


    // initialize previousStrain(3,3)
    for(int i_ps=0; i_ps!=2; i_ps++){
	for(int j_ps=0; j_ps!=2; j_ps++){
		curr_strain_rate(i_ps+1,j_ps+1) = 0;	// calculate strain by rate, July 8, 2013
	}
    }
    
    readSample(iniptclfile);	// read particle information


/*
// find the xmin, xmax, ymin, ymax, zmin, zmax for all the particles
REAL xmin = (*ParticleVec.begin())->getCurrPosition().getx(); REAL xmax = (*ParticleVec.begin())->getCurrPosition().getx();
REAL ymin = (*ParticleVec.begin())->getCurrPosition().gety(); REAL ymax = (*ParticleVec.begin())->getCurrPosition().gety();
REAL zmin = (*ParticleVec.begin())->getCurrPosition().getz(); REAL zmax = (*ParticleVec.begin())->getCurrPosition().getz();
for(std::vector<particle*>::const_iterator pt=ParticleVec.begin(); pt!=ParticleVec.end(); pt++){
    if(xmin > (*pt)->getCurrPosition().getx()){
	xmin = (*pt)->getCurrPosition().getx();
    }
    if(xmax < (*pt)->getCurrPosition().getx()){
	xmax = (*pt)->getCurrPosition().getx();
    }

    if(ymin > (*pt)->getCurrPosition().gety()){
	ymin = (*pt)->getCurrPosition().gety();
    }
    if(ymax < (*pt)->getCurrPosition().gety()){
	ymax = (*pt)->getCurrPosition().gety();
    }

    if(zmin > (*pt)->getCurrPosition().getz()){
	zmin = (*pt)->getCurrPosition().getz();
    }
    if(zmax < (*pt)->getCurrPosition().getz()){
	zmax = (*pt)->getCurrPosition().getz();
    }
}
std::cout << "xmin: " << xmin << "  xmax: " << xmax << std::endl;
std::cout << "ymin: " << ymin << "  ymax: " << ymax << std::endl;
std::cout << "zmin: " << zmin << "  zmax: " << zmax << std::endl;
*/

    // find particles in RVE, store the numbers of these particles in NOvec_part
    std::vector<particle*>::const_iterator  it_p;
    RVE_number = 0;	// the number of particles in RVE
    REAL px, py, pz;
    ParticleVec_part.clear();
    for (it_p=ParticleVec.begin();it_p!=ParticleVec.end();++it_p)  {	
	px = (*it_p)->getCurrPosition().getx();
	py = (*it_p)->getCurrPosition().gety();
	pz = (*it_p)->getCurrPosition().getz();
	if(px>=RVE_xmin && px<=RVE_xmax && py>=RVE_ymin && py<=RVE_ymax && pz>=RVE_zmin && pz<=RVE_zmax){	// means this particle is in RVE
		ParticleVec_part.push_back(*it_p);
		RVE_number++;
	}
    
    }

    if(RVE_number<5){	// Qhull needs more than 4 particles
	progressinf << "\n the number of particles in this RVE becomes: " << RVE_number << std::endl;
	return;
    }

    // first tessellation
    createInputForQhull();
    callQhull();
//    readTesse("tess_info");
    readTesse_finite("tess_info");
    cellVec_init = cellVec;	// keep cellVec_init for lagrangian strain unchanged 

    REAL* time_interval = new REAL[nstep];
    time_interval[0] = 60;	// 1 minutes
    time_interval[1] = 60;
    time_interval[2] = 90;
    time_interval[3] = 90;
    time_interval[4] = 60;
    time_interval[5] = 120;
    time_interval[6] = 120;
    time_interval[7] = 180;
//    time_interval[8] = 330;



    // calculate granular strain measures
    for(int i=0; i<nstep; i++){

std::cout << "begin step " << i+1 << "!" << std::endl;

	updateParticle(i, time_interval[i]);

//	resetStartPosition();	// reset initial position

    	RVE_number = 0;	// the number of particles in RVE
    	ParticleVec_part.clear();
    	for (it_p=ParticleVec.begin();it_p!=ParticleVec.end();++it_p)  {	
	    px = (*it_p)->getCurrPosition().getx();
	    py = (*it_p)->getCurrPosition().gety();
	    pz = (*it_p)->getCurrPosition().getz();
	    if(px>=RVE_xmin && px<=RVE_xmax && py>=RVE_ymin && py<=RVE_ymax && pz>=RVE_zmin && pz<=RVE_zmax){	// means this particle is in RVE
		ParticleVec_part.push_back(*it_p);
		RVE_number++;
	    }
    
    	}

      	if(RVE_number<5){	// Qhull needs more than 4 particles
	    progressinf << "\n the number of particles in this RVE becomes: " << RVE_number << std::endl;
	    return;
    	}

	// tessellate again
	createInputForQhull();
	callQhull();
	readTesse_finite("tess_info");

	// calculate strain by rate
	prev_strain_rate.clear();
	prev_strain_rate = curr_strain_rate;	// set previous strain
	
	// calculate spatial velocity gradient tensor
	spatial_dvdx.clear();
	spatial_dvdx = getAverage_dvdx();
	curr_dvdx.clear();
	curr_dvdx = spatial_dvdx;
	
	curr_rate = 0.5*(curr_dvdx+curr_dvdx.getTrans());	// current spatial deformation rate
	curr_strain_rate = prev_strain_rate+(curr_rate-curr_dvdx.getTrans()*prev_strain_rate-prev_strain_rate*curr_dvdx)*time_interval[i];

	progressinf << setw(OWID) << spatial_dvdx(1,1) << setw(OWID) << spatial_dvdx(1,2) << setw(OWID) << spatial_dvdx(1,3)
		    << setw(OWID) << spatial_dvdx(2,1) << setw(OWID) << spatial_dvdx(2,2) << setw(OWID) << spatial_dvdx(2,3)
		    << setw(OWID) << spatial_dvdx(3,1) << setw(OWID) << spatial_dvdx(3,2) << setw(OWID) << spatial_dvdx(3,3)
		    << setw(OWID) << curr_strain_rate(1,1) << setw(OWID) << curr_strain_rate(1,2) << setw(OWID) << curr_strain_rate(1,3)
		    << setw(OWID) << curr_strain_rate(2,1) << setw(OWID) << curr_strain_rate(2,2) << setw(OWID) << curr_strain_rate(2,3)
		    << setw(OWID) << curr_strain_rate(3,1) << setw(OWID) << curr_strain_rate(3,2) << setw(OWID) << curr_strain_rate(3,3)
		    << setw(OWID) << RVE_number
	            << endl;

    } // end for

    progressinf.close();


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
	ParticleVec_part.clear();
    	cellVec.clear();
    	cellVec_init.clear();

} // end calculateGranularStrain


// create input file for Qhull
void assembly::createInputForQhull() const{

        int OPREC = 9;
        int OWID = 20;

	std::ofstream ofs("input_for_Qhull");
	if(!ofs){
		cout << "Stream error when create input for Qhull!" << endl;
		exit (-1);
	}
	ofs.setf(std::ios::scientific, std::ios::floatfield);
  	ofs.precision(OPREC);
	ofs << setw(OWID) << "3" << endl
	    << setw(OWID) << RVE_number << endl;
	vec tmp;
  	std::vector<particle*>::const_iterator  it;
  	for (it=ParticleVec_part.begin();it!=ParticleVec_part.end();++it)  {
    		tmp=(*it)->getCurrPosition();
    		ofs << setw(OWID) << tmp.getx()
		    << setw(OWID) << tmp.gety()
		    << setw(OWID) << tmp.getz() << endl;
	}

    	ofs.close();
} // end createInputForQhull


// call Qhull to tessellate and ouput "tess_info"
void assembly::callQhull() const{
	std::cout << "Checking if processor is available..." << std::endl;
	if(system(NULL)) std::cout << "Ok!" << std::endl;
	else exit(EXIT_FAILURE);
	system("./qdelaunay Qt i < input_for_Qhull TO tess_info");	// call the external command qdelaunay
} // end callQhull


// read cell information into std::vector<cell> cellVec from Qhull output file, for finite granular strain calculation
// written on March 26, 2013
void assembly::readTesse_finite(const char* str){
 	REAL length_scale = 1.0e-9;	// microns in input file

	std::ifstream ifs(str);
	std::cout << "Read cell information begin!" << std::endl;
	if(!ifs) {
		std::cout << "stream error!" << std::endl; exit(-1);
	}
	int m, n, i, j;
	int totalNum;

	ifs>>totalNum;
	cellVec.clear();
	for(int it=0; it!=totalNum; it++){
		ifs>>m>>n>>i>>j;	// read nodes in each line
		m = m+1;
		n = n+1;
		i = i+1;
		j = j+1;	// the ID from Qhull is starting from 0
		cell* ptCell = new cell();
		ptCell->setNodes(m,n,i,j);
		cellVec.push_back(ptCell);
	}
	setNumberingOrder();	// this is important
} // end readTesse_finite


// finite granular strain calculation, March 26, 2013
void assembly::setNumberingOrder() {
	int nm, nn, ni, nj;	// 4 nodes of a tetrahedron
	int n1, n2 , n3, n4;	// the final odering node 
	std::vector<particle*>::const_iterator it;
//	REAL x1, y1, z1;	// the current position of node 1
//	REAL x2, y2, z2;	// the current position of node 2
//	REAL x3, y3, z3;	// the current position of node 3
//	REAL x4, y4, z4;	// the current position of node 4
	REAL cell_volume;	// area of triangle, volume of tet
	for(std::vector<cell*>::const_iterator iter=cellVec.begin(); iter!=cellVec.end(); iter++){
		nm = (*iter)->getm();
		nn = (*iter)->getn();
		ni = (*iter)->geti();
		nj = (*iter)->getj();
		// 1, choose m, n as node 1, 2
		n1 = nm;
		n2 = nn;
		n3 = ni;
		n4 = nj;
		(*iter)->setNodes(n1,n2,n3,n4);
		cell_volume = getCellVolume(**iter);

		if(cell_volume < 0){	// swap n2 and n3
			n2 = ni;
			n3 = nn;
			(*iter)->setNodes(n1,n2,n3,n4);
			cell_volume = getCellVolume(**iter);
		}
		
//		// get coordinates for node 1, 2
//		it = ParticleVec.begin();
//		for(int i=0; i!=n1-1; i++)
//			it++;	// go to the n1 particle
//		x1 = ((*it)->getCurrPosition()).getx();
//		y1 = ((*it)->getCurrPosition()).gety();
//		z1 = ((*it)->getCurrPosition()).getz();
//
//		it = ParticleVec.begin();
//		for(int i=0; i!=n2-1; i++)
//			it++;	// go to the n2 particle
//		x2 = ((*it)->getCurrPosition()).getx();
//		y2 = ((*it)->getCurrPosition()).gety();
//		z2 = ((*it)->getCurrPosition()).getz();
//		
//		// 2, choose i as node 3 and test if Area n1,n2,n3 is positive
//		n3 = ni;
//		n4 = nj;
//		// get coordinates for node 3
//		it = ParticleVec.begin();
//		for(int i=0; i!=n3-1; i++)
//			it++;	// go to the n3 particle
//		x3 = ((*it)->getCurrPosition()).getx();
//		y3 = ((*it)->getCurrPosition()).gety();
//		z3 = ((*it)->getCurrPosition()).getz();
//		
//		// calculate area
//		area = 0.5*( x1*y2+x3*y1+x2*y3-x3*y2-x1*y3-x2*y1);
//		if(area < 0){	// means the previous three nodes are not numbered counter-clockwise
//			n3 = nj;
//			n4 = ni;
//			// get coordinates for node 3
//			it = ParticleVec.begin();
//			for(int i=0; i!=n3-1; i++)
//				it++;	// go to the n3 particle
//			x3 = ((*it)->getCurrPosition()).getx();
//			y3 = ((*it)->getCurrPosition()).gety();
//			z3 = ((*it)->getCurrPosition()).getz();
//		}
//		// get coordinates for node 4
//		it = ParticleVec.begin();
//		for(int i=0; i!=n4-1; i++)
//			it++;	// go to the n3 particle
//		x4 = ((*it)->getCurrPosition()).getx();
//		y4 = ((*it)->getCurrPosition()).gety();
//		z4 = ((*it)->getCurrPosition()).getz();
//		// test if volume is positive
//		cell_volume = x1*y2*z3+x2*y3*z4+x3*y4*z1+x4*y1*z2
//			     -x1*y4*z3-x2*y1*z4-x3*y2*z1-x4*y3*z2;
		if(cell_volume < 0){	// means n1,n2,n3,n4 are not numbered counter-clockwise
			std::cout << "Error: nodes are not numbered counter-clockwise!" << std::endl;
		}

		(*iter)->setInitialCellVolume(cell_volume);
		matrix bigB(4,3);
		bigB = getBigB(**iter);
		(*iter)->setInitialBigB(bigB);
	}
} // end setNumberingOrder


// caclulate average spatial velocity gradient tensor, June 24, 2013
matrix assembly::getAverage_dvdx() const{
	matrix average_dvdx(3,3);
	REAL curr_totalVolume = 0;	// the current summed volume of all cells
	// initialize average_dvdx(3,3)
	for(int ir=1; ir!=4;ir++){
		for(int ic=1; ic!=4; ic++){
			average_dvdx(ir,ic) = 0;
		}
	}
	matrix temp_dvdx(3,3);
	REAL temp_volume;
	for(std::vector<cell*>::const_iterator iter=cellVec.begin(); iter!=cellVec.end(); iter++){
		temp_volume = getCellVolume(**iter);
		temp_dvdx = getdvdx_curr(**iter);

		average_dvdx += temp_volume*temp_dvdx;	
		curr_totalVolume += temp_volume;
	}
	average_dvdx = average_dvdx/curr_totalVolume;	
	return average_dvdx;
} // end getAverage_dvdx


// get volume of cell
REAL assembly::getCellVolume(const cell& tempCell) const{
	int n1,n2,n3,n4;
	REAL x1, y1, z1;	// the current position of node 1
	REAL x2, y2, z2;	// the current position of node 2
	REAL x3, y3, z3;	// the current position of node 3
	REAL x4, y4, z4;	// the current position of node 4
	// get ID for 4 nodes
	n1 = tempCell.getm();
	n2 = tempCell.getn();
	n3 = tempCell.geti();
	n4 = tempCell.getj();
	
	// get coordinates for the 4 nodes
	std::vector<particle*>::const_iterator it;
	
	it = ParticleVec_part.begin();
	for(int i=0; i!=n1-1; i++)
		it++;	// go to the n1 particle
	x1 = ((*it)->getCurrPosition()).getx();
	y1 = ((*it)->getCurrPosition()).gety();
	z1 = ((*it)->getCurrPosition()).getz();
	
	it = ParticleVec_part.begin();
	for(int i=0; i!=n2-1; i++)
		it++;	// go to the n2 particle
	x2 = ((*it)->getCurrPosition()).getx();
	y2 = ((*it)->getCurrPosition()).gety();
	z2 = ((*it)->getCurrPosition()).getz();

	it = ParticleVec_part.begin();
	for(int i=0; i!=n3-1; i++)
		it++;	// go to the n3 particle
	x3 = ((*it)->getCurrPosition()).getx();
	y3 = ((*it)->getCurrPosition()).gety();
	z3 = ((*it)->getCurrPosition()).getz();

	it = ParticleVec_part.begin();
	for(int i=0; i!=n4-1; i++)
		it++;	// go to the n4 particle
	x4 = ((*it)->getCurrPosition()).getx();
	y4 = ((*it)->getCurrPosition()).gety();
	z4 = ((*it)->getCurrPosition()).getz();

	// get volume of this cell
	REAL cell_volume;
	cell_volume = (x2*y3*z4-x2*y4*z3-x3*y2*z4+x3*y4*z2+x4*y2*z3-x4*y3*z2
		      -x1*y3*z4+x1*y4*z3+x3*y1*z4-x3*y4*z1-x4*y1*z3+x4*y3*z1
		      +x1*y2*z4-x1*y4*z2-x2*y1*z4+x2*y4*z1+x4*y1*z2-x4*y2*z1
		      -x1*y2*z3+x1*y3*z2+x2*y1*z3-x2*y3*z1-x3*y1*z2+x3*y2*z1)/6.0;
	return cell_volume;
} // end getCellVolume


// calculate spatial velocity gradient tensor for each tet, June 24, 2013
matrix assembly::getdvdx_curr(const cell& tempCell) const{
	int n1,n2,n3,n4;

	// get ID for 4 nodes
	n1 = tempCell.getm();
	n2 = tempCell.getn();
	n3 = tempCell.geti();
	n4 = tempCell.getj();
	// get velocities for 4 nodes
	matrix v1(3,1), v2(3,1), v3(3,1), v4(3,1);
	std::vector<particle*>::const_iterator it;
	// go to particle n1
	it = ParticleVec_part.begin();
	for(int i=0; i!=n1-1; i++)
		it++;	// go to the n1 particle
	v1(1,1) = ((*it)->getCurrVelocity()).getx();
	v1(2,1) = ((*it)->getCurrVelocity()).gety();
	v1(3,1) = ((*it)->getCurrVelocity()).getz();
	// go to particle n2
	it = ParticleVec_part.begin();
	for(int i=0; i!=n2-1; i++)
		it++;	// go to the n2 particle
	v2(1,1) = ((*it)->getCurrVelocity()).getx();
	v2(2,1) = ((*it)->getCurrVelocity()).gety();
	v2(3,1) = ((*it)->getCurrVelocity()).getz();
	// go to particle n3
	it = ParticleVec_part.begin();
	for(int i=0; i!=n3-1; i++)
		it++;	// go to the n3 particle
	v3(1,1) = ((*it)->getCurrVelocity()).getx();
	v3(2,1) = ((*it)->getCurrVelocity()).gety();
	v3(3,1) = ((*it)->getCurrVelocity()).getz();
	// go to particle n4
	it = ParticleVec_part.begin();
	for(int i=0; i!=n4-1; i++)
		it++;	// go to the n4 particle
	v4(1,1) = ((*it)->getCurrVelocity()).getx();
	v4(2,1) = ((*it)->getCurrVelocity()).gety();
	v4(3,1) = ((*it)->getCurrVelocity()).getz();

	// get volume of this cell
	REAL cell_volume;
	cell_volume = getCellVolume(tempCell);
	if(cell_volume < 0){
		std::cout << "Error: volume is negative in getdvdx_curr()!" << std::endl;
	}

	matrix bigV(3,4), bigB(4,3);
	// assemble bigU
	bigV(1,1) = v1(1,1);
	bigV(2,1) = v1(2,1);
	bigV(3,1) = v1(3,1);

	bigV(1,2) = v2(1,1);
	bigV(2,2) = v2(2,1);
	bigV(3,2) = v2(3,1);

	bigV(1,3) = v3(1,1);
	bigV(2,3) = v3(2,1);
	bigV(3,3) = v3(3,1);

	bigV(1,4) = v4(1,1);
	bigV(2,4) = v4(2,1);
	bigV(3,4) = v4(3,1);

	bigB = getBigB(tempCell);
	matrix result;
	result = 1/(6.0*cell_volume)*bigV*bigB;
	return result;
} // end getdvdx_curr


matrix assembly::getBigB(const cell& tempCell) const{
	int n1,n2,n3,n4;
	REAL x1, y1, z1;	// the current position of node 1
	REAL x2, y2, z2;	// the current position of node 2
	REAL x3, y3, z3;	// the current position of node 3
	REAL x4, y4, z4;	// the current position of node 4
	// get ID for 4 nodes
	n1 = tempCell.getm();
	n2 = tempCell.getn();
	n3 = tempCell.geti();
	n4 = tempCell.getj();

	// get coordinates for the 4 nodes
	std::vector<particle*>::const_iterator it;
	
	it = ParticleVec_part.begin();
	for(int i=0; i!=n1-1; i++)
		it++;	// go to the n1 particle
	x1 = ((*it)->getCurrPosition()).getx();
	y1 = ((*it)->getCurrPosition()).gety();
	z1 = ((*it)->getCurrPosition()).getz();
	
	it = ParticleVec_part.begin();
	for(int i=0; i!=n2-1; i++)
		it++;	// go to the n2 particle
	x2 = ((*it)->getCurrPosition()).getx();
	y2 = ((*it)->getCurrPosition()).gety();
	z2 = ((*it)->getCurrPosition()).getz();

	it = ParticleVec_part.begin();
	for(int i=0; i!=n3-1; i++)
		it++;	// go to the n3 particle
	x3 = ((*it)->getCurrPosition()).getx();
	y3 = ((*it)->getCurrPosition()).gety();
	z3 = ((*it)->getCurrPosition()).getz();

	it = ParticleVec_part.begin();
	for(int i=0; i!=n4-1; i++)
		it++;	// go to the n4 particle
	x4 = ((*it)->getCurrPosition()).getx();
	y4 = ((*it)->getCurrPosition()).gety();
	z4 = ((*it)->getCurrPosition()).getz();

	// get intermediate variables
	REAL a1, a2, a3, a4;
	REAL b1, b2, b3, b4;
	REAL c1, c2, c3, c4;
	a1 = y2*(z4-z3)-y3*(z4-z2)+y4*(z3-z2);
	a2 = -y1*(z4-z3)+y3*(z4-z1)-y4*(z3-z1);
	a3 = y1*(z4-z2)-y2*(z4-z1)+y4*(z2-z1);
	a4 = -y1*(z3-z2)+y2*(z3-z1)-y3*(z2-z1);

	b1 = -x2*(z4-z3)+x3*(z4-z2)-x4*(z3-z2);
	b2 = x1*(z4-z3)-x3*(z4-z1)+x4*(z3-z1);
	b3 = -x1*(z4-z2)+x2*(z4-z1)-x4*(z2-z1);
	b4 = x1*(z3-z2)-x2*(z3-z1)+x3*(z2-z1);

	c1 = x2*(y4-y3)-x3*(y4-y2)+x4*(y3-y2);
	c2 = -x1*(y4-y3)+x3*(y4-y1)-x4*(y3-y1);
	c3 = x1*(y4-y2)-x2*(y4-y1)+x4*(y2-y1);
	c4 = -x1*(y3-y2)+x2*(y3-y1)-x3*(y2-y1);

	// assemble bigB
	matrix bigB(4,3);
	bigB(1,1) = a1;
	bigB(2,1) = a2;
	bigB(3,1) = a3;
	bigB(4,1) = a4;

	bigB(1,2) = b1;
	bigB(2,2) = b2;
	bigB(3,2) = b3;
	bigB(4,2) = b4;
	
	bigB(1,3) = c1;
	bigB(2,3) = c2;
	bigB(3,3) = c3;
	bigB(4,3) = c4;
	
	return bigB;
} // end getBigB


void assembly::updateParticle(int step, REAL time){
    for(std::vector<particle*>::iterator it=ParticleVec.begin();it!=ParticleVec.end();++it){
	(*it)->update(step, time);
    }
}


void assembly::readSample(const char* str){
    std::ifstream ifs(str);
    if(!ifs) {
	cout << "stream error!" << endl; exit(-1);
    }
    REAL temp;
    for(int i=0; i<nstep; i++){
    	ifs >> temp;
    }


    ParticleVec.clear();

    REAL a, b, c, px,py,pz,dax,day,daz,dbx,dby,dbz,dcx,dcy,dcz;
    REAL vx,vy,vz,omx,omy,omz,fx,fy,fz,mx,my,mz;
    for (int i=0;i<TotalNum;i++){
	ifs >> temp;

	REAL* posi_x = new REAL[nstep];
        REAL* posi_y = new REAL[nstep];
    	REAL* posi_z = new REAL[nstep];
	for(int j=0; j<nstep; j++){
	    ifs >> posi_x[j] >> posi_y[j] >> posi_z[j];
	}	
	
	particle* pt= new particle(i+1,posi_x,posi_y,posi_z);
	ParticleVec.push_back(pt);
    }
    ifs.close();
} // end readSample

 
