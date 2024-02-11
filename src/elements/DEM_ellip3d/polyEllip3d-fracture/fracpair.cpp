#include "fracpair.h"
#include "matrix.h"
#include <iostream>
#include <cstdlib>

//#define K_MAX 1E6
//#define DIST_COMPRESSION_MAX -0.05
//#define DIST_TENSION_MAX 0.1

// October 17, 2013
namespace dem{


fracpair::fracpair(){
    p1 = NULL;
    p2 = NULL;
    type = 0;
    num_springs = 0;
//    normal_strength = 4;
//    shear_strength = 4;
//    strength = 0;	// not distinguish shear and normal at present
//    stiffness = 0;

    num_broken = 0;	// number of springs broken

    area = 0;

}

// create a fracture pair based on the different break planes
// by definition, t1 is the sub-particle that is in positive direction
// while t2 is the one in negative direction
fracpair::fracpair(particle* t1, particle* t2, int break_plane){
    p1 = t1;
    p2 = t2;
    type = break_plane;

    // calculate local points in p1 and p2
    REAL l;
    switch (type){
	case 1:	// break plane is ab-plane
	    local1_points[0] = vec(-p1->getAminus(), 0, -p1->getCminus());
	    local1_points[1] = vec(0, -p1->getBminus(), -p1->getCminus());
	    local1_points[2] = vec(p1->getAplus(), 0, -p1->getCminus());
	    local1_points[3] = vec(0, p1->getBplus(), -p1->getCminus());
	    
	    local2_points[0] = vec(-p2->getAminus(), 0, p2->getCplus());
	    local2_points[1] = vec(0, -p2->getBminus(), p2->getCplus());
	    local2_points[2] = vec(p2->getAplus(), 0, p2->getCplus());
	    local2_points[3] = vec(0, p2->getBplus(), p2->getCplus());

	    zr = 0.5*(p1->getAplus() + p1->getAminus() + p1->getBplus() + p1->getBminus());
	    // get spring stiffness
	    area = 0.25*PI*( p1->getAplus()*p1->getBplus() + p1->getAplus()*p1->getBminus() + p1->getAminus()*p1->getBplus() + p1->getAminus()*p1->getBminus() );
	    l = p1->getCminus() + p1->getCplus() + p2->getCminus() + p2->getCplus();
	    k0_spring[0] = p1->getYoung()*area/4.0/l;	// 4.0 means 4 spring points
//	    kt_spring[0] = 0.5*p1->getYoung()/(1+p1->getPoisson())*area/4.0/l;

	    break;
	case 2: // break plane is ac-plane
	    local1_points[0] = vec(-p1->getAminus(), -p1->getBminus(), 0);
	    local1_points[1] = vec(0, -p1->getBminus(), -p1->getCminus());
	    local1_points[2] = vec(p1->getAplus(), -p1->getBminus(), 0);
	    local1_points[3] = vec(0, -p1->getBminus(), p1->getCplus());

	    local2_points[0] = vec(-p2->getAminus(), p2->getBplus(), 0);
	    local2_points[1] = vec(0, p2->getBplus(), -p2->getCminus());
	    local2_points[2] = vec(p2->getAplus(), p2->getBplus(), 0);
	    local2_points[3] = vec(0, p2->getBplus(), p2->getCplus());

	    zr = 0.5*(p1->getAplus() + p1->getAminus() + p1->getCplus() + p1->getCminus());
	    // get spring stiffness
	    area = 0.25*PI*( p1->getAplus()*p1->getCplus() + p1->getAplus()*p1->getCminus() + p1->getAminus()*p1->getCplus() + p1->getAminus()*p1->getCminus() );
	    l = p1->getBminus() + p1->getBplus() + p2->getBminus() + p2->getBplus();
	    k0_spring[0] = p1->getYoung()*area*0.25/l;
//	    kt_spring[0] = 0.5*p1->getYoung()/(1+p1->getPoisson())*area/4.0/l;

	    break;
	case 3:	// break plane is bc-plane
	    local1_points[0] = vec(-p1->getAminus(), -p1->getBminus(), 0);
	    local1_points[1] = vec(-p1->getAminus(), 0, -p1->getCminus());
	    local1_points[2] = vec(-p1->getAminus(), p1->getBplus(), 0);
	    local1_points[3] = vec(-p1->getAminus(), 0, p1->getCplus());

	    local1_points[0] = vec(p2->getAplus(), -p2->getBminus(), 0);
	    local1_points[1] = vec(p2->getAplus(), 0, -p2->getCminus());
	    local1_points[2] = vec(p2->getAplus(), p2->getBplus(), 0);
	    local1_points[3] = vec(p2->getAplus(), 0, p2->getCplus());

	    zr = 0.5*(p1->getBplus() + p1->getBminus() + p1->getCplus() + p1->getCminus());
	    // get spring stiffness
	    area = 0.25*PI*( p1->getBplus()*p1->getCplus() + p1->getBplus()*p1->getCminus() + p1->getBminus()*p1->getCplus() + p1->getBminus()*p1->getCminus() );
	    l = p1->getAminus() + p1->getAplus() + p2->getAminus() + p2->getAplus();
	    k0_spring[0] = p1->getYoung()*area/4.0/l;
//	    kt_spring[0] = 0.5*p1->getYoung()/(1+p1->getPoisson())*area/4.0/l;

	    break;

	default:
	    std::cout << "The fracture type should be be 1/2/3!" << std::endl;
	    exit(-1);
    }// end of switch

//    strength = 2e5;
//    stiffness = 2e6;
//    stiffness = 2e7;

    // initialize k_springs
    k0_spring[1] = k0_spring[0];
    k0_spring[2] = k0_spring[0];
    k0_spring[3] = k0_spring[0];
//    kt_spring[1] = kt_spring[0];
//    kt_spring[2] = kt_spring[0];
//    kt_spring[3] = kt_spring[0];

    // initialize num_springs
    num_springs = 4;
//    d[0]=0; d[1]=0; d[2]=0; d[3]=0;
    fc_init[0]=0; fc_init[1]=0; fc_init[2]=0; fc_init[3]=0;	// the initial cohesive foce will be calculated next step
    isInitialForce = false;					// after calculate contact forces and before update, since 
    num_broken = 0;	// number of broken springs		// we need the separate forces acting on this part
    isLive[0]=true; isLive[1]=true;
    isLive[2]=true; isLive[3]=true;

    // get initial spring distance, since some numerical problem, these distance are not zero exactly.
    // then we need to store the initial values to be compared
    vec global1_points[4];
    vec global2_points[4];

    global1_points[0] = p1->globalVec(local1_points[0]) + p1->getCurrPosition();
    global1_points[1] = p1->globalVec(local1_points[1]) + p1->getCurrPosition();
    global1_points[2] = p1->globalVec(local1_points[2]) + p1->getCurrPosition();
    global1_points[3] = p1->globalVec(local1_points[3]) + p1->getCurrPosition();
//    global1_points[3] = 0.5*(global1_points[0]+global1_points[1]);
//    global1_points[4] = 0.5*(global1_points[1]+global1_points[2]);
//    global1_points[5] = 0.5*(global1_points[0]+global1_points[2]);
//    global1_points[6] = 1.0/3.0*(global1_points[0]+global1_points[1]+global1_points[2]);

    global2_points[0] = p2->globalVec(local2_points[0]) + p2->getCurrPosition();
    global2_points[1] = p2->globalVec(local2_points[1]) + p2->getCurrPosition();
    global2_points[2] = p2->globalVec(local2_points[2]) + p2->getCurrPosition();
    global2_points[3] = p2->globalVec(local2_points[3]) + p2->getCurrPosition();
//    global2_points[3] = 0.5*(global2_points[0]+global2_points[1]);
//    global2_points[4] = 0.5*(global2_points[1]+global2_points[2]);
//    global2_points[5] = 0.5*(global2_points[0]+global2_points[2]);
//    global2_points[6] = 1.0/3.0*(global2_points[0]+global2_points[1]+global2_points[2]);

    for(int i=0; i!=4; i++){
	l0_spring[i] = vfabs(global2_points[i] - global1_points[i]);
    }

}

// jugde if particle t1 belongs to this fracture pair
bool fracpair::isIn(particle * t1){
    if(t1==p1 || t1==p2)
	return true;
    else
	return false;
}

void fracpair::calcInitialCohesiveForce(){
    dem::vec fc1 = p1->calculateInitialCohesiveForce();	// total (by all four springs) initial cohesive force acting on p1 by p2 
    fc_init[0] = 0.25*fc1; fc_init[1] = 0.25*fc1;	// currently we don't consider the moment should applied by these
    fc_init[2] = 0.25*fc1; fc_init[2] = 0.25*fc1;	// four springs, these moment should be provided by the difference between forces
							// since we believe when the particle is broken, their rotation should not be very large
} // applyInitialCohesiveForce


// calculate k for each spring and add force/ moments to p1 and p2
void fracpair::calculateResultant(REAL &fracForce){
//    REAL max_radius = p1->getMaxRadius();	// as threshold for broken
//    if(p2->getMaxRadius()>max_radius)
//	max_radius = p2->getMaxRadius();

    vec global1_points[4];
    vec global2_points[4];

    global1_points[0] = p1->globalVec(local1_points[0]) + p1->getCurrPosition();
    global1_points[1] = p1->globalVec(local1_points[1]) + p1->getCurrPosition();
    global1_points[2] = p1->globalVec(local1_points[2]) + p1->getCurrPosition();
    global1_points[3] = p1->globalVec(local1_points[3]) + p1->getCurrPosition();
//    global1_points[3] = 0.5*(global1_points[0]+global1_points[1]);
//    global1_points[4] = 0.5*(global1_points[1]+global1_points[2]);
//    global1_points[5] = 0.5*(global1_points[0]+global1_points[2]);
//    global1_points[6] = 1.0/3.0*(global1_points[0]+global1_points[1]+global1_points[2]);

    global2_points[0] = p2->globalVec(local2_points[0]) + p2->getCurrPosition();
    global2_points[1] = p2->globalVec(local2_points[1]) + p2->getCurrPosition();
    global2_points[2] = p2->globalVec(local2_points[2]) + p2->getCurrPosition();
    global2_points[3] = p2->globalVec(local2_points[3]) + p2->getCurrPosition();
//    global2_points[3] = 0.5*(global2_points[0]+global2_points[1]);
//    global2_points[4] = 0.5*(global2_points[1]+global2_points[2]);
//    global2_points[5] = 0.5*(global2_points[0]+global2_points[2]);
//    global2_points[6] = 1.0/3.0*(global2_points[0]+global2_points[1]+global2_points[2]);

    vec l_spring;	// vector of spring, pointing from point 1 to point 2
    REAL delta_x;	// l-l0 of spring

    fracForce = 0;
    for(int i=0; i!=4; i++){
	// for non-broken springs
	if(isLive[i]){
	    // get normalized direction of springs
	    l_spring = (global2_points[i] - global1_points[i]);
	    delta_x = vfabs(l_spring) - l0_spring[i];

	    if(delta_x<=0){	// spring is in compression, no cohesive force
		continue;
	    }

	    REAL fc_norm = vfabs(fc_init[i])-k0_spring[i]*delta_x;	// cohesive force in spring

	    vec fc = fc_norm*normalize(fc_init[i]);	// pointing from point 1 to point 2
	    if(fc_norm<=0){
		isLive[i]=false;	// no spring force, break the spring
		num_broken++;
		num_springs--;
		continue;
	    }

	    p1->addForce(fc);			// don't need to consider these forces in the average stress, since 
	    p2->addForce(-fc);			// when there are springs for this particle, it can not be sub-divided further
	    p1->addMoment( (global1_points[i]-p1->getCurrCenterMass())* fc );
	    p2->addMoment( (global2_points[i]-p2->getCurrCenterMass())* (-fc) ); 

	    fracForce += vfabs(fc);

	}// end if
    }// end for

if(num_broken+num_springs != 4){
std::cout << "summation of non-broken and broken springs are not 4..." << std::endl;
}

}//end of calculateResultant()


}// end of namespace dem





