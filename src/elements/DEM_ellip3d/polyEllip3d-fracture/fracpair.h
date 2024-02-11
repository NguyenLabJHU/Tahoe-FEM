#ifndef FRACPAIR_H
#define FRACPAIR_H

#include "realtypes.h"
#include "vec.h"
#include "particle.h"
#include "parameter.h"

// fracture pair store the two particles that are attached to the same common plane with 7 springs.
// And calculate the spring based on compression and tension displacement and to apply these spring resultants.
// October 17, 2013

// the four points of the springs are located in the four end points. For example, if the break plane is plane-ab,
// then point 1 is aminus on the principle-a axle, point 2 is bminus on the principle-b axle, point 3 is aplus on
// the principle-a axle and point 4 is bplus on the principle-b axle
//        4
//  	  |
//    1-------3
//	  |
//        2

namespace dem{

class fracpair{

    private:
	particle * p1;	// particle 1, the positive sub-particle
	particle * p2;	// particle 2, the negative sub-particle
	int type;
	vec local1_points[4];	// the three spring points on particle 1
	vec local2_points[4];	// the three spring points on particle 2, the other four points can be generated based on the three
	REAL k0_spring[4];	// the initial stiffness of the springs, k0=EA/nsl, we do not distinguish normal or tangential
//	REAL d[4];		// breakage parameter of each spring
	REAL l0_spring[4]; 	// the initial spring vector w.r.t. p1p2
   	dem::vec fc_init[4];	// the initial cohesive force in each spring
				// force in the springs are in the same direction as acting on particle p1, 
				// i.e., when calculate spring forces on p1, then just need to add the spring forces
				// while, when calculate spring forces on p2, then need to add the repulsant forces -fc
				// that's to say, force acting on p1, f1/2 (force acting on p1 by p2) is fc while f2/1 is -fc
	REAL zr;		// particle size, used to calculate the time used by crack to propagate through the particle	
	REAL area;		// the area of this fracture face
	int num_springs;	// the number of non-broken springs
	int num_broken;	// the number of springs broken 
	bool isInitialForce;	// is applied initial force or not, false: not, true: yes 
        bool isLive[4];	// spring is still live or not

//	REAL shear_strength;
//	REAL normal_strength;
//	REAL strength;	// not distinguish shear and normal
//	REAL stiffness;

    public:
	fracpair();
	fracpair(particle* t1, particle* t2, int break_plane);	// the break_plane means which is the break plane
//	~fracpair();	

	particle * getP1() const {return p1;}
	particle * getP2() const {return p2;}
	int getNumSprings() const {return num_springs;}
	int getNumBroken() const {return num_broken;}
	bool isIn(particle * t1);	// jugde if particle t1 is in this fracture pair based on the value of address
	bool getIsInitialForce() const {return isInitialForce;}	// return isInitialForce
	void setIsInitialForceTrue() {isInitialForce=true;}	// set isInitialForce = true
	void calcInitialCohesiveForce();	// calculate and apply initial cohesive force on p1 & p2
	void calculateResultant(REAL&);	// 1, calculate k_spring for those non-zero k based on displacement, compression and tension,
					// of the corresponding points. If compression or tension reaches threshold, then k=0 and 
					// num_springs -= 1
					// 2, calculate resultant forces and moments for each particle based on k and displacement


};


}// end of namespace dem

#endif


