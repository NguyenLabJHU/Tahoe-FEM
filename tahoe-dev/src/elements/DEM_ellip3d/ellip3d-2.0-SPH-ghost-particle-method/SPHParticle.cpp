// Function Definitions
#include "SPHParticle.h"
#include "particle.h"

namespace sph {

    SPHParticle::SPHParticle(REAL m, REAL rho, REAL x, REAL y, REAL z, int t){
	mass = m;
    	density = rho;
    	volume = mass/density;
    	pressure = dem::P0*(std::pow(density/dem::SPHInitialDensity, dem::gamma)-1);;
    	mu = density*dem::nu;
    	densityDot = 0;

    	initial_X = dem::vec(x,y,z);
	curr_x = initial_X;
    	velocity = 0;
    	velocityDot = 0;

	local_X = 0;	// free SPH point
	
	if(t!=1 && t!=3){
	    std::cout << "Error in creating SPH free/boundary particle!" << std::endl;
	    exit(-1);
	}
	type = t;
	demParticle = NULL;
	
    } // end SPHParticle()

    SPHParticle::SPHParticle(REAL m, REAL rho, REAL x, REAL y, REAL z, dem::vec local, dem::particle* p, int t){
	mass = m;
    	density = rho;
    	volume = mass/density;
    	pressure = dem::P0*(std::pow(density/dem::SPHInitialDensity, dem::gamma)-1);;
    	mu = density*dem::nu;
    	densityDot = 0;

    	initial_X = dem::vec(x,y,z);
	curr_x = initial_X;
    	velocity = 0;
    	velocityDot = 0;

	local_X = local;	// ghost SPH point

	if(t!=2){
	    std::cout << "Error in creating ghost SPH particle!" << std::endl;
	    exit(-1);
	}
	type = 2;
	demParticle = p;
	
    } // end SPHParticle()
	
    void SPHParticle::initial(){

    	volume = mass/density;
    	pressure = dem::P0*(std::pow(density/dem::SPHInitialDensity, dem::gamma)-1);;
    	mu = density*dem::nu;
    	densityDot = 0;

	curr_x = initial_X;
    	velocity = 0;
    	velocityDot = 0;

    } // end initial()

    void SPHParticle::updateParticle(){	// here the forward Euler time integration is used

//if(fabs(densityDot)>100){
//    densityDot = 0.1*densityDot;
//}
	density = density+densityDot*(dem::TIMESTEP);
	velocity = velocity+velocityDot*(dem::TIMESTEP)-(dem::sphDamping)*velocity;
	curr_x = curr_x+(velocity+velocityCorrection)*(dem::TIMESTEP);
//	curr_x = curr_x+velocity*(dem::TIMESTEP);

    } // end updateParticle()

    void SPHParticle::updateParticleDensity(){	// here the forward Euler time integration is used

//if(fabs(densityDot)>100){
//    densityDot = 0.1*densityDot;
//}
	density = density+densityDot*(dem::TIMESTEP);

    } // end updateParticle()


    void SPHParticle::updateParticlePositionDensityLeapFrog(){	// update position and density based on equation (4.1)

	curr_x = curr_x+(velocity+velocityCorrection)*(dem::TIMESTEP);
	density = density+densityDot*(dem::TIMESTEP);

    } // end updateParticlePositionDensityLeapFrog()

    void SPHParticle::updateParticleVelocityLeapFrog(){	// update velocity based on equation (4.2)

	velocity = velocity+velocityDot*(dem::TIMESTEP);

    } // end updateParticleVelocityLeapFrog()

    void SPHParticle::initialParticleVelocityLeapFrog(){	// update velocity based on equation (4.3)

	velocity = velocity+velocityDot*(dem::TIMESTEP)*0.5;

    } // end initialParticleVelocityLeapFrog()



} // end namespace sph

