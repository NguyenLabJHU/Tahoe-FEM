/* $Id: ConstrainedPressureT.cpp,v 1.1.2.1 2003-09-18 21:03:38 cjkimme Exp $ */
#include "ConstrainedPressureT.h"
#include "ArrayT.h"
#include <iostream.h>
#include "ifstreamT.h"
#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"
#include "RaggedArray2DT.h"
#include "ParticlePropertyT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* initialize static variables */
double ConstrainedPressureT::staticChi = 0.;

/* constructor */
ConstrainedPressureT::ConstrainedPressureT(ifstreamT& in, const ElementSupportT& support,
	const dArray2DT& dynStress):
	fElementSupport(support),
	ThermostatBaseT(in, support.NumSD(), support.TimeStep()),
	fDynStress(dynStress)
{
	
}

/* write properties to output */
void ConstrainedPressureT::Write(ostream& out) const
{
	ThermostatBaseT::Write(out);	
}

/* restart files */
void ConstrainedPressureT::WriteRestart(ostream& out) const
{
	/* Base class */
	ThermostatBaseT::WriteRestart(out);
}

void ConstrainedPressureT::ReadRestart(istream& in) 
{
	/* Base class */
	ThermostatBaseT::ReadRestart(in);
}

void ConstrainedPressureT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
	/* Get the pressure */
	double pressure = fTemperatureSchedule->Value()*fTemperatureScale;

	double denom = 9.*pressure*fV;
	double num = 0.;
	double* v_j;
	double* f_j;
	int tag_j, currType, natoms;
	double mass;
	
	/* calculate drag coefficient */
	if (fNodes.Length() == 0)
	{ // All the nodes are damped, use neighbors
		currType = types[*neighbors(0)];
		mass = particleProperties[currType]->Mass();
		natoms = neighbors.MajorDim();
		for (int j = 0; j < natoms; j++) 
		{
			tag_j = *neighbors(j);
	    	v_j = (*velocities)(tag_j);
	 		f_j = forces(tag_j);
			if (types[tag_j] != currType)
			{
				currType = types[tag_j];
				mass = particleProperties[currType]->Mass();
			}
				
			for (int i = 0; i < fSD; i++)
			{
				denom += 2.*mass*(*v_j)*(*v_j);
				num += 2.*(*f_j++)*(*v_j++);
			}
			
			denom += fDynStress(tag_j,1);
			num -= fDynStress(tag_j,0);
		}
	}
	else
	{
		currType = types[fNodes[0]];
		mass = particleProperties[currType]->Mass();
		natoms = fNodes.Length();
		for (int j = 0; j < natoms; j++)
		{ 
			tag_j = fNodes[j];
			v_j = (*velocities)(tag_j);
			f_j = forces(tag_j);
			if (types[tag_j] != currType)
			{
				currType = types[tag_j];
				mass = particleProperties[currType]->Mass();
			}
		
			for (int i = 0; i < fSD; i++)
			{
				denom += 2.*mass*(*v_j)*(*v_j); 	
				num += 2.*(*f_j++)*(*v_j++); 
			}
			
			denom += fDynStress(tag_j,1);
			num -= fDynStress(tag_j,0);
	    }
	}
	
	/* compute damping coefficient */
	if (fabs(denom) > kSmall)
		fBeta = num/denom;
	else
		fBeta = 0.;
		
	/* set parameter for velocity corrector */
	staticChi = fBeta;
		
	fVdot = 3*fV*fBeta;
	
	/* get scale factor to multiply periodic lengths by */
	double vScale = 1 + fVdot/fV;
	fV += fVdot*fTimeStep;
	// Set new PBCs	 


//	cout <<" temp = "<< denom/natoms/fSD/fkB<<"\n";

	
	ThermostatBaseT::ApplyDamping(neighbors,velocities,forces,
							types,particleProperties);
}

/*
 *  static functions below needed to modify velocity integration. Called by
 *  a KBCController.
 *
 */

double ConstrainedPressureT::VelocityCorrector(void)
{
	return staticChi;
}
	
