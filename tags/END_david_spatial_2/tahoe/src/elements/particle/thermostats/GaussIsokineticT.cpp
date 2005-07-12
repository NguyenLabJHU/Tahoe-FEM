/* $Id: GaussIsokineticT.cpp,v 1.9.20.1 2005-06-05 06:23:55 paklein Exp $ */
#include "GaussIsokineticT.h"

#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"
#include "RaggedArray2DT.h"
#include "ParticlePropertyT.h"
#include "BasicSupportT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
GaussIsokineticT::GaussIsokineticT(const BasicSupportT& support):
	ThermostatBaseT(support)
{
	SetName("Gauss_isokinetic");
}

void GaussIsokineticT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, ArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
	/* Get the temperature */
	fTemperature = fTemperatureSchedule->Value()*fTemperatureScale;
	if (fTemperature < 0.)
		ExceptionT::GeneralFail("LangevinT::ApplyDamping","schedule generated negative temperature");

	int nsd = fSupport.NumSD();
	double denom = 0.;
	double num = 0.;
	const double* v_j;
	double* f_j;
	int tag_j, currType, natoms;
	double mass;
	
	/* calculate drag coefficient */
	if (fAllNodes)
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
				
			for (int i = 0; i < nsd; i++)
			{
				denom += mass*(*v_j)*(*v_j);
				num += (*f_j++)*(*v_j++);
			}
		}
	}
	else if (fNodes.Length() > 0)
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
		
			for (int i = 0; i < nsd; i++)
			{
				denom += mass*(*v_j)*(*v_j); 	
				num += (*f_j++)*(*v_j++); 
			}
	    }
	}
	
	/* compute damping coefficient */
	if (fabs(denom) > kSmall)
		fBeta = num/denom;
	else
		fBeta = 0.; 

//	cout <<" temp = "<< denom/natoms/fSD/fkB<<"\n";
	
	ThermostatBaseT::ApplyDamping(neighbors,velocities,forces,
							types,particleProperties);
}			
	
