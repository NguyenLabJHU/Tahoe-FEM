/* $Id: NPT_EnsembleT.cpp,v 1.1.2.1 2003-09-18 21:03:38 cjkimme Exp $ */
#include "NPT_EnsembleT.h"
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

/* constructor */
NPT_EnsembleT::NPT_EnsembleT(ifstreamT& in, const ElementSupportT& support,
	const dArray2DT& dynStress):
	fElementSupport(support),
	ThermostatBaseT(in, support.NumSD(), support.TimeStep()),
	fDynStress(dynStress)
{
	
}

/* write properties to output */
void NPT_EnsembleT::Write(ostream& out) const
{
	ThermostatBaseT::Write(out);	
}

/* restart files */
void NPT_EnsembleT::WriteRestart(ostream& out) const
{
	/* Base class */
	ThermostatBaseT::WriteRestart(out);
}

void NPT_EnsembleT::ReadRestart(istream& in) 
{
	/* Base class */
	ThermostatBaseT::ReadRestart(in);
}

void NPT_EnsembleT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
	/* Get the temperature */
	double pressure = fTemperatureSchedule->Value()*fTemperatureScale;

	double denom = 0.;
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
				denom += mass*(*v_j)*(*v_j);
				num += (*f_j++)*(*v_j++);
			}
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
				denom += mass*(*v_j)*(*v_j); 	
				num += (*f_j++)*(*v_j++); 
			}
	    }
	}
	
	/* compute damping coefficient */
	if (fabs(denom) > kSmall)
		fBeta = -num/denom;
	else
		fBeta = 0.; 

//	cout <<" temp = "<< denom/natoms/fSD/fkB<<"\n";
	/* compute change in volume */

	
	ThermostatBaseT::ApplyDamping(neighbors,velocities,forces,
							types,particleProperties);
}			
	
