/* $Id: GaussIsokineticT.cpp,v 1.3 2003-04-24 20:43:20 cjkimme Exp $ */
#include "GaussIsokineticT.h"
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
GaussIsokineticT::GaussIsokineticT(ifstreamT& in, const int& nsd, const double& dt):
	ThermostatBaseT(in, nsd, dt)
{
	
}

/* write properties to output */
void GaussIsokineticT::Write(ostream& out) const
{
	ThermostatBaseT::Write(out);	
}

/* restart files */
void GaussIsokineticT::WriteRestart(ostream& out) const
{
	/* Base class */
	ThermostatBaseT::WriteRestart(out);
}

void GaussIsokineticT::ReadRestart(istream& in) 
{
	/* Base class */
	ThermostatBaseT::ReadRestart(in);
}

void GaussIsokineticT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
	/* Get the temperature */
	fTemperature = fTemperatureSchedule->Value()*fTemperatureScale;
	if (fTemperature < 0.)
		ExceptionT::GeneralFail("LangevinT::ApplyDamping","schedule generated negative temperature");

	double denom = 0.;
	double num = 0.;
	double* v_j;
	double* f_j;
	int tag_j, currType;
	double mass;
	
	/* calculate drag coefficient */
	if (fNodes.Length() == 0)
	{ // All the nodes are damped, use neighbors
		currType = types[*neighbors(0)];
		mass = particleProperties[currType]->Mass();
		for (int j = 0; j < neighbors.MajorDim(); j++) 
		{
			tag_j = *neighbors(j);
	    	v_j = (*velocities)(tag_j);
			f_j = forces(j);
			if (types[tag_j] != currType)
			{
				currType = types[tag_j];
				mass = particleProperties[currType]->Mass();
			}
				
			for (int i = 0; i < fSD; i++)
			{
				denom += (*v_j)*(*v_j);
				num += mass*(*f_j++)*(*v_j++);
			}
		}
	}
	else
	{
		currType = types[fNodes[0]];
		mass = particleProperties[currType]->Mass();
		for (int j = 0; j < fNodes.Length(); j++)
		{ 
			tag_j = fNodes[j];
			v_j = (*velocities)(tag_j);
			f_j = forces(j);
			if (types[tag_j] != currType)
			{
				currType = types[tag_j];
				mass = particleProperties[currType]->Mass();
			}
			
			for (int i = 0; i < fSD; i++)
			{
				denom += (*v_j)*(*v_j); 	
				num += mass*(*f_j++)*(*v_j++); 
			}
	    }
	}
	
	/* compute damping coefficient */
	if (abs(denom) > kSmall)
		fBeta = num/denom;
	else
		fBeta = 0.; 
	
	ThermostatBaseT::ApplyDamping(neighbors,velocities,forces,
							types,particleProperties);
}			
	
