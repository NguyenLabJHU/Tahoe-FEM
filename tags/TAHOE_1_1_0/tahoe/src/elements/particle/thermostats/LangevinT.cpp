/* $Id: LangevinT.cpp,v 1.2 2003-04-18 19:01:56 cjkimme Exp $ */
#include "LangevinT.h"
#include "ArrayT.h"
#include <iostream.h>
#include "ifstreamT.h"
#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"
#include "ParticlePropertyT.h"
#include "RaggedArray2DT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
LangevinT::LangevinT(ifstreamT& in, int nsd, double dt):
	ThermostatBaseT(in,nsd,dt)
{
	in >> fTemperature;
	fAmp = sqrt(2.*fBeta*fkB*fTemperature/fTimeStep);
}

/* write properties to output */
void LangevinT::Write(ostream& out) const
{
	ThermostatBaseT::Write(out);
	
	out << " Temperature . . . . . . . . . . . . . . . . . . = " << fTemperature << '\n';
}

/* restart files */
void LangevinT::WriteRestart(ostream& out) const
{
	/* Base class */
	ThermostatBaseT::WriteRestart(out);
}

void LangevinT::ReadRestart(istream& in) 
{
	/* Base class */
	ThermostatBaseT::ReadRestart(in);
}

void LangevinT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{

	dArrayT rArray(fSD); // random force
	double* rf_i;

	if (fNodes.Length() == 0)
	{ // All the nodes are damped, use neighbors
		int currType = types[*neighbors(0)];
		double mass = particleProperties[currType]->Mass();
		double amp = fAmp*sqrt(mass);
		double beta = fBeta*mass;
		for (int j = 0; j < neighbors.MajorDim(); j++) 
		{
			int tag_j = *neighbors(j);
			double* f_j = forces(j);
	    	double* v_j = (*velocities)(tag_j);
			if (types[tag_j] != currType)
			{
				currType = types[tag_j];
				mass = particleProperties[currType]->Mass();
				beta = fBeta*mass;
				amp = fAmp*sqrt(mass);
			}
				
			fRandom->RandomArray(rArray);
			rArray *= amp;
			rf_i = rArray.Pointer();
	    				
			for (int i = 0; i < fSD; i++)
				*f_j++ -= beta*(*v_j++) - *rf_i++;
		}
	}
	else
	{
		double currType = types[fNodes[0]];
		double mass = particleProperties[currType]->Mass();
		double amp = fAmp*sqrt(mass);
		double beta = fBeta*mass;
		for (int j = 0; j < fNodes.Length(); j++)
		{ 
			int tag_j = fNodes[j];
			double* f_j = forces(j);
			double* v_j = (*velocities)(tag_j);
			if (types[tag_j] != currType)
			{
				currType = types[tag_j];
				mass = particleProperties[currType]->Mass();
				beta = fBeta*mass;
				amp = fAmp*sqrt(mass);
			}

			fRandom->RandomArray(rArray).Pointer();
			rArray *= amp;
			rf_i = rArray.Pointer();
	    				
			for (int i = 0; i < fSD; i++)
				*f_j++ -= beta*(*v_j++) - *rf_i++;	
	    }
	}
}			
