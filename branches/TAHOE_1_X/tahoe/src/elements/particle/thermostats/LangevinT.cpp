/* $Id: LangevinT.cpp,v 1.6 2003-11-21 22:47:11 paklein Exp $ */
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
LangevinT::LangevinT(ifstreamT& in, const int& nsd, const double& dt):
	ThermostatBaseT(in,nsd,dt)
{
	SetName("Langevin");
//	in >> fTemperature;
//	fAmp = sqrt(2.*fBeta*fkB*fTemperature/fTimeStep);
}

LangevinT::LangevinT(void)
{
	SetName("Langevin");
}

/* write properties to output */
void LangevinT::Write(ostream& out) const
{
	ThermostatBaseT::Write(out);
	
	out << " Temperature . . . . . . . . . . . . . . . . . . = " << fTemperatureSchedule->Value()*fTemperatureScale << '\n';
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
	
	fTemperature = fTemperatureSchedule->Value()*fTemperatureScale;
	if (fTemperature < 0.)
		ExceptionT::GeneralFail("LangevinT::ApplyDamping","schedule generated negative temperature");
	fAmp = sqrt(2.*fBeta*fkB*fTemperature/fTimeStep);

	if (fNodes.Length() == 0)
	{ // All the nodes are damped, use neighbors
		int currType = types[*neighbors(0)];
		double mass = particleProperties[currType]->Mass();
		double amp = fAmp*sqrt(mass);
		double beta = fBeta*mass;
		for (int j = 0; j < neighbors.MajorDim(); j++) 
		{
			int tag_j = *neighbors(j);
			double* f_j = forces(tag_j);
	    	const double* v_j = (*velocities)(tag_j);
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
		int currType = types[fNodes[0]];
		double mass = particleProperties[currType]->Mass();
		double amp = fAmp*sqrt(mass);
		double beta = fBeta*mass;
		for (int j = 0; j < fNodes.Length(); j++)
		{ 
			int tag_j = fNodes[j];
			double* f_j = forces(tag_j);
			const double* v_j = (*velocities)(tag_j);
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
