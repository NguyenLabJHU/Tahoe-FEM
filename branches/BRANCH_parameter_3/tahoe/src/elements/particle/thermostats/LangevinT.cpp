/* $Id: LangevinT.cpp,v 1.6.20.1 2004-05-25 16:36:43 paklein Exp $ */
#include "LangevinT.h"

#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"
#include "ParticlePropertyT.h"
#include "RaggedArray2DT.h"
#include "BasicSupportT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
#if 0
LangevinT::LangevinT(ifstreamT& in, const int& nsd, const double& dt):
	ThermostatBaseT(in,nsd,dt)
{
	SetName("Langevin");
//	in >> fTemperature;
//	fAmp = sqrt(2.*fBeta*fkB*fTemperature/fTimeStep);
}
#endif

LangevinT::LangevinT(const BasicSupportT& support)
{
	SetName("Langevin");
}

void LangevinT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
	int nsd = fSupport.NumSD();
	dArrayT rArray(nsd); // random force
	double* rf_i;
	
	fTemperature = fTemperatureSchedule->Value()*fTemperatureScale;
	if (fTemperature < 0.)
		ExceptionT::GeneralFail("LangevinT::ApplyDamping","schedule generated negative temperature");
	double amp = sqrt(2.*fBeta*fkB*fTemperature/fTimeStep);

	if (fAllNodes)



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
