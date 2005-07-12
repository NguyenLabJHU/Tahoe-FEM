/* $Id: LangevinT.cpp,v 1.7.20.1 2005-06-05 06:23:55 paklein Exp $ */
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
LangevinT::LangevinT(const BasicSupportT& support):
	ThermostatBaseT(support)
{
	SetName("Langevin");
}

void LangevinT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, ArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
	int nsd = fSupport.NumSD();
	double dt = fSupport.TimeStep();

	/* quick exit */
	if (fabs(dt) < kSmall) return;
	
	dArrayT rArray(nsd); // random force
	double* rf_i;
	
	fTemperature = fTemperatureSchedule->Value()*fTemperatureScale;
	if (fTemperature < 0.)
		ExceptionT::GeneralFail("LangevinT::ApplyDamping","schedule generated negative temperature");
	double amp_m = sqrt(2.*fBeta*fkB*fTemperature/dt);

	if (fAllNodes)
	{ // All the nodes are damped, use neighbors
		int currType = types[*neighbors(0)];
		double mass = particleProperties[currType]->Mass();
		double amp = amp_m*sqrt(mass);
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
				amp = amp_m*sqrt(mass);
			}
				
			fRandom.RandomArray(rArray);
			rArray *= amp;
			rf_i = rArray.Pointer();
	    				
			for (int i = 0; i < nsd; i++)
				*f_j++ -= beta*(*v_j++) - *rf_i++;
		}
	}
	else if (fNodes.Length() > 0)
	{
		int currType = types[fNodes[0]];
		double mass = particleProperties[currType]->Mass();
		double amp = amp_m*sqrt(mass);
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
				amp = amp_m*sqrt(mass);
			}

			fRandom.RandomArray(rArray).Pointer();
			rArray *= amp;
			rf_i = rArray.Pointer();
	    				
			for (int i = 0; i < nsd; i++)
				*f_j++ -= beta*(*v_j++) - *rf_i++;	
	    }
	}
}			

/* describe the parameters needed by the interface */
void LangevinT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ThermostatBaseT::DefineParameters(list);
	
	/* random number seed */
	list.AddParameter(ParameterT::Integer, "random_seed");
}

/* accept parameter list */
void LangevinT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ThermostatBaseT::TakeParameterList(list);

	/* set the seed */
	int seed = list.GetParameter("random_seed");
	fRandom.sRand(seed);
}
