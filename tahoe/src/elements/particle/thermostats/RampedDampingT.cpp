/* $Id: RampedDampingT.cpp,v 1.5.30.1 2005-07-25 02:37:14 paklein Exp $ */
#include "RampedDampingT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "BasicSupportT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
RampedDampingT::RampedDampingT(const BasicSupportT& support):
	ThermostatBaseT(support),
	fBeta(0.0)
{
	SetName("ramped_damping");
}

void RampedDampingT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, ArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
	int nsd = fSupport.NumSD();

#pragma unused(neighbors)
#pragma unused(types)
#pragma unused(particleProperties)
	if (qNodesInRegion && fNodes.Length() == 0)
	{ 
		ExceptionT::GeneralFail("RampedDampingT::ApplyDamping","Where have all the nodes gone?");
	}
	else
	{
		for (int j = 0; j < fNodes.Length(); j++)
		{ 
			int tag_j = fNodes[j];
			double* f_j = forces(j);
			const double* v_j = (*velocities)(tag_j);

			for (int i = 0; i < nsd; i++)
				*f_j++ -= fBeta*(*v_j++); 	
	    }
	}
}				




