/* $Id: RampedDampingT.cpp,v 1.4 2003-11-21 22:47:11 paklein Exp $ */
#include "RampedDampingT.h"
#include "ArrayT.h"
#include <iostream.h>
#include "ifstreamT.h"
//#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"
#include "AutoArrayT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
RampedDampingT::RampedDampingT(ifstreamT& in, const int& nsd, const double& dt):
	ThermostatBaseT(in,nsd,dt)
{
	SetName("ramped_damping");
	// Not yet
}

RampedDampingT::RampedDampingT(void)
{
	SetName("ramped_damping");
}

/* write properties to output */
void RampedDampingT::Write(ostream& out) const
{
	ThermostatBaseT::Write(out);
}

/* restart files */
void RampedDampingT::WriteRestart(ostream& out) const
{
	/* Base class */
	ThermostatBaseT::WriteRestart(out);
}

void RampedDampingT::ReadRestart(istream& in) 
{
	/* Base class */
	ThermostatBaseT::ReadRestart(in);
}

void RampedDampingT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
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

			for (int i = 0; i < fSD; i++)
				*f_j++ -= fBeta*(*v_j++); 	
	    }
	}
}				




