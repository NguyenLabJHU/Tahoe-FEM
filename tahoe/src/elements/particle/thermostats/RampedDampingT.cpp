/* $Id: RampedDampingT.cpp,v 1.2.12.1 2003-09-18 21:03:38 cjkimme Exp $ */
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
	// Damping over a region with a variable damping coefficient 
	// Region is defined in ParticleT::EchoDamping after construction here
	
	// Direction it varies in 
	in >> nRampedDOF;
	if (nRampedDOF > nsd) ExceptionT::GeneralFail("RampedDampingT::RampedDampingT","damping DOF > num spatial DOFs");
	nRampedDOF--;
	
	// Values at 2 endpoints
	in >> fLLval >> fURval;
	fLLval *= fBeta;
	fURval *= fBeta;
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

	for (int j = 0; j < fNodes.Length(); j++)
	{ 
		int tag_j = fNodes[j];
		double* f_j = forces(j);
		double* v_j = (*velocities)(tag_j);

		for (int i = 0; i < fSD; i++)
			*f_j++ -= fBeta*(*v_j++); 	
    }

}				
