/* $Id: NoseHooverT.cpp,v 1.1 2003-04-16 18:15:54 cjkimme Exp $ */
#include "NoseHooverT.h"
#include "ArrayT.h"
#include <iostream.h>
#include "ifstreamT.h"
#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
NoseHooverT::NoseHooverT(ifstreamT& in, int nsd, double dt):
	ThermostatBaseT(in, nsd, dt)
{
	// not yet	
}

/* write properties to output */
void NoseHooverT::Write(ostream& out) const
{
	ThermostatBaseT::Write(out);	
}

/* restart files */
void NoseHooverT::WriteRestart(ostream& out) const
{
#pragma unused(out)
	// Do nothing
}

void NoseHooverT::ReadRestart(istream& in) 
{
#pragma unused(in)
	// Do nothing
}

void NoseHooverT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
					dArray2DT& forces)
{
	if (fNodes.Length() == 0)
	{ // All the nodes are damped, use neighbors
		for (int j = 0; j < neighbors.MajorDim(); j++) 
		{
			int tag_j = *neighbors(j);
			double* f_j = forces(j);
	    	double* v_j = (*velocities)(tag_j);
				
			for (int i = 0; i < fSD; i++)
				*f_j++ -= fBeta*(*v_j++);
		}
	}
	else
	{
		for (int j = 0; j < fNodes.Length(); j++)
		{ 
			int tag_j = fNodes[j];
			double* f_j = forces(j);
			double* v_j = (*velocities)(tag_j);

			for (int i = 0; i < fSD; i++)
				*f_j++ -= fBeta*(*v_j++); 	
	    }
	}
}			
	
