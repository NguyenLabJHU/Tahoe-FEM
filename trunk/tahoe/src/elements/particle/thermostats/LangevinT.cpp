/* $Id: LangevinT.cpp,v 1.1 2003-04-16 18:15:54 cjkimme Exp $ */
#include "LangevinT.h"
#include "ArrayT.h"
#include <iostream.h>
#include "ifstreamT.h"
#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"

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
#pragma unused(out)
	// Do nothing
}

void LangevinT::ReadRestart(istream& in) 
{
#pragma unused(in)
	// Do nothing
}

void LangevinT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
					dArray2DT& forces)
{

	dArrayT rArray(fSD); // random force
	double *rf_i;

	if (fNodes.Length() == 0)
	{ // All the nodes are damped, use neighbors
		for (int j = 0; j < neighbors.MajorDim(); j++) 
		{
			int tag_j = *neighbors(j);
			double* f_j = forces(j);
	    	double* v_j = (*velocities)(tag_j);
				
			fRandom->RandomArray(rArray);
			rArray *= fAmp;
			rf_i = rArray.Pointer();
	    				
			for (int i = 0; i < fSD; i++)
				*f_j++ -= fBeta*(*v_j++) - *rf_i++;
		}
	}
	else
	{
		for (int j = 0; j < fNodes.Length(); j++)
		{ 
			int tag_j = fNodes[j];
			double* f_j = forces(j);
			double* v_j = (*velocities)(tag_j);

			fRandom->RandomArray(rArray);
			rArray *= fAmp;
			rf_i = rArray.Pointer();
	    				
			for (int i = 0; i < fSD; i++)
				*f_j++ -= fBeta*(*v_j++) - *rf_i++;	
	    }
	}
}			
