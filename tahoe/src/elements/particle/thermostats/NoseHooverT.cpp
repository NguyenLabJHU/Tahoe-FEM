/* $Id: NoseHooverT.cpp,v 1.3 2003-04-22 01:23:16 cjkimme Exp $ */
#include "NoseHooverT.h"
#include "ArrayT.h"
#include <iostream.h>
#include "ifstreamT.h"
#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"
#include "RaggedArray2DT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
NoseHooverT::NoseHooverT(ifstreamT& in, const int& nsd, const double& dt):
	ThermostatBaseT(in, nsd, dt)
{
	fBetaOrig = fBeta;
	in >> fTemperature; 
	if (fTemperature < 0.)
		ExceptionT::BadInputValue("NoseHooverT::NoseHooverT","Negative Absolute T");
}

/* write properties to output */
void NoseHooverT::Write(ostream& out) const
{
	ThermostatBaseT::Write(out);	
}

/* restart files */
void NoseHooverT::WriteRestart(ostream& out) const
{
	/* Base class */
	ThermostatBaseT::WriteRestart(out);
	
	out << fBeta;
	out << fEta;
	out << fEtaDot;
}

void NoseHooverT::ReadRestart(istream& in) 
{
	/* Base class */
	ThermostatBaseT::ReadRestart(in);
	
	in >> fBeta;
	in >> fEta;
	in >> fEtaDot;
}

void NoseHooverT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
	/* calculate current temperature */
	double kineticTemp = 0.;
	int nDOF = 0;
	if (fNodes.Length() == 0)
	{ // All the nodes are damped, use neighbors
		nDOF = fSD*neighbors.MajorDim();
		for (int j = 0; j < neighbors.MajorDim(); j++) 
		{
			int tag_j = *neighbors(j);
	    	double* v_j = (*velocities)(tag_j);
				
			for (int i = 0; i < fSD; i++)
				kineticTemp += (*v_j++)*(*v_j);
		}
	}
	else
	{
		nDOF = fSD*fNodes.Length();
		for (int j = 0; j < fNodes.Length(); j++)
		{ 
			int tag_j = fNodes[j];
			double* v_j = (*velocities)(tag_j);

			for (int i = 0; i < fSD; i++)
				kineticTemp += (*v_j++)*(*v_j); 	
	    }
	}
	kineticTemp /= fkB*nDOF;
	fEtaDot = (kineticTemp-fTemperature)/fBetaOrig;
	fBeta = fEta += fEtaDot*fTimeStep;
	
	ThermostatBaseT::ApplyDamping(neighbors,velocities,forces,
						types,particleProperties);
}			
	
