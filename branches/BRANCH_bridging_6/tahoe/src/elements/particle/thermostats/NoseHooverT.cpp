/* $Id: NoseHooverT.cpp,v 1.7 2003-11-21 22:47:11 paklein Exp $ */
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
	ThermostatBaseT(in, nsd, dt),
	fEta(0.)
{
	SetName("NoseHoover");
	fBetaOrig = fBeta;
}

NoseHooverT::NoseHooverT(void)
{
	SetName("NoseHoover");
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

/* accept parameter list */
void NoseHooverT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ThermostatBaseT::TakeParameterList(list);
	fBetaOrig = fBeta;
}

void NoseHooverT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
	/* Get temperature */
	fTemperature = fTemperatureSchedule->Value()*fTemperatureScale;
	if (fTemperature < 0.)
		ExceptionT::GeneralFail("LangevinT::ApplyDamping","schedule generated negative temperature");
	
	/* calculate current temperature */
	double kineticTemp = 0.;
	int nDOF = 0;
	if (fNodes.Length() == 0)
	{ // All the nodes are damped, use neighbors
		nDOF = fSD*neighbors.MajorDim();
		for (int j = 0; j < neighbors.MajorDim(); j++) 
		{
			int tag_j = *neighbors(j);
	    	const double* v_j = (*velocities)(tag_j);
				
			for (int i = 0; i < fSD; i++, *v_j++)
				kineticTemp += (*v_j)*(*v_j);
		}
	}
	else
	{
		nDOF = fSD*fNodes.Length();
		for (int j = 0; j < fNodes.Length(); j++)
		{ 
			int tag_j = fNodes[j];
			const double* v_j = (*velocities)(tag_j);

			for (int i = 0; i < fSD; i++, *v_j++)
				kineticTemp += (*v_j)*(*v_j); 	
	    }
	}
	kineticTemp /= fkB*nDOF;
	fEtaDot = (kineticTemp-fTemperature)/fBetaOrig;
	fBeta = fEta += fEtaDot*fTimeStep;
	
//	cout << " NoseHooverT::ApplyDamping KE = " << kineticTemp << " " << fBeta << "\n";
	
	ThermostatBaseT::ApplyDamping(neighbors,velocities,forces,
						types,particleProperties);
}			
	
