/* $Id: NoseHooverT.cpp,v 1.8 2004-07-15 08:29:54 paklein Exp $ */
#include "NoseHooverT.h"
#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"
#include "RaggedArray2DT.h"
#include "BasicSupportT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
NoseHooverT::NoseHooverT(const BasicSupportT& support):
	ThermostatBaseT(support),
	fBetaOrig(0.0),
	fEta(0.0),
	fEtaDot(0.0)
{
	SetName("Nose-Hoover");
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
	
	int nsd = fSupport.NumSD();
	double dt = fSupport.TimeStep();
	
	/* calculate current temperature */
	double kineticTemp = 0.;
	int nDOF = 0;
	if (fAllNodes)
	{ // All the nodes are damped, use neighbors
		nDOF = nsd*neighbors.MajorDim();
		for (int j = 0; j < neighbors.MajorDim(); j++) 
		{
			int tag_j = *neighbors(j);
	    	const double* v_j = (*velocities)(tag_j);
				
			for (int i = 0; i < nsd; i++, *v_j++)
				kineticTemp += (*v_j)*(*v_j);
		}
	}
	else if (fNodes.Length() > 0)
	{
		nDOF = nsd*fNodes.Length();
		for (int j = 0; j < fNodes.Length(); j++)
		{ 
			int tag_j = fNodes[j];
			const double* v_j = (*velocities)(tag_j);

			for (int i = 0; i < nsd; i++, *v_j++)
				kineticTemp += (*v_j)*(*v_j); 	
	    }
	}
	kineticTemp /= fkB*nDOF;
	fEtaDot = (kineticTemp-fTemperature)/fBetaOrig;
	fBeta = fEta += fEtaDot*dt;
	
//	cout << " NoseHooverT::ApplyDamping KE = " << kineticTemp << " " << fBeta << "\n";
	
	ThermostatBaseT::ApplyDamping(neighbors,velocities,forces,
						types,particleProperties);
}			
	
