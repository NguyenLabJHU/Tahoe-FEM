/* $Id: NoseHooverT.cpp,v 1.9.4.1 2005-06-05 06:23:55 paklein Exp $ */
#include "NoseHooverT.h"
#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"
#include "RaggedArray2DT.h"
#include "ParticlePropertyT.h"
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
			dArray2DT& forces, ArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
	/* Get temperature */
	fTemperature = fTemperatureSchedule->Value()*fTemperatureScale;
	if (fTemperature < 0.)
		ExceptionT::GeneralFail("LangevinT::ApplyDamping","schedule generated negative temperature");
	
	int nsd = fSupport.NumSD();
	double dt = fSupport.TimeStep();
	
	/* calculate current temperature */
	double kineticTemp = 0., mass;
	int nDOF = 0, currType;
	if (fAllNodes) { // All the nodes are damped, use neighbors
	  int currType = types[*neighbors(0)];
	  double mass = particleProperties[currType]->Mass();
	  nDOF = nsd*neighbors.MajorDim();
	  for (int j = 0; j < neighbors.MajorDim(); j++) {
	    int tag_j = *neighbors(j);
	    const double* v_j = (*velocities)(tag_j);
	    if (types[tag_j] != currType) {
	      currType = types[tag_j];
	      mass = particleProperties[currType]->Mass();
	    }

	    for (int i = 0; i < nsd; i++, *v_j++)
	      kineticTemp += mass*(*v_j)*(*v_j);
	  }
	}
	else if (fNodes.Length() > 0) {
	  currType = types[fNodes[0]];
	  mass = particleProperties[currType]->Mass();
	  nDOF = nsd*fNodes.Length();
	  for (int j = 0; j < fNodes.Length(); j++) { 
	    int tag_j = fNodes[j];
	    const double* v_j = (*velocities)(tag_j);

	    if (types[tag_j] != currType) {
	      currType = types[tag_j];
	      mass = particleProperties[currType]->Mass();
	    }

	    for (int i = 0; i < nsd; i++, *v_j++)
	      kineticTemp += mass*(*v_j)*(*v_j); 	
	  }
	}
	kineticTemp /= fkB*nDOF;
	fEtaDot = (kineticTemp-fTemperature)/fBetaOrig;
	fBeta = fEta += fEtaDot*dt;
	
//cout << " NoseHooverT::ApplyDamping KE = " << kineticTemp << " " << fBeta << " " << nDOF << "\n";
//cout << fEtaDot*fEtaDot/fEta/fEta/2. << " " <<  nDOF*fkB*fTemperature << " " << fEta << "\n";	
	ThermostatBaseT::ApplyDamping(neighbors,velocities,forces,
						types,particleProperties);
}			
	
