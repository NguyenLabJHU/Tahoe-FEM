/* $Id: ConfigurationalT.cpp,v 1.1.2.1 2003-09-18 21:03:38 cjkimme Exp $ */
#include "ConfigurationalT.h"
#include "ArrayT.h"
#include <iostream.h>
#include "ifstreamT.h"
#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"
#include "RaggedArray2DT.h"
#include "ElementSupportT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
ConfigurationalT::ConfigurationalT(ifstreamT& in, const ElementSupportT& support, 
		const dArrayT& delDotF):
	ThermostatBaseT(in, support.NumSD(), support.TimeStep()),
	fEta(0.),
	fDelDotF(delDotF)
{
	fBetaOrig = fBeta;
}

/* write properties to output */
void ConfigurationalT::Write(ostream& out) const
{
	ThermostatBaseT::Write(out);	
}

/* restart files */
void ConfigurationalT::WriteRestart(ostream& out) const
{
	/* Base class */
	ThermostatBaseT::WriteRestart(out);
	
	out << fBeta;
	out << fEta;
	out << fEtaDot;
}

void ConfigurationalT::ReadRestart(istream& in) 
{
	/* Base class */
	ThermostatBaseT::ReadRestart(in);
	
	in >> fBeta;
	in >> fEta;
	in >> fEtaDot;
}

void ConfigurationalT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
	/* Get temperature */
	fTemperature = fTemperatureSchedule->Value()*fTemperatureScale;
	if (fTemperature < 0.)
		ExceptionT::GeneralFail("LangevinT::ApplyDamping","schedule generated negative temperature");
	
	/* calculate current temperature */
	double configTemp = 0.;
	double configTempDenom = 0.; // divergence of force term
	double configTempNum = 0.; // force dot force term
	if (fNodes.Length() == 0)
	{ // All the nodes are damped, use neighbors
		for (int j = 0; j < neighbors.MajorDim(); j++) 
		{
			int tag_j = *neighbors(j);
			double* f_j = forces(tag_j);

//			configTempDenom += fDelDotF[tag_j];
				
			for (int i = 0; i < fSD; i++, *f_j++)
				configTempNum += (*f_j)*(*f_j);
		}
	}
	else
	{
		for (int j = 0; j < fNodes.Length(); j++)
		{ 
			int tag_j = fNodes[j];
			double* f_j = forces(tag_j);

//			configTempDenom += fDelDotF[tag_j];

			for (int i = 0; i < fSD; i++, *f_j++)
				configTempNum += (*f_j)*(*f_j); 	
	    }
	}

	/* now for the temperature */
	configTemp = configTempNum/configTempDenom/fkB;

	fEtaDot = (configTemp-fTemperature)/fBetaOrig;
	fBeta = fEta += fEtaDot*fTimeStep;
	
	cout << " ConfigurationalT::ApplyDamping KE = " << configTemp << " " << fBeta << "\n";
	
	// Can't apply damping here. It has to modify integration of positions
	/*	ThermostatBaseT::ApplyDamping(neighbors,velocities,forces,
		types,particleProperties);*/
}			
	
