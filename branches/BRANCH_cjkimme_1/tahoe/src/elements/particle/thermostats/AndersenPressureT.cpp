/* $Id: AndersenPressureT.cpp,v 1.1.2.1 2003-09-18 21:03:38 cjkimme Exp $ */
#include "AndersenPressureT.h"
#include "ArrayT.h"
#include <iostream.h>
#include "ifstreamT.h"
#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"
#include "RaggedArray2DT.h"
#include "ParticlePropertyT.h"
#include "CommManagerT.h"

const double fkB = 0.00008617385;
const double fVOrder = 3;

using namespace Tahoe;

/* initialize static variables */
double AndersenPressureT::MdP = 0.;
dArrayT* AndersenPressureT::fV_field = NULL;

/* constructor */
AndersenPressureT::AndersenPressureT(ifstreamT& in, const ElementSupportT& support,
	const double* virial, const double boxLength):
	fElementSupport(support),
	ThermostatBaseT(in, support.NumSD(), support.TimeStep()),
	fVirial(virial),
	fV_local(fVOrder+1)
{	
	fV_local = 0.; // zero all the derivatives of the system volume
	fV_local[0] = pow(boxLength,fSD);
	fV_field = &fV_local; // set the static'n	
}

/* write properties to output */
void AndersenPressureT::Write(ostream& out) const
{
	ThermostatBaseT::Write(out);	
}

GlobalT::RelaxCodeT AndersenPressureT::RelaxSystem(void)
{
	double len = pow(fV_local[0],1./double(fSD));
	
	for (int i = 0; i < fSD; i++)
		fElementSupport.CommManager().SetPeriodicBoundaries(i,len,len);
		
	return GlobalT::kReEQ;
}

/* restart files */
void AndersenPressureT::WriteRestart(ostream& out) const
{
	/* Base class */
	ThermostatBaseT::WriteRestart(out);
}

void AndersenPressureT::ReadRestart(istream& in) 
{
	/* Base class */
	ThermostatBaseT::ReadRestart(in);
}

void AndersenPressureT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
#pragma unused(forces)
	/* Get the prescribed pressure */
	double pressure = fTemperatureSchedule->Value()*fTemperatureScale;

	double* v_j;
	int tag_j, currType, natoms;
	double mass;
	double KE = 0.;
	
	/* calculate current pressure */
	if (fNodes.Length() == 0)
	{ // All the nodes are damped, use neighbors
		currType = types[*neighbors(0)];
		mass = particleProperties[currType]->Mass();
		natoms = neighbors.MajorDim();
		for (int j = 0; j < natoms; j++) 
		{
			tag_j = *neighbors(j);
	    	v_j = (*velocities)(tag_j);

			if (types[tag_j] != currType)
			{
				currType = types[tag_j];
				mass = particleProperties[currType]->Mass();
			}
				
			for (int i = 0; i < fSD; i++, v_j++)
			{
				KE += mass*(*v_j)*(*v_j);
			}
		}
	}
	else
	{
		currType = types[fNodes[0]];
		mass = particleProperties[currType]->Mass();
		natoms = fNodes.Length();
		for (int j = 0; j < natoms; j++)
		{ 
			tag_j = fNodes[j];
			v_j = (*velocities)(tag_j);
			
			if (types[tag_j] != currType)
			{
				currType = types[tag_j];
				mass = particleProperties[currType]->Mass();
			}
		
			for (int i = 0; i < fSD; i++, v_j++)
			{
				KE += mass*(*v_j)*(*v_j);
			}
	    }
	}
	
	/* compute damping coefficient */
	KE += *fVirial;
	MdP = (KE-pressure)/fBeta;
	
	/* still have to get scale factor to multiply periodic lengths by */
	cout << "AndersentPressureT:: Volume is " << fV_local[0] <<" Mdp is "<< MdP <<"\n";

	/* No damping to apply. The integrator will handle everything. */
	/* Alternatively, this could handle the V terms in the corrector. What's best? */
	/*ThermostatBaseT::ApplyDamping(neighbors,velocities,forces,
						types,particleProperties);*/
}

/*
 *  static function below needed to modify integration. Called by
 *  Gear4 integrator.
 *
 */

double AndersenPressureT::VelocityCorrector(void)
{
	return MdP;
}
	
