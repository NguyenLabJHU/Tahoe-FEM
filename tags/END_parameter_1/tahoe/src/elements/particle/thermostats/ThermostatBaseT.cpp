/* $Id: ThermostatBaseT.cpp,v 1.4 2003-04-24 20:43:20 cjkimme Exp $ */
#include "ThermostatBaseT.h"
#include "ArrayT.h"
#include <iostream.h>
#include "ifstreamT.h"
//#include <math.h>
#include "dArrayT.h"
#include "dArray2DT.h"
#include "AutoArrayT.h"
#include "RaggedArray2DT.h"
#include "ParticlePropertyT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
ThermostatBaseT::ThermostatBaseT(ifstreamT& in, const int& nsd, 
	const double& dt):
	fNodes(),
	fTemperature(-1.),
	fSD(nsd),
	fTimeStep(dt),
	fTemperatureSchedule(NULL)
{
	in >> fBeta;
}

/* write properties to output */
void ThermostatBaseT::Write(ostream& out) const
{
	out << " Beta. . . . . . . . . . . . . . . . . . . . . . = " << fBeta << '\n';
}

/* restart files */
void ThermostatBaseT::WriteRestart(ostream& out) const
{
#pragma unused(out)
	// Do nothing
}

void ThermostatBaseT::ReadRestart(istream& in) 
{
#pragma unused(in)
	// Do nothing
}

void ThermostatBaseT::ApplyDamping(const RaggedArray2DT<int>& neighbors, const dArray2DT* velocities,
			dArray2DT& forces, AutoArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
	double* v_j;
	double* f_j;
	int tag_j, currType;
	double mass, beta;
	if (fNodes.Length() == 0)
	{ // All the nodes are damped, use neighbors
		currType = types[*neighbors(0)];
		mass = particleProperties[currType]->Mass();
		beta = fBeta*mass;
		for (int j = 0; j < neighbors.MajorDim(); j++) 
		{
			tag_j = *neighbors(j);
			f_j = forces(j);
	    	v_j = (*velocities)(tag_j);
	    	if (types[tag_j] != currType)
			{
				currType = types[tag_j];
				mass = particleProperties[currType]->Mass();
				beta = fBeta*mass;
			}
				
			for (int i = 0; i < fSD; i++)
				*f_j++ -= beta*(*v_j++);
		}
	}
	else
	{
		currType = types[fNodes[0]];
		mass = particleProperties[currType]->Mass();
		beta = fBeta*mass;
		for (int j = 0; j < fNodes.Length(); j++)
		{ 
			tag_j = fNodes[j];
			f_j = forces(j);
			v_j = (*velocities)(tag_j);
			if (types[tag_j] != currType)
			{
				currType = types[tag_j];
				mass = particleProperties[currType]->Mass();
				beta = fBeta*mass;
			}

			for (int i = 0; i < fSD; i++)
				*f_j++ -= beta*(*v_j++); 	
	    }
	}
}			

void ThermostatBaseT::InitRegion(ifstreamT& in, const dArray2DT& coords,	
					const ArrayT<int>* partition_nodes)
{
	fxmin.Dimension(fSD);
	fxmax.Dimension(fSD);

	in >> fxmin;
	in >> fxmax;
	in >> nIncs; 
	if (nIncs < 0) 
		ExceptionT::BadInputValue("ThermostatBaseT::InitRegion","Bad increment value");

	for (int i = 0; i < fSD; i++)
		if (fxmin[i] >= fxmax[i])
			ExceptionT::BadInputValue("ThermostatBaseT::InitRegion","Bad bounding box coordinates");
	
	/* get the nodes in the region */
	NodesInRegion(coords, partition_nodes);
}	

void ThermostatBaseT::NodesInRegion(const dArray2DT& coords,	
					const ArrayT<int>* partition_nodes)
{
	if (fxmin.Length() != coords.MinorDim())
		ExceptionT::GeneralFail("ThermostattedRegionT::NodesInRegion",
				"Dimension mismatch between coords and bounding box");
	AutoArrayT<int> tmpList;

	double* xmin = fxmin.Pointer();
	double* xmax = fxmax.Pointer(); 
	double* x_i;
	int ihits = 0;
	bool isSerial = !partition_nodes;
	int nnd = isSerial ? coords.MajorDim() : partition_nodes->Length();
	for (int i = 0; i < nnd; i++)
	{
		bool inBox = true;
		if (isSerial)
			x_i = coords(i);
		else
			x_i = coords((*partition_nodes)[i]);
		for (int j = 0; inBox && j < fSD; j++)
		{
			inBox = (xmin[j] < x_i[j]) && (x_i[j] < xmax[j]);
		}
		if (inBox)
		{
			tmpList.Append(i);
			ihits++;
		}
	}
	fNodes.Dimension(ihits);
	tmpList.CopyInto(fNodes);
}

void ThermostatBaseT::SetTemperatureSchedule(const ScheduleT* schedule, const double& value)
{
	fTemperatureSchedule = schedule;
	fTemperatureScale = value;
}	

namespace Tahoe {

/* stream extraction operator */
istream& operator>>(istream& in, ThermostatBaseT::ThermostatT& property)
{
	int i_property;
	in >> i_property;
	switch (i_property)
	{
		case ThermostatBaseT::kFree:
		{
			property = ThermostatBaseT::kFree;
			break;
		}
		case ThermostatBaseT::kDamped:
		{
			property = ThermostatBaseT::kDamped;
			break;
		}
		case ThermostatBaseT::kLangevin:
		{
			property = ThermostatBaseT::kLangevin;
			break;
		}
		case ThermostatBaseT::kNoseHoover:
		{
			property = ThermostatBaseT::kNoseHoover;
			break;
		}
		case ThermostatBaseT::kGaussIsokinetic:
		{
			property = ThermostatBaseT::kGaussIsokinetic;
			break;
		}
		case ThermostatBaseT::kRampedDamping:
		{
			property = ThermostatBaseT::kRampedDamping;
			break;
		}
		default:
			ExceptionT::BadInputValue("operator>>ThermostatBaseT::ThermostatT", 
				"unknown code: %d", i_property);
	}
	return in;
}

} /* namespace Tahoe */
