/* $Id: ThermostatBaseT.cpp,v 1.1 2003-04-16 18:15:54 cjkimme Exp $ */
#include "ThermostatBaseT.h"
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
ThermostatBaseT::ThermostatBaseT(ifstreamT& in, int nsd, double dt):
	fNodes(),
	fTemperature(-1.),
	fSD(nsd),
	fTimeStep(dt)
{
	in >> fBeta;
//	if (QLangevin)
//	{
//		in >> fTemperature;
//		fAmp = sqrt(2.*fBeta*fkB*fTemperature/
//						fTimeStep);
//	}
}

/* write properties to output */
void ThermostatBaseT::Write(ostream& out) const
{
	out << " Beta. . . . . . . . . . . . . . . . . . . . . . = " << fBeta << '\n';
//	if (QLangevin)
//		out << " Temperature . . . . . . . . . . . . . . . . . . = " << fTemperature << '\n';
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
		for (int j = 0; inBox && j < fSD; i++)
		{
			inBox = (xmin[j] < *x_i) && (*x_i++ < xmax[j]);
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
		case ThermostatBaseT::kGaussianIsokinetic:
		{
			property = ThermostatBaseT::kGaussianIsokinetic;
			break;
		}
		default:
			ExceptionT::BadInputValue("operator>>ThermostatBaseT::ThermostatT", 
				"unknown code: %d", i_property);
	}
	return in;
}

} /* namespace Tahoe */
