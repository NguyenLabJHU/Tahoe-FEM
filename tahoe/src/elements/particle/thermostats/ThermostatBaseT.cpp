/* $Id: ThermostatBaseT.cpp,v 1.11 2003-12-28 23:37:27 paklein Exp $ */
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
#include "ModelManagerT.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
ThermostatBaseT::ThermostatBaseT(ifstreamT& in, const int& nsd, 
	const double& dt):
	ParameterInterfaceT("thermostat"),
	fBeta(0.0),
	fTemperature(-1.),
	fTimeStep(dt),
	fSD(nsd),
	fTemperatureSchedule(NULL)
{
	in >> fBeta;
}

ThermostatBaseT::ThermostatBaseT(void):
	ParameterInterfaceT("thermostat"),
	fBeta(0.0),
	fTemperature(0.0),
	fSD(0),
	fTimeStep(0.0),
	fTemperatureSchedule(NULL)
{
	SetName("thermostat");
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
	const double* v_j;
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
			f_j = forces(tag_j);
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

void ThermostatBaseT::InitNodeSets(ifstreamT& in, ModelManagerT& model)
{
	/* Get node sets */
	char peekahead = in.next_char();
	int allOrSome;
	if (peekahead == '-')
		in >> allOrSome;
	else
		allOrSome = atoi(&peekahead);
	if (!allOrSome)
	{   // use all the nodes
		in >> allOrSome; // fast-forward the stream
	}
	else	
		if (allOrSome > 0)
		{   // Read in Node Sets and use them

            /* read node set ids */
            ArrayT<StringT> ids;
            model.NodeSetList(in, ids);
            iArrayT tags;
            model.ManyNodeSets(ids, fNodes);        
		}
		else
		{   // Read in Node Sets and use all nodes but theirs
			ArrayT<StringT> notTheseSets;
			int numNotSets = -allOrSome;
			numNotSets = abs(numNotSets);
			notTheseSets.Dimension(numNotSets); 
			
			for (int i=0; i < numNotSets; i++)
			{
	  			StringT& name = notTheseSets[i];
	  			in >> name;
	  			int index = model.NodeSetIndex(name);
	  			if (index < 0) {
	  				cout << "\n ScaledVelocityT::Initialize: error retrieving node set " << name << endl;
	  				throw ExceptionT::kDatabaseFail;
	  			}
			}
			
			// get all the nodes we don't want
			iArrayT fNotNodes;
			model.ManyNodeSets(notTheseSets, fNotNodes);
			// get all the nodes
			fNodes.Dimension(model.NumNodes());
			fNodes.SetValueToPosition();
			// Take the complement
			for (int i = 0; i < fNotNodes.Length(); i++)
				fNodes[fNotNodes[i]] = -fNodes[fNotNodes[i]] - 1; // if the node is to be deleted, make it < 0
			fNodes.SortDescending();
			fNodes.Resize(fNodes.Length() - fNotNodes.Length());
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

/* describe the parameters needed by the interface */
void ThermostatBaseT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT beta(fBeta, "beta");
	beta.AddLimit(0.0, LimitT::LowerInclusive);
	list.AddParameter(beta);
}

/* accept parameter list */
void ThermostatBaseT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	fBeta = list.GetParameter("beta");
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
	const double* x_i;
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
