/* $Id: ThermostatBaseT.cpp,v 1.12.20.1 2005-06-05 06:23:55 paklein Exp $ */
#include "ThermostatBaseT.h"

#include "BasicSupportT.h"

#include "dArrayT.h"
#include "dArray2DT.h"
#include "RaggedArray2DT.h"
#include "ParticlePropertyT.h"
#include "ModelManagerT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* copy behavior for arrays ThermostatBaseT*'s */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ThermostatBaseT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
ThermostatBaseT::ThermostatBaseT(const BasicSupportT& support):
	ParameterInterfaceT("velocity_damping"),
	fSupport(support),
	fBeta(0.0),
	fTemperature(0.0),
	fAllNodes(false),
	fTemperatureSchedule(NULL),
	fTemperatureScale(0.0)
{

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
			dArray2DT& forces, ArrayT<int>& types,
			ArrayT<ParticlePropertyT*>& particleProperties)
{
	int nsd = fSupport.NumSD();
	const double* v_j;
	double* f_j;
	int tag_j, currType;
	double mass, beta;
	if (fAllNodes)
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
				
			for (int i = 0; i < nsd; i++)
				*f_j++ -= beta*(*v_j++);
		}
	}
	else if (fNodes.Length() > 0)
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

			for (int i = 0; i < nsd; i++)
				*f_j++ -= beta*(*v_j++); 	
	    }
	}
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

/* information about subordinate parameter lists */
void ThermostatBaseT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* prescribed temperature */
	sub_list.AddSub("thermostat_temperature");

	/* method for defining affected particles */
	sub_list.AddSub("particle_pick_choice", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ThermostatBaseT::NewSub(const StringT& name) const
{
	if (name == "thermostat_temperature")
	{
		ParameterContainerT* temp = new ParameterContainerT(name);	
		temp->AddParameter(ParameterT::Integer, "schedule");
		temp->AddParameter(ParameterT::Double, "value");	
		return temp;
	}
	else if (name == "particle_pick_choice") 
	{
		ParameterContainerT* pick = new ParameterContainerT(name);
		pick->SetSubSource(this);
		pick->SetListOrder(ParameterListT::Choice);

		ParameterContainerT pick_all("pick_all");
		pick->AddSub(pick_all);
		pick->AddSub("pick_by_list");
		pick->AddSub("pick_by_region");
		return pick;
	}
	else if (name == "pick_by_list")
	{
		ParameterContainerT* pick = new ParameterContainerT(name);
		
		/* these nodes or all but these */
		ParameterT use_or_no(ParameterT::Enumeration, "selection");
		use_or_no.AddEnumeration("include_these", 0);
		use_or_no.AddEnumeration("exclude_these", 1);
		pick->AddParameter(use_or_no);	

		/* the node list */
		pick->AddSub("node_ID_list");
		return pick;
	}
	else if (name == "pick_by_region")
	{
		ParameterContainerT* pick = new ParameterContainerT(name);
		pick->SetDescription("provide bounds for each direction");

		ParameterT search_inc(ParameterT::Integer, "search_increment");
		search_inc.AddLimit(0, LimitT::LowerInclusive);
		pick->AddParameter(search_inc);
				
		ParameterContainerT bounds("pick_bounds");
		ParameterT direction(ParameterT::Integer, "direction");
		direction.AddLimit(1, LimitT::LowerInclusive);
		bounds.AddParameter(direction);
		bounds.AddParameter(ParameterT::Double, "x_min");
		bounds.AddParameter(ParameterT::Double, "x_max");
		pick->AddSub(bounds, ParameterListT::OnePlus);
		
		return pick;
	}
	else /* inherited */
		return ParameterInterfaceT::NewSub(name);
}

/* accept parameter list */
void ThermostatBaseT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "ThermostatBaseT::TakeParameterList";

	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* temperature schedule */
	const ParameterListT& temp_schedule = list.GetList("thermostat_temperature");
	fTemperatureScale = temp_schedule.GetParameter("value");
	int schedule = temp_schedule.GetParameter("schedule");
	fTemperatureSchedule = fSupport.Schedule(--schedule);

	/* damping parameter */
	fBeta = list.GetParameter("beta");
	
	/* method for specifying affected nodes */
	const char choice[] = "particle_pick_choice";
	const ParameterListT& pick = list.GetListChoice(*this, choice);
	fAllNodes = false;
	if (pick.Name() == "pick_all")
		fAllNodes = true;
	else if (pick.Name() == "pick_by_list")
		InitNodeSets(pick);
	else if (pick.Name() == "pick_by_region")
		InitRegion(pick);
	else
		ExceptionT::GeneralFail(caller, "unrecognized pick method \"%s\"", pick.Name().Pointer());
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void ThermostatBaseT::InitNodeSets(const ParameterListT& pick_nodes)
{
	const char caller[] = "ThermostatBaseT::InitNodeSets";

	/* model information */
	ModelManagerT& model = fSupport.ModelManager();

	/* determine selection method */
	const StringT& selection = pick_nodes.GetParameter("selection");
	if (selection == "include_these")
	{
		/* read node set ids */
		ArrayT<StringT> ids;
		StringListT::Extract(pick_nodes.GetList("node_ID_list"), ids);
		model.ManyNodeSets(ids, fNodes);        
	}
	else if (selection == "exclude_these")
	{
		/* read sets of nodes to omit */
		ArrayT<StringT> ids;
		StringListT::Extract(pick_nodes.GetList("node_ID_list"), ids);
		iArrayT not_nodes;
		model.ManyNodeSets(ids, not_nodes);

		/* get all the nodes */
		fNodes.Dimension(model.NumNodes());
		fNodes.SetValueToPosition();

		/* take the complement */
		for (int i = 0; i < not_nodes.Length(); i++)
			fNodes[not_nodes[i]] = -fNodes[not_nodes[i]] - 1; /* if the node is to be deleted, make it < 0 */
		fNodes.SortDescending();
		fNodes.Resize(fNodes.Length() - not_nodes.Length());
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized selection method \"%s\"",
			selection.Pointer());
}	

void ThermostatBaseT::InitRegion(const ParameterListT& pick_region)
{
	const char caller[] = "ThermostatBaseT::InitRegion";

	/* search increment */
	nIncs = pick_region.GetParameter("search_increment");
	
	/* bounds */
	int nsd = fSupport.NumSD();
	const char bounds[] = "pick_bounds";
	if (pick_region.NumLists(bounds) != nsd)
		ExceptionT::GeneralFail(caller, "expecting %d \"%s\" not %d in \"%s\"",
			nsd, bounds, pick_region.NumLists(bounds), pick_region.Name().Pointer());

	fxmin.Dimension(nsd);
	fxmax.Dimension(nsd);
	for (int i = 0; i < nsd; i++) {
		const ParameterListT& bound = pick_region.GetList(bounds,i);
		int direction = bound.GetParameter("direction");
		direction--;
		fxmin[direction] = bound.GetParameter("x_min");
		fxmax[direction] = bound.GetParameter("x_max");

		/* check */
		if (fxmin[direction] >= fxmax[direction])
			ExceptionT::BadInputValue(caller,"bad bounding box coordinates");
	}
	
	/* get the nodes in the region */
	NodesInRegion(fSupport.CurrentCoordinates(), fSupport.PartitionNodes());
}	

void ThermostatBaseT::NodesInRegion(const dArray2DT& coords, const ArrayT<int>* partition_nodes)
{
	if (fxmin.Length() != coords.MinorDim())
		ExceptionT::GeneralFail("ThermostattedRegionT::NodesInRegion",
				"Dimension mismatch between coords and bounding box");
	AutoArrayT<int> tmpList;

	int nsd = fSupport.NumSD();
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
		for (int j = 0; inBox && j < nsd; j++)
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
