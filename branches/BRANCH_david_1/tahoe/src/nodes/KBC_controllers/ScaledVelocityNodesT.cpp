/* $Id: ScaledVelocityNodesT.cpp,v 1.8 2004-07-15 08:31:21 paklein Exp $ */
#include "ScaledVelocityNodesT.h"
#include "NodeManagerT.h"
#include "FEManagerT.h"
#include "ModelManagerT.h"

#include "RandomNumberT.h"
#include "MessageT.h"
#include "CommunicatorT.h"
#include "iArrayT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"

const double fkB = 0.00008617385;

using namespace Tahoe;

/* constructor */
ScaledVelocityNodesT::ScaledVelocityNodesT(const BasicSupportT& support, BasicFieldT& field):
	KBC_ControllerT(support),
	fField(field),
	fDummySchedule(1.0),
	qIConly(false),
	qFirstTime(false),
	qAllNodes(false),
	fIncs(0),
	fIncCt(0),
	fTempSchedule(NULL),
	fRandom(RandomNumberT::kParadynGaussian)
{
	SetName("scaled_velocity");
}

void ScaledVelocityNodesT::InitStep(void)
{
	/* really bad, this */
	if (qIConly)
	{
		if (!qFirstTime)
			qFirstTime = true;
		else
			qFirstTime = false;
	}
	else
	{
		fIncCt++;
		if (fIncCt == fIncs) // time to rescale velocities
		{
			fIncCt = 0;
			SetBCCards();
		}
	}
					
	/* inherited */
	KBC_ControllerT::InitStep();
	
	if (qIConly && !qFirstTime)
		fKBC_Cards.Dimension(0);
}

void ScaledVelocityNodesT::InitialCondition(void)
{
	if (!qIConly)
	{
		fT_0 = fTempScale*fTempSchedule->Value();
	}
	
	/* number of scaled nodes */
	int n_scaled = 0;
	int ndof = fField.NumDOF();
	
	/* workspace to generate velocities */
	dArray2DT velocities(fSupport.NumNodes(), ndof); 

	/* generate gaussian dist of random vels */
	fRandom.RandomArray(velocities);
	
	/* get MPI stuff */
	const CommunicatorT& communicator = fSupport.Communicator();
	int nProcs = fSupport.Size();
	int thisProc = fSupport.Rank();
	const ArrayT<int>* pMap = fSupport.ProcessorMap();
	
	double tKE = 0.;
	dArrayT vCOM(ndof);
	vCOM = 0.;
	//double totalMass = fMass*n_scaled;
	/* only change velocities for nodes on this processor */
	iArrayT myNodes; 
	
	if (fSupport.Size() == 1 || !pMap)
		myNodes.Set(fNodes.Length(), fNodes.Pointer());
	else
	{
		for (int i = 0; i < myNodes.Length(); i++)
			if ((*pMap)[fNodes[i]] == thisProc)
				myNodes[i] = fNodes[i];
			else
				myNodes[i] = -1;
	}
	
	for (int i = 0; i < myNodes.Length(); i++)
	{	
		if (myNodes[i] >= 0)
		{
			double* v_i = velocities(myNodes[i]);
		
			for (int j = 0; j < ndof; j++)
			{	
				tKE += (*v_i)*(*v_i);
				vCOM[j] += *v_i++;
			} 
			
			n_scaled++;
		}
	}
	int n_total = communicator.Sum(n_scaled);
	double KE_total = communicator.Sum(tKE);
	dArray2DT vCOM_all(nProcs, ndof);
	communicator.AllGather(vCOM, vCOM_all);
	
	for (int j = 0; j < ndof; j++)
		vCOM[j] = vCOM_all.ColumnSum(j);
	
	vCOM /= n_total;

	for (int i = 0; i < myNodes.Length(); i++)
	{
		if (myNodes[i] >= 0)
			for (int j = 0; j < ndof; j++)
				velocities(myNodes[i],j) -= vCOM[j];
	}
	
	/* adjust KE to COM frame  and convert it to a temperature */
	for (int j = 0; j < ndof; j++)
		KE_total -= n_total*vCOM[j]*vCOM[j];
	KE_total *= fMass/n_total/ndof/fkB;
	
	/* set the temperature */
	velocities *= sqrt(fT_0/KE_total);

	/* generate BC cards */
	fKBC_Cards.Dimension(n_scaled*ndof);
	KBC_CardT* pcard = fKBC_Cards.Pointer();
	for (int i = 0; i < myNodes.Length(); i++)
	{
		if (myNodes[i] >= 0)
		{
			double* v_i = velocities(myNodes[i]);	
			
	    	for (int j = 0; j < ndof; j++)
			{	
				/* set values */
				pcard->SetValues(myNodes[i], j, KBC_CardT::kVel, &fDummySchedule, *v_i++);
				pcard++;
			}
		} 
	}

}

/* describe the parameters needed by the interface */
void ScaledVelocityNodesT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	KBC_ControllerT::DefineParameters(list);

	/* mass parameter */
	ParameterT mass(ParameterT::Double, "mass");
	mass.AddLimit(0.0, LimitT::Lower);
	list.AddParameter(mass);

	/* random number seed */
	list.AddParameter(ParameterT::Integer, "seed");
}
	
/* information about subordinate parameter lists */
void ScaledVelocityNodesT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	KBC_ControllerT::DefineSubs(sub_list);

	/* when to apply scaling */
	sub_list.AddSub("BC_or_IC", ParameterListT::Once, true);

	/* method for defining affected particles */
	sub_list.AddSub("node_pick_choice", ParameterListT::Once, true);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* ScaledVelocityNodesT::NewSub(const StringT& name) const
{
	if (name == "BC_or_IC")
	{
		ParameterContainerT* bc_or_ic = new ParameterContainerT(name);
		bc_or_ic->SetListOrder(ParameterListT::Choice);
		
		ParameterContainerT bc("scale_as_BC");
		bc.AddParameter(ParameterT::Double, "scale");
		bc.AddParameter(ParameterT::Integer, "schedule");
		bc.AddParameter(ParameterT::Integer, "increments");
		bc_or_ic->AddSub(bc);

		ParameterContainerT ic("scale_as_IC");
		ic.AddParameter(ParameterT::Double, "temperature");
		bc_or_ic->AddSub(ic);
				
		return bc_or_ic;
	}
	else if (name == "node_pick_choice") 
	{
		ParameterContainerT* pick = new ParameterContainerT(name);
		pick->SetSubSource(this);
		pick->SetListOrder(ParameterListT::Choice);

		ParameterContainerT pick_all("pick_all_nodes");
		pick->AddSub(pick_all);
		pick->AddSub("pick_nodes_by_list");
//		pick->AddSub("pick_nodes_by_region");
		return pick;
	}
	else if (name == "pick_nodes_by_list")
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
	else if (name == "pick_nodes_by_region")
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
		return KBC_ControllerT::NewSub(name);
}

/* accept parameter list */
void ScaledVelocityNodesT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "ScaledVelocityNodesT::TakeParameterList";

	/* inherited */
	KBC_ControllerT::TakeParameterList(list);

	/* simple parameters */
	fMass = list.GetParameter("mass");
	int rseed = list.GetParameter("seed");
	fRandom.sRand(rseed);

	/* resolve IC/BC */
	const ParameterListT& bc_or_ic = list.GetListChoice(*this, "BC_or_IC");
	if (bc_or_ic.Name() == "scale_as_BC") {
		fTempScale = bc_or_ic.GetParameter("scale");
		int schedule = bc_or_ic.GetParameter("schedule");
		schedule--;
		fTempSchedule = fSupport.Schedule(schedule);
		if (!fTempSchedule) ExceptionT::BadInputValue(caller);
		qIConly = false;
	}
	else if (bc_or_ic.Name() == "scale_as_IC") {
		fT_0 = bc_or_ic.GetParameter("temperature");
		qIConly = true;
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized \"BC_or_IC\" \"%s\"",
			bc_or_ic.Name().Pointer());

	/* method for specifying affected nodes */
	const char choice[] = "node_pick_choice";
	const ParameterListT& pick = list.GetListChoice(*this, choice);
	qAllNodes = false;
	if (pick.Name() == "pick_all_nodes") {		
		qAllNodes = true;
		fNodes.Dimension(fSupport.NumNodes());
		fNodes.SetValueToPosition();
	}
	else if (pick.Name() == "pick_nodes_by_list")
		InitNodeSets(pick);
//	else if (pick.Name() == "pick_nodes_by_region")
//		InitRegion(*pick);
	else
		ExceptionT::GeneralFail(caller, "unrecognized pick method \"%s\"",
			pick.Name().Pointer());
}

/**********************************************************************
 * Protected
 **********************************************************************/

void ScaledVelocityNodesT::InitNodeSets(const ParameterListT& pick_nodes)
{
	const char caller[] = "ScaledVelocityNodesT::InitNodeSets";

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

/* initialize the current step */
void ScaledVelocityNodesT::SetBCCards(void)
{
	/* number of scaled nodes */
	int n_scaled = 0; 
	int ndof = fField.NumDOF();

	/* grab the velocities */
	const dArray2DT* velocities = NULL;
	if (fField.Order() > 0)
		velocities = &fField[1];		
	if (!velocities)
		ExceptionT::GeneralFail("ScaledVelocityNodesT::SetBCCards","Cannot get velocity field ");
		
	/* get MPI stuff */
	const CommunicatorT& communicator = fSupport.Communicator();
	int nProcs = fSupport.Size();
	int thisProc = fSupport.Rank();
	const ArrayT<int>* pMap = fSupport.ProcessorMap();	
	
	/* 	assume uniform mass for now */

	/* figure out which nodes to affect */
	iArrayT myNodes;

	if (fSupport.Size() == 1 || !pMap)
		myNodes.Set(fNodes.Length(), fNodes.Pointer());
	else
	{
		for (int i = 0; i < myNodes.Length(); i++)
			if ((*pMap)[fNodes[i]] == thisProc)
				myNodes[i] = fNodes[i];
			else
				myNodes[i] = -1;
	}


	/* calculate CM velocity and temperature */
	dArrayT vCOM(ndof);
	double tKE = 0.; // total kinetic energy
	double vscale; // scale velocity by this after subtracting off vCOM
	if (myNodes.Count(-1) != myNodes.Length())
	{
		vCOM = 0.;
		//double totalMass = fMass*n_scaled;
		for (int i = 0; i < myNodes.Length(); i++)
		{	
			if (myNodes[i] >= 0)
			{
				const double* v_i = (*velocities)(myNodes[i]);
			
				for (int j = 0; j < ndof; j++)
				{	
					tKE += (*v_i)*(*v_i);
					vCOM[j] += *v_i++;
				} 
				
				n_scaled++;
			}
		}
		int n_total = communicator.Sum(n_scaled);
		double KE_total = communicator.Sum(tKE);
		dArray2DT vCOM_all(nProcs, ndof);
		communicator.AllGather(vCOM, vCOM_all);
		
		for (int j = 0; j < ndof; j++)
			vCOM[j] = vCOM_all.ColumnSum(j);
		
		vCOM /= n_total;
		
		/* adjust KE to COM frame  and convert it to a temperature */
		for (int j = 0; j < ndof; j++)
			KE_total -= n_total*vCOM[j]*vCOM[j];
		KE_total *= fMass;
			
		/* want new KE to be ndof/2*n_scaled * kT */
		vscale = sqrt(ndof*n_total*fkB*fTempScale*(fTempSchedule->Value())/KE_total);
		
		/* generate BC cards */
		fKBC_Cards.Dimension(n_scaled*ndof);
		KBC_CardT* pcard = fKBC_Cards.Pointer();
		for (int i = 0; i < myNodes.Length(); i++)
		{
			if (myNodes[i] >= 0)
			{
				const double* v_i = (*velocities)(myNodes[i]);	
				
		    	for (int j = 0; j < ndof; j++)
				{	
					/* set values */
					pcard->SetValues(myNodes[i], j, KBC_CardT::kVel, &fDummySchedule, (*v_i++-vCOM[j])*vscale);
					pcard++;
				}
			} 
		}
	}
}

