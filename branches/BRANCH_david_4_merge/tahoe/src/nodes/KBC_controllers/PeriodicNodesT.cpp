/* $Id: PeriodicNodesT.cpp,v 1.6 2004-09-01 10:46:27 paklein Exp $ */
#include "PeriodicNodesT.h"
#include "BasicSupportT.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

/* constructor */
PeriodicNodesT::PeriodicNodesT(const BasicSupportT& support, BasicFieldT& field):
	TiedNodesT(support, field)
{
	SetName("periodic_nodes");
}

/* information about subordinate parameter lists */
void PeriodicNodesT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	TiedNodesT::DefineSubs(sub_list);

	/* periodic strides */
	sub_list.AddSub("periodic_stride", ParameterListT::OnePlus);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* PeriodicNodesT::NewSub(const StringT& name) const
{
	if (name == "periodic_stride")
	{
		ParameterContainerT* periodic_stride = new ParameterContainerT(name);
	
		periodic_stride->AddParameter(ParameterT::Integer, "direction");

		ParameterT stride(ParameterT::Double, "stride");
		stride.AddLimit(0.0, LimitT::Lower);
		periodic_stride->AddParameter(stride);
	
		return periodic_stride;
	}
	else /* inherited */
		return TiedNodesT::NewSub(name);
}

/* accept parameter list */
void PeriodicNodesT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "PeriodicNodesT::TakeParameterList";

	/* extract periodic strides before calling inherited method */
	int nsd = fSupport.NumSD();
	fIsPeriodic.Dimension(nsd);
	fPeriodicStride.Dimension(nsd);
	fIsPeriodic = false;
	fPeriodicStride = 0.0;
	int num_stride = list.NumLists("periodic_stride");
	if (num_stride > nsd) 
		ExceptionT::GeneralFail(caller, "expecting at most %d not %d \"periodic_stride\"",
			nsd, num_stride);
	for (int i = 0; i < num_stride; i++)
	{
		const ParameterListT& periodic_stride = list.GetList("periodic_stride", i);
		
		/* extract direction */
		int direction = periodic_stride.GetParameter("direction"); direction--;
		if (direction < 0 || direction >= nsd)
			ExceptionT::GeneralFail(caller, "direction %d out of range {1,%d}",
				direction+1, nsd);
		
		/* store */
		fIsPeriodic[direction] = true;
		fPeriodicStride[direction] = periodic_stride.GetParameter("stride");
	}
	
	/* inherited */
	TiedNodesT::TakeParameterList(list);
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* set initial tied node pairs */
void PeriodicNodesT::InitTiedNodePairs(const iArrayT& leader_nodes, 
	iArrayT& follower_nodes)
{
	/* coordinates */
	const dArray2DT& coords = fSupport.InitialCoordinates();
	
	/* get processor number */
	int np = fSupport.Rank();
	const ArrayT<int>* pMap = fSupport.ProcessorMap();

	/* dumb search */
	int nsd = coords.MinorDim();

	/* Length may change during search if external nodes are removed */
	int FLength = follower_nodes.Length()-1;
	int LLength = leader_nodes.Length()-1;

	/* num of followers and leaders to be removed */	
	int Fct = 0, Lct = 0;
	
	/* periodic strides */
	dArrayT dx(nsd);
	for (int i = 0; i < dx.Length(); i++)
//		dx[i] = (fIsPeriodic[i]) ? fPeriodicStride[i] : 0.0; // for some reason this doesn't
		if (fIsPeriodic[i])                                  // compile correctly on GNU-Darwin
			dx[i] = fPeriodicStride[i];
		else
			dx[i] = 0.0;

	for (int i = 0; i <= FLength; i++)
	{
		const double* x_f = coords(follower_nodes[i]);
		
		/*If a follower is external, flag it for removal from the list*/
		if (pMap && (*pMap)[follower_nodes[i]] != np)
		{
			fPairStatus[i] = kChangeF;
		}
		
		for (int j = 0; j < leader_nodes.Length(); j++)
		{
			const double* x_l = coords(leader_nodes[j]);
			bool OK = true;
			for (int k = 0; OK && k < nsd; k++)
			{
				/* try all periodic images */
				OK = (fabs(x_f[k] - x_l[k]) < kSmall) ||
					 (fabs(x_f[k] - x_l[k] + dx[k]) < kSmall) ||
					 (fabs(x_f[k] - x_l[k] - dx[k]) < kSmall);
			}
				
			/* found */
			if (OK) 
			{
				fNodePairs(i,1) = leader_nodes[j];
				if (pMap && (*pMap)[follower_nodes[i]] != np && 
					(*pMap)[leader_nodes[j]] != np)
				{
					/* Flag the pair as external */
					fPairStatus[i] = kTiedExt;
				}
				else
					fPairStatus[i] = kTied;
			}
		}
			
		if (fPairStatus[i] > kTied)
		{
			/* If something will get removed, send it to the
			 * end of the list and resize after the search */
			follower_nodes[i] = follower_nodes[FLength];
			Fct++;
			fPairStatus[i] = kFree;
			fNodePairs(i,0) = fNodePairs(FLength,0);
			fNodePairs(i,1) = fNodePairs(FLength,1);
			if (i != follower_nodes.Length()-1) 
				i--;
			FLength--;
		}
	}

	if (Fct > 0)
	{
		fPairStatus.Resize(fPairStatus.Length()-Fct);
		follower_nodes.Resize(follower_nodes.Length()-Fct);
		fNodePairs.Resize(fNodePairs.MajorDim()-Fct);
	}
}
