/* $Id: PeriodicNodesT.cpp,v 1.4 2003-11-21 22:47:59 paklein Exp $ */
#include "PeriodicNodesT.h"
#include "NodeManagerT.h"
#include "ifstreamT.h"

using namespace Tahoe;

/* constructor */
PeriodicNodesT::PeriodicNodesT(NodeManagerT& node_manager, BasicFieldT& field):
	TiedNodesT(node_manager, field),
	fIsPeriodic(fNodeManager.NumSD()),
	fPeriodicStride(fNodeManager.NumSD())
{
	fIsPeriodic = false;
	fPeriodicStride = 0.0;
}

void PeriodicNodesT::WriteParameters(ostream& out) const
{
	/* inherited */
	TiedNodesT::WriteParameters(out);

	/* periodic information */
	int d_width = OutputWidth(out, fPeriodicStride.Pointer());
	out << " Periodic boundaries:\n";
	out << setw(kIntWidth) << "dir" << setw(d_width) << "length" << '\n';
	for (int i = 0; i < fPeriodicStride.Length(); i++)
	{
		out << setw(kIntWidth) << i+1;
		if (fIsPeriodic[i]) 
			out << setw(d_width) << fPeriodicStride[i] << '\n';
		else
			out << setw(d_width) << "-" << '\n';
	}
	out.flush();
}

/**********************************************************************
 * Protected
 **********************************************************************/

void PeriodicNodesT::ReadParameters(ifstreamT& in)
{
	/* inherited */
	TiedNodesT::ReadParameters(in);
	
	/* periodic information */
	for (int i = 0; i < fIsPeriodic.Length(); i++)
	{
		int tf;
		double s;
		in >> tf >> s;
		if (tf == 1) {
			fPeriodicStride[i] = s;
			fIsPeriodic[i] = true;
			if (fPeriodicStride[i] < kSmall) throw ExceptionT::kBadInputValue;
		}
	}
}

/* set initial tied node pairs */
void PeriodicNodesT::InitTiedNodePairs(const iArrayT& leader_nodes, 
	iArrayT& follower_nodes)
{
	/* coordinates */
	const dArray2DT& coords = fNodeManager.InitialCoordinates();
	
	/* get processor number */
	int np = fNodeManager.Rank();
	const ArrayT<int>* pMap = fNodeManager.ProcessorMap();

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
