/* $Id: CommManagerT.cpp,v 1.1.2.1 2002-12-05 21:48:17 paklein Exp $ */
#include "CommManagerT.h"
#include "CommunicatorT.h"
#include "ModelManagerT.h"
#include "PartitionT.h"
#include "AllGatherT.h"
#include <float.h>

using namespace Tahoe;

CommManagerT::CommManagerT(CommunicatorT& comm, ModelManagerT& model_manager):
	fComm(comm),
	fModelManager(model_manager),
	fPartition(NULL),
	fPeriodicBoundaries(0,2),
	fFirstConfigure(true)
{

}

/* set boundaries */
void CommManagerT::SetPeriodicBoundaries(int i, double x_i_min, double x_i_max)
{
	/* dimension dynamically */
	if (i >= fIsPeriodic.Length())
	{
		int len = i + 1;
		fIsPeriodic.Resize(len, false);
		fPeriodicBoundaries.Resize(len, 0.0);
	}

	if (x_i_min > x_i_max) ExceptionT::GeneralFail();
	fIsPeriodic[i] = true;
	fPeriodicBoundaries(i,0) = x_i_min;
	fPeriodicBoundaries(i,1) = x_i_max;

//TEMP
ExceptionT::Stop("CommManagerT::SetPeriodicBoundaries", "not supported yet");
}

/* unset boundaries */
void CommManagerT::ClearPeriodicBoundaries(int i) { fIsPeriodic[i] = false; }

/* configure local environment */
void CommManagerT::Configure(double range)
{
	/* nothing to do */
	if (fComm.Size() == 1)
	{
		bool has_periodic = false;
		for (int i = 0; i < fIsPeriodic.Length(); i++)
			has_periodic = has_periodic || fIsPeriodic[i];
		if (!has_periodic) return;
	}
	if (!fPartition || fPartition->DecompType() != PartitionT::kAtom) return;

	/* first time through */
	if (fFirstConfigure) {
		FirstConfigure();
		fFirstConfigure = false;
	}

	/* determine the current processor bounds */
	//GetBounds(current_coords, , fBounds);
}

/* perform actions needed the first time CommManagerT::Configure is called. */
void CommManagerT::FirstConfigure(void)
{
	/* nothing to do in serial */
	if (!fPartition) return;
	
	if (fPartition->DecompType() == PartitionT::kAtom)
	{
		/* number of nodes owned by this partition */
		int npn = fPartition->NumPartitionNodes();

		/* total number of nodes */
		int nsd = fModelManager.NumDimensions();
		AllGatherT all_gather(fComm);
		all_gather.Initialize(npn*nsd);
		int ntn = all_gather.Total()/nsd;

		/* collect global reference coordinate list */
		const dArray2DT& reference_coords = fModelManager.Coordinates();
		if (reference_coords.MajorDim() != ntn)
		{
			/* collect */
			dArray2DT coords_all(ntn, nsd);		
			all_gather.AllGather(reference_coords, coords_all);

			/* reset model manager */
			fModelManager.UpdateNodes(coords_all, true);
		}

		/* translate element sets to global numbering */
		const ArrayT<StringT>& el_id = fModelManager.ElementGroupIDs();
		for (int i = 0; i < el_id.Length(); i++)
		{
			/* copy existing connectivities */
			iArray2DT elems = fModelManager.ElementGroup(el_id[i]);

			/* map scope of nodes in connectivities */
			fPartition->SetNodeScope(PartitionT::kGlobal, elems);
			
			/* re-set the connectivities */
			fModelManager.UpdateElementGroup(el_id[i], elems, true);
		}

		/* translate node sets to global numbering */
		const ArrayT<StringT>& nd_id = fModelManager.NodeSetIDs();
		for (int i = 0; i < nd_id.Length(); i++)
		{
			/* copy existing connectivities */
			iArrayT nodes = fModelManager.NodeSet(nd_id[i]);

			/* map scope of nodes in connectivities */
			fPartition->SetNodeScope(PartitionT::kGlobal, nodes);
			
			/* re-set the connectivities */
			fModelManager.UpdateNodeSet(nd_id[i], nodes, true);		
		}

		/* translate side sets to global numbering */	
		const ArrayT<StringT>& ss_id = fModelManager.SideSetIDs();
		//TEMP
		if (ss_id.Length() > 0) cout << "\n CommManagerT::FirstConfigure: skipping size sets" << endl;
	}
}

/*************************************************************************
 * Protected
 *************************************************************************/

/** determine the local coordinate bounds */
void CommManagerT::GetBounds(const dArray2DT& coords_all, const iArrayT& local, 
	dArray2DT& bounds) const
{
	/* dimensions */
	int nsd = coords_all.MinorDim();
	bounds.Dimension(nsd, 2);

	/* initialize */
	if (local.Length() == 0)
	{
		bounds = 0.0;
		return;
	}
	else
	{
		bounds.SetColumn(0, DBL_MAX);
		bounds.SetColumn(1, DBL_MIN);
	}

	/* resolve by spatial dimension */
	if (nsd == 1)
	{
		int* p_i = local.Pointer();
		double& x_min = bounds(0,0);
		double& x_max = bounds(0,1);
		for (int i = 0; i < local.Length(); i++)
		{
			double* x = coords_all(*p_i++);
			x_min = (*x < x_min) ? *x : x_min;
			x_max = (*x > x_max) ? *x : x_max;
		}
	}
	else if (nsd == 2)
	{
		int* p_i = local.Pointer();
		double& x_min = bounds(0,0);
		double& x_max = bounds(0,1);
		double& y_min = bounds(1,0);
		double& y_max = bounds(1,1);
		for (int i = 0; i < local.Length(); i++)
		{
			double* x = coords_all(*p_i++);
			x_min = (*x < x_min) ? *x : x_min;
			x_max = (*x > x_max) ? *x : x_max;
			x++;
			y_min = (*x < y_min) ? *x : y_min;
			y_max = (*x > y_max) ? *x : y_max;
		}
	}
	else if (nsd == 3)
	{
		int* p_i = local.Pointer();
		double& x_min = bounds(0,0);
		double& x_max = bounds(0,1);
		double& y_min = bounds(1,0);
		double& y_max = bounds(1,1);
		double& z_min = bounds(2,0);
		double& z_max = bounds(2,1);
		for (int i = 0; i < local.Length(); i++)
		{
			double* x = coords_all(*p_i++);
			x_min = (*x < x_min) ? *x : x_min;
			x_max = (*x > x_max) ? *x : x_max;
			x++;
			y_min = (*x < y_min) ? *x : y_min;
			y_max = (*x > y_max) ? *x : y_max;
			x++;
			z_min = (*x < z_min) ? *x : z_min;
			z_max = (*x > z_max) ? *x : z_max;
		}
	}
	else ExceptionT::OutOfRange("CommManagerT::SetBounds");
}
