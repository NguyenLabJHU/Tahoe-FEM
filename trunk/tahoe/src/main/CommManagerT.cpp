/* $Id: CommManagerT.cpp,v 1.1 2002-12-05 08:31:13 paklein Exp $ */
#include "CommManagerT.h"
#include "CommunicatorT.h"
#include "ModelManagerT.h"
#include "PartitionT.h"
#include <float.h>

using namespace Tahoe;

CommManagerT::CommManagerT(CommunicatorT& comm, int nsd):
	fComm(comm),
	fIsPeriodic(nsd),
	fPeriodicBoundaries(nsd, 2),
	fBounds(nsd, 2),
	fPartition(NULL),
	fModelManager(NULL)
{
	fIsPeriodic = false;
	fPeriodicBoundaries = 0.0;
}

/* set boundaries */
void CommManagerT::SetPeriodicBoundaries(int i, double x_i_min, double x_i_max)
{
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
	if (!fModelManager) ExceptionT::GeneralFail("CommManagerT::Configure");

	/* nothing to do */
	if (fComm.Size() == 1)
	{
		bool has_periodic = false;
		for (int i = 0; i < fIsPeriodic.Length(); i++)
			has_periodic = has_periodic || fIsPeriodic[i];
		if (!has_periodic) return;
	}
	if (!fPartition || fPartition->DecompType() != PartitionT::kAtom) return;

	/* number of nodes owned by this partition */
	int num_part_nodes = fPartition->NumPartitionNodes();

	/* reference coordinates */
	const dArray2DT& reference_coords = fModelManager->Coordinates();

	/* local node manager does not have the total number of nodes */
	if (num_part_nodes != reference_coords.MajorDim())
	{
	
	
	}

	/* determine the current processor bounds */
	//GetBounds(current_coords, , fBounds);
}

/*************************************************************************
 * Protected
 *************************************************************************/

/** determine the local coordinate bounds */
void CommManagerT::GetBounds(const dArray2DT& coords_all, const iArrayT& local, 
	dArray2DT& bounds) const
{
	/* dimensions */
	int nsd = bounds.MajorDim();

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
