/* $Id: SpatialGridT.cpp,v 1.1 2004-10-06 21:07:10 paklein Exp $ */
#include "SpatialGridT.h"
#include "dArray2DT.h"
#include "iArrayT.h"
#include "CommunicatorT.h"

#include <float.h>

using namespace Tahoe;

/* constructor */
SpatialGridT::SpatialGridT(GridBoundT grid_bound):
	fGridBound(grid_bound)
{

}

/* set number of grid cells in each dimension */
void SpatialGridT::Dimension(const iArrayT& num_cells)
{
	const char caller[] = "SpatialGridT::Dimension";

	/* grid dimensions */
	fNx = num_cells;
	if (fNx.Min() < 1) ExceptionT::GeneralFail(caller, "bad dimension %d", fNx.Min());
	
	/* set numbering indexing offsets */
	fN_offset.Dimension(fNx.Length() - 1);
	int offset = 1;
	for (int i = fN_offset.Length() - 1; i > -1; i--) {
		offset *= fNx[i+1];
		fN_offset[i] = offset;
	}
	
	/* dimension other work space */
	fMinMax.Dimension(fNx.Length(), 2);
	fMinMax = 0.0;
	fdx.Dimension(fNx.Length());
	fdx = 0.0;
}

/* set bounds */
void SpatialGridT::SetBounds(const dArray2DT& min_max)
{
	const char caller[] = "SpatialGridT::SetBounds";

	/* check dimensions */
	if (min_max.MajorDim() != fMinMax.MajorDim() || min_max.MinorDim() != 2)
	    ExceptionT::SizeMismatch(caller);

	/* copy */
	fMinMax = min_max;
	
	/* grid sizes */
	for (int i = 0; i < fMinMax.MajorDim(); i++) {	
		double min = fMinMax(i,0);
		double max = fMinMax(i,1);
		if (max < min) ExceptionT::GeneralFail(caller, "max < min: %g < %g", max, min);
	
		fdx[i] = (max - min)/fNx[i];	
	}
}

/* set bounds based on points */
void SpatialGridT::SetBounds(CommunicatorT& comm, const dArray2DT& points, const ArrayT<int>* points_used)
{
	const char caller[] = "SpatialGridT::SetBounds";
	if (points.MinorDim() != fMinMax.MajorDim()) ExceptionT::SizeMismatch(caller);

	/* initialize */
	dArray2DT min_max(points.MinorDim(), 2);
	min_max.SetColumn(0, DBL_MAX);
	min_max.SetColumn(1,-DBL_MAX);

	/* search */
	int npt = (points_used) ? points_used->Length() : points.MajorDim();
	int nsd = points.MinorDim();
	for (int i = 0; i < nsd; i++) 
	{	
		/* limits */
		double& min = min_max(i,0);
		double& max = min_max(i,1);

		for (int j = 0; j < npt; j++)
		{
			int index = (points_used) ? (*points_used)[j] : j;		
			double x = points(index,i);
			if (x < min)
				min = x;
			else if (x > max)
				max = x;
		}	
	}

	/* globalize */
	for (int i = 0; i < min_max.MajorDim(); i++) {
		min_max(i,0) = comm.Min(min_max(i,0));
		min_max(i,1) = comm.Max(min_max(i,1));
	}

	/* set limits */
	SetBounds(min_max);
}

/* sort coordinates into bins */
void SpatialGridT::Bin(const dArray2DT& points, iArrayT& bin, iArrayT& bin_counts, const ArrayT<int>* points_used)
{
	const char caller[] = "SpatialGridT::Bin";
	if (points.MinorDim() != fMinMax.MajorDim()) ExceptionT::SizeMismatch(caller);

	/* dimensions */
	int npt = (points_used) ? points_used->Length() : points.MajorDim();
	int nsd = points.MinorDim();
	bin.Dimension(npt);
	bin_counts.Dimension(fNx.Product());
	
	/* initialize */
	bin = -1;
	bin_counts = 0;	

	if (nsd == 2)
		Bin2D(points, bin, bin_counts, points_used);
	else if (nsd == 3)
		Bin3D(points, bin, bin_counts, points_used);
	else
		ExceptionT::GeneralFail(caller, "unsupported dimension %d", nsd);
}

void SpatialGridT::Bin2D(const dArray2DT& points, iArrayT& bin, iArrayT& bin_counts,  
	const ArrayT<int>* points_used)
{
	const char caller[] = "SpatialGridT::Bin2D";

	/* grid dimensions */
	int nx = fNx[0];
	int ny = fNx[1];
	
	/* cell dimensions */
	double dx = fdx[0];
	double dy = fdx[1];
	
	/* grid lower bounds */
	double xmin = fMinMax(0,0);
	double ymin = fMinMax(1,0);

	/* sort into bins */
	int npt = bin.Length();
	for (int i = 0; i < npt; i++)
	{
		int index = (points_used) ? (*points_used)[i] : i;
		const double* px = points(index);

		/* grid index */
		int ix = int((px[0] - xmin)/dx);
		int iy = int((px[1] - ymin)/dy);
		
		/* beyond grid bounds */
		bool in_bounds = true;
		if (ix < 0 || ix >= nx) {
			if (fGridBound == kCutOff)
				in_bounds = false;
			else if (fGridBound == kExtended) {
				ix = (ix < 0) ? 0 : ix;
				ix = (ix >= nx) ? nx - 1: ix;
			}
			else /* fGridBound = kError */
				ExceptionT::GeneralFail(caller);
		}
		if (iy < 0 || iy >= ny) {
			if (fGridBound == kCutOff)
				in_bounds = false;
			else if (fGridBound == kExtended) {
				iy = (iy < 0) ? 0 : iy;
				iy = (iy >= ny) ? ny - 1: iy;
			}
			else /* fGridBound = kError */
				ExceptionT::GeneralFail(caller);
		}
		
		/* bin index */
		if (in_bounds) {
			int nbin = ix*ny + iy;
			bin[index] = nbin;
			bin_counts[nbin]++;
		}
	}	
}

void SpatialGridT::Bin3D(const dArray2DT& points, iArrayT& bin, iArrayT& bin_counts,  
	const ArrayT<int>* points_used)
{
	const char caller[] = "SpatialGridT::Bin3D";

	/* grid dimensions */
	int nx = fNx[0];
	int ny = fNx[1];
	int nz = fNx[2];
	
	int ix_jump = nz*ny;
	int iy_jump = nz;
	
	/* cell dimensions */
	double dx = fdx[0];
	double dy = fdx[1];
	double dz = fdx[2];
	
	/* grid lower bounds */
	double xmin = fMinMax(0,0);
	double ymin = fMinMax(1,0);
	double zmin = fMinMax(2,0);

	/* sort into bins */
	int npt = bin.Length();
	for (int i = 0; i < npt; i++)
	{
		int index = (points_used) ? (*points_used)[i] : i;		
		const double* px = points(index);

		/* grid index */
		int ix = int((px[0] - xmin)/dx);
		int iy = int((px[1] - ymin)/dy);
		int iz = int((px[2] - zmin)/dz);
		
		/* beyond grid bounds */
		bool in_bounds = true;
		if (ix < 0 || ix >= nx) {
			if (fGridBound == kCutOff)
				in_bounds = false;
			else if (fGridBound == kExtended) {
				ix = (ix < 0) ? 0 : ix;
				ix = (ix >= nx) ? nx - 1: ix;
			}
			else /* fGridBound = kError */
				ExceptionT::GeneralFail(caller);
		}
		if (iy < 0 || iy >= ny) {
			if (fGridBound == kCutOff)
				in_bounds = false;
			else if (fGridBound == kExtended) {
				iy = (iy < 0) ? 0 : iy;
				iy = (iy >= ny) ? ny - 1: iy;
			}
			else /* fGridBound = kError */
				ExceptionT::GeneralFail(caller);
		}
		if (iz < 0 || iz >= nz) {
			if (fGridBound == kCutOff)
				in_bounds = false;
			else if (fGridBound == kExtended) {
				iz = (iz < 0) ? 0 : iz;
				iz = (iz >= nz) ? nz - 1: iz;
			}
			else /* fGridBound = kError */
				ExceptionT::GeneralFail(caller);
		}

		/* bin index */
		if (in_bounds) {
			int nbin = ix*ix_jump + iy*iy_jump + iz;
			bin[index] = nbin;
			bin_counts[nbin]++;
		}
	}
}
