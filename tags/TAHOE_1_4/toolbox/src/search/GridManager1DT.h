/* $Id: GridManager1DT.h,v 1.11 2004-03-18 01:15:39 paklein Exp $ */
#ifndef _GRIDMANAGER1D_T_H_
#define _GRIDMANAGER1D_T_H_

#include "Environment.h"
#include "toolboxConstants.h"

/* language support */
#include <iostream.h>

/* direct members */
#include "AutoArrayT.h"
#include "pArrayT.h"
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class iArrayT;

template <class sTYPE>
class GridManager1DT
{
public:

	/* constructors */
	GridManager1DT(double xmin, double xmax, int nx);
	GridManager1DT(int nx, const dArray2DT& coords,
		const ArrayT<int>* nodes_used);
	
	/* destructor */
	~GridManager1DT(void);
	
	/* empty grid */
	void Reset(void);
	void Reset(const dArray2DT& coords, const ArrayT<int>* nodes_used);

	/* insert data into the grid */
	void Add(const sTYPE& data);

	/* closest point - error if grid is empty */
	const sTYPE& Closest(double* target);
	
	/* return list of data falling within the defined region */
	const AutoArrayT<sTYPE>& HitsInRegion(const double* coords, double distance);
	const AutoArrayT<sTYPE>& HitsInRegion(const double* coords, int cell_span);
	const AutoArrayT<sTYPE>& HitsInRegion(const double* coords, const ArrayT<double>& dist_xy);

	/* the distance covered by the given cell span */
	double CellSpan(int cell_span) const;

	/* grid statistics */
	void WriteStatistics(ostream& out) const;

	/** the grid */
	const ArrayT<AutoArrayT<sTYPE>*>& Grid(void) const { return fGrid; };

protected:

	/* return pointer to the content list for the given coords */
	AutoArrayT<sTYPE>** FetchGrid(const double* coords);

protected:

	/* grid dimensions */
	double fxmin;
	double fxmax;
	double fdx;
	int	   fnx;
	
	/* array of pointers to lists */
	pArrayT<AutoArrayT<sTYPE>*> fGrid;
	
	/* hit list */	
	AutoArrayT<sTYPE> fHits;

};

/*************************************************************************
 * Implementation
 *************************************************************************/

/* constructors */
template <class sTYPE>
GridManager1DT<sTYPE>::GridManager1DT(double xmin, double xmax, int nx):
	fxmin(xmin),
	fxmax(xmax),
	fnx(nx),
	fGrid(fnx)
{
	/* consistency */
	if (fxmax <= fxmin) throw ExceptionT::kGeneralFail;

	/* grid spacings */
	fdx = (fxmax - fxmin)/fnx;
}	

template <class sTYPE>
GridManager1DT<sTYPE>::GridManager1DT(int nx, const dArray2DT& coords,
	const ArrayT<int>* nodes_used):
	fnx(nx)
{
	/* initialize grid data */
	Reset(coords, nodes_used);
}	

template <class sTYPE>
GridManager1DT<sTYPE>::~GridManager1DT(void)
{

}	

/* empty grid */
template <class sTYPE>
void GridManager1DT<sTYPE>::Reset(void)
{
	AutoArrayT<sTYPE>** pgrid = fGrid.Pointer();

	for (int i = 0; i < fGrid.Length(); i++)
	{
		if (*pgrid) (*pgrid)->Dimension(0);
		pgrid++;
	}

}

template <class sTYPE>
void GridManager1DT<sTYPE>::Reset(const dArray2DT& coords,
	const ArrayT<int>* nodes_used)
{
	/* empty grid */
	Reset();

	int n_max;
	if (nodes_used == NULL) /* all coordinates */
	{
		n_max = coords.MajorDim();
		if (coords.MajorDim() > 0)
		{
			/* initialize limits */
			fxmin = fxmax = coords(0,0);
			
			for (int i = 0; i < coords.MajorDim(); i++)
			{
				const double* p = coords(i);
			
				if (p[0] < fxmin) fxmin = p[0];
				else if (p[0] > fxmax) fxmax = p[0];
			}
		}
		else
		{
			fxmin = fxmax = 0.0;
		}
	}
	else /* selected nodes */
	{
		n_max = nodes_used->Length();
		if (nodes_used->Length() > 0)
		{
			const int* dex = nodes_used->Pointer();

			/* initialize limits */
			fxmin = fxmax = coords(*dex, 0);
		
			for (int i = 0; i < nodes_used->Length(); i++)
			{
				const double* p = coords(*dex++);
		
				if (p[0] < fxmin) fxmin = p[0];
				else if (p[0] > fxmax) fxmax = p[0];
			}
		}
		else
		{
			fxmin = fxmax = 0.0;
		}		
	}

	/* min n_max */
	n_max = (n_max < 1) ? 1 : n_max;

	/* add a little space */
	fdx = 0.0001*(fxmax - fxmin);
	
	/* check for zero width */
	fdx = (fdx < kSmall) ? 1.0 : fdx;
	
	fxmax += fdx;
	fxmin -= fdx;

	/* grid spacings */
	fdx = (fxmax - fxmin)/fnx;

	/* set grid parameters */
	fGrid.Dimension(fnx);
}	

/* insert data into the grid */
template <class sTYPE>
void GridManager1DT<sTYPE>::Add(const sTYPE& data)
{
	/* get content list */
	AutoArrayT<sTYPE>** griddata = FetchGrid(data.Coords());

	/* initialize new list */
	if (!(*griddata))
	{
		*griddata = new AutoArrayT<sTYPE>;
		if (!*griddata) throw ExceptionT::kOutOfMemory;
	}

	/* append value */
	(*griddata)->Append(data);	
}

/* closest point */
template <class sTYPE>
const sTYPE& GridManager1DT<sTYPE>::Closest(double* target)
{
	/* broaden search */
	double distance = fdx/4.0; //so initial is 1/2 of spacing
	do {
	
		distance *= 2;
		HitsInRegion(target, distance);
	
	} while (fHits.Length() == 0 && (target[0] - distance > fxmin ||
	                                 target[0] + distance < fxmax));

	/* grid is empty ! */
	if (fHits.Length() == 0)
	{
		/* make blank reference */
		fHits[0].Clear();
		return fHits[0];
	}

	/* find closest to target */
	sTYPE*   pdata = fHits.Pointer();
	double* coords = pdata->Coords();

	double dx = coords[0] - target[0];

	double minsqr = dx*dx;
	int    mindex = 0;
	for (int i = 1; i < fHits.Length(); i++)
	{
		coords = (++pdata)->Coords();
		dx = coords[0] - target[0];
		double dsqr = dx*dx;
	
		/* found new min */
		if (dsqr < minsqr)
		{
			minsqr = dsqr;
			mindex = i;
		}	
	}
	
	return fHits[mindex];
}

/* the distance covered by the given cell span */
template <class sTYPE>
inline double GridManager1DT<sTYPE>::CellSpan(int cell_span) const
{
	/* max cell dimension */
	double dmax = fdx;
	return (cell_span*dmax)/2.0;
}

/* return list of data falling within the defined region */
template <class sTYPE>
inline const AutoArrayT<sTYPE>& GridManager1DT<sTYPE>::
	HitsInRegion(const double* coords, int cellspan)
{
	return HitsInRegion(coords, CellSpan(cellspan));
}	

/* return list of data falling within the defined region */
template <class sTYPE>
const AutoArrayT<sTYPE>& GridManager1DT<sTYPE>::
	HitsInRegion(const double* coords, double distance)
{
  /* NOT FINISHED CHANGING YET! */
	/* empty hit list */
	fHits.Dimension(0);

	/* grid indices */
	int ixstart = int((coords[0] - fxmin - distance)/fdx);
	int ixstop  = int((coords[0] - fxmin + distance)/fdx);

	/* keep within grid */
	ixstart = (ixstart < 0) ? 0 : ixstart;
	ixstop = (ixstop >= fnx) ? fnx - 1 : ixstop;
	bool out_of_range = ixstart > ixstop;
	if (!out_of_range)
	{
		/* scan section of grid */
		AutoArrayT<sTYPE>** griddata = fGrid.Pointer(ixstart);
		for (int ix = ixstart; ix <= ixstop; ix++)
		{
			if (*griddata) fHits.Append(**griddata);
			griddata++;
		}
	}	
	return fHits;
}	

template <class sTYPE>
const AutoArrayT<sTYPE>& GridManager1DT<sTYPE>::
	HitsInRegion(const double* coords, const ArrayT<double>& dist_x)
{
	/* empty hit list */
	fHits.Dimension(0);

	/* grid indices */
	int ixstart = int((coords[0] - fxmin - dist_x[0])/fdx);
	int ixstop  = int((coords[0] - fxmin + dist_x[0])/fdx);

	/* keep within grid */
	ixstart = (ixstart < 0) ? 0 : ixstart;
	ixstop = (ixstop >= fnx) ? fnx - 1 : ixstop;
	bool out_of_range = ixstart > ixstop;
	if (!out_of_range)
	{
		/* scan section of grid */
		AutoArrayT<sTYPE>** griddata = fGrid.Pointer(ixstart);
		for (int ix = ixstart; ix <= ixstop; ix++)
		{
			if (*griddata) fHits.Append(**griddata);
			griddata++;
		}
	}
	return fHits;
}	

/* write grid statistics */
template <class sTYPE>
void GridManager1DT<sTYPE>::WriteStatistics(ostream& out) const
{
	/* dimensions */
	out << "\n Search grid statistics:\n";
	out << " Number of grid cells. . . . . . . . . . . . . . = " << fGrid.Length() << '\n';
	out << "     dx = " << fdx << " (" << fnx << ")" <<'\n';

	/* occupancy */
	int min_count =-1;
	int max_count = 0;
	int null_count = 0;
	int tot_count  = 0;
	AutoArrayT<sTYPE>** ppgrid = (AutoArrayT<sTYPE>**) fGrid.Pointer();
	for (int i = 0; i < fGrid.Length(); i++)
	{
		AutoArrayT<sTYPE>* pgrid = *ppgrid++;
	
		if (pgrid)
		{
			int length = pgrid->Length();
			if (length == 0) null_count++;
		
			/* intialized mininum */
			if (min_count == -1) min_count = length;
			
			/* max/min */
			if (length > max_count)
				max_count = length;
			else if (length < min_count)
				min_count = length;
			
			/* total occupancy */
			tot_count += length;
		}
		else
			null_count++;
	}

	/* write results */
	out << " Number of NULL or empty cells . . . . . . . . . = " << null_count << '\n';
	out << " Minimum number of occupants per cell. . . . . . = " << min_count  << '\n';
	out << " Maximum number of occupants per cell. . . . . . = " << max_count  << '\n';
	out << " Average number of occupants per cell. . . . . . = ";
	out << tot_count/fGrid.Length() << '\n';
}

/***********************************************************************
* Protected
***********************************************************************/

/* return pointer to the content list for the given coords */
template <class sTYPE>
AutoArrayT<sTYPE>** GridManager1DT<sTYPE>::FetchGrid(const double* coords)
{
	/* grid indices */
	int ix = int((coords[0] - fxmin)/fdx);
	
	/* range check */
	if (ix < 0 || ix >= fnx) throw ExceptionT::kGeneralFail;		
	
	/* stored column major */
	return fGrid.Pointer(ix);
}

} /* namespace Tahoe */

#endif /* _GRIDMANAGER1D_T_H_ */
