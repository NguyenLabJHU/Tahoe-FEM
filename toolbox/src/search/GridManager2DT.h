/* $Id: GridManager2DT.h,v 1.10 2004-03-18 01:15:39 paklein Exp $ */
/* created: paklein (12/06/1997) */
#ifndef _GRIDMANAGER2D_T_H_
#define _GRIDMANAGER2D_T_H_

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

/** interface for regular 2D rectangular search and storage grid */
template <class sTYPE>
class GridManager2DT
{
public:

	/* constructors */
	GridManager2DT(double xmin, double xmax, int nx,
	             double ymin, double ymax, int ny);
	GridManager2DT(int nx, int ny, const dArray2DT& coords,
		const ArrayT<int>* nodes_used);
	
	/* destructor */
	~GridManager2DT(void);
	
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
	
	double fymin;
	double fymax;
	double fdy;
	int	   fny;

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
GridManager2DT<sTYPE>::GridManager2DT(double xmin, double xmax, int nx,
	double ymin, double ymax, int ny):
	fxmin(xmin),
	fxmax(xmax),
	fnx(nx),
	fymin(ymin),
	fymax(ymax),
	fny(ny),
	fGrid(fnx*fny)
{
	/* consistency */
	if (fxmax <= fxmin || fymax <= fymin) throw ExceptionT::kGeneralFail;

	/* grid spacings */
	fdx = (fxmax - fxmin)/fnx;
	fdy = (fymax - fymin)/fny;
}	

template <class sTYPE>
GridManager2DT<sTYPE>::GridManager2DT(int nx, int ny, const dArray2DT& coords,
	const ArrayT<int>* nodes_used):
	fnx(nx),
	fny(ny)
{
	/* initialize grid data */
	Reset(coords, nodes_used);
}	

template <class sTYPE>
GridManager2DT<sTYPE>::~GridManager2DT(void)
{

}	

/* empty grid */
template <class sTYPE>
void GridManager2DT<sTYPE>::Reset(void)
{
	AutoArrayT<sTYPE>** pgrid = fGrid.Pointer();

	for (int i = 0; i < fGrid.Length(); i++)
	{
		if (*pgrid) (*pgrid)->Dimension(0);
		pgrid++;
	}

}

template <class sTYPE>
void GridManager2DT<sTYPE>::Reset(const dArray2DT& coords,
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
			fymin = fymax = coords(0,1);
	
			for (int i = 0; i < coords.MajorDim(); i++)
			{
				const double* p = coords(i);
			
				if (p[0] < fxmin) fxmin = p[0];
				else if (p[0] > fxmax) fxmax = p[0];
		
				if (p[1] < fymin) fymin = p[1];
				else if (p[1] > fymax) fymax = p[1];
			}
		}
		else
		{
			fxmin = fxmax = 0.0;
			fymin = fymax = 0.0;
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
			fymin = fymax = coords(*dex, 1);

			for (int i = 0; i < nodes_used->Length(); i++)
			{
				const double* p = coords(*dex++);
		
				if (p[0] < fxmin) fxmin = p[0];
				else if (p[0] > fxmax) fxmax = p[0];
		
				if (p[1] < fymin) fymin = p[1];
				else if (p[1] > fymax) fymax = p[1];
			}
		}
		else
		{
			fxmin = fxmax = 0.0;
			fymin = fymax = 0.0;
		}		
	}

	/* min n_max */
	n_max = (n_max < 1) ? 1 : n_max;

	/* add a little space */
	fdx = 0.0001*(fxmax - fxmin);
	fdy = 0.0001*(fymax - fymin);
	
	/* check for zero width */
	fdx = (fdx < kSmall) ? 1.0 : fdx;
	fdy = (fdy < kSmall) ? 1.0 : fdy;
	
	fxmax += fdx;
	fymax += fdy;
	fxmin -= fdx;
	fymin -= fdy;

	/* grid spacings */
	fdx = (fxmax - fxmin)/fnx;
	fdy = (fymax - fymin)/fny;

	/* limit grid stretch */
	double stretch_limit = 2.0;
	if (fdx < fdy && fdy/fdx > stretch_limit)
	{
		double r = fdy/fdx/stretch_limit;
		fny = int((2.0*r*fny + 1.0)/2);
		fny = (fny > n_max) ? n_max : fny;
		fdy = (fymax - fymin)/fny;
	}
	else if (fdy < fdx && fdx/fdy > stretch_limit)
	{
		double r = fdx/fdy/stretch_limit;
		fnx = int((2.0*r*fnx + 1.0)/2);
		fnx = (fnx > n_max) ? n_max : fnx;
		fdx = (fxmax - fxmin)/fnx;
	}

	/* set grid parameters */
	fGrid.Dimension(fnx*fny);
}	

/* insert data into the grid */
template <class sTYPE>
void GridManager2DT<sTYPE>::Add(const sTYPE& data)
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
const sTYPE& GridManager2DT<sTYPE>::Closest(double* target)
{
	/* broaden search */
	double distance = (fdx + fdy)/8.0; //so initial is 1/2 of average
	do {
	
		distance *= 2;
		HitsInRegion(target, distance);
	
	} while (fHits.Length() == 0 && (target[0] - distance > fxmin ||
	                                 target[0] + distance < fxmax ||
	                                 target[1] - distance > fymin ||
	                                 target[1] + distance < fymax) );

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
	double dy = coords[1] - target[1];

	double minsqr = dx*dx + dy*dy;
	int    mindex = 0;
	for (int i = 1; i < fHits.Length(); i++)
	{
		coords = (++pdata)->Coords();
		dx = coords[0] - target[0];
		dy = coords[1] - target[1];
		double dsqr = dx*dx + dy*dy;
	
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
inline double GridManager2DT<sTYPE>::CellSpan(int cell_span) const
{
	/* max cell dimension */
	double dmax = (fdx > fdy) ? fdx : fdy;
	return (cell_span*dmax)/2.0;
}

/* return list of data falling within the defined region */
template <class sTYPE>
inline const AutoArrayT<sTYPE>& GridManager2DT<sTYPE>::
	HitsInRegion(const double* coords, int cellspan)
{
	return HitsInRegion(coords, CellSpan(cellspan));
}	

/* return list of data falling within the defined region */
template <class sTYPE>
const AutoArrayT<sTYPE>& GridManager2DT<sTYPE>::
	HitsInRegion(const double* coords, double distance)
{
	/* empty hit list */
	fHits.Dimension(0);

	/* grid indices */
	int ixstart = int((coords[0] - fxmin - distance)/fdx);
	int iystart = int((coords[1] - fymin - distance)/fdy);
	int ixstop  = int((coords[0] - fxmin + distance)/fdx);
	int iystop  = int((coords[1] - fymin + distance)/fdy);

	/* keep within grid */
	ixstart = (ixstart < 0) ? 0 : ixstart;
	iystart = (iystart < 0) ? 0 : iystart;
	ixstop = (ixstop >= fnx) ? fnx - 1 : ixstop;
	iystop = (iystop >= fny) ? fny - 1 : iystop;

	bool out_of_range = ixstart > ixstop ||
	                    iystart > iystop;
	if (!out_of_range)
	{
		/* scan section of grid */
		for (int ix = ixstart; ix <= ixstop; ix++)
		{
			/* column top */
			AutoArrayT<sTYPE>** griddata = fGrid.Pointer(ix*fny + iystart);
		
			/* copy contents from cells */
			for (int iy = iystart; iy <= iystop; iy++)
			{
				if (*griddata) fHits.Append(**griddata);
				griddata++;
			}
		}
	}
	
	return fHits;
}	

template <class sTYPE>
const AutoArrayT<sTYPE>& GridManager2DT<sTYPE>::
	HitsInRegion(const double* coords, const ArrayT<double>& dist_xy)
{
	/* empty hit list */
	fHits.Dimension(0);

	/* grid indices */
	int ixstart = int((coords[0] - fxmin - dist_xy[0])/fdx);
	int iystart = int((coords[1] - fymin - dist_xy[1])/fdy);
	int ixstop  = int((coords[0] - fxmin + dist_xy[0])/fdx);
	int iystop  = int((coords[1] - fymin + dist_xy[1])/fdy);

	/* keep within grid */
	ixstart = (ixstart < 0) ? 0 : ixstart;
	iystart = (iystart < 0) ? 0 : iystart;
	ixstop = (ixstop >= fnx) ? fnx - 1 : ixstop;
	iystop = (iystop >= fny) ? fny - 1 : iystop;

	bool out_of_range = ixstart > ixstop ||
	                    iystart > iystop;
	if (!out_of_range)
	{
		/* scan section of grid */
		for (int ix = ixstart; ix <= ixstop; ix++)
		{
			/* column top */
			AutoArrayT<sTYPE>** griddata = fGrid.Pointer(ix*fny + iystart);
		
			/* copy contents from cells */
			for (int iy = iystart; iy <= iystop; iy++)
			{
				if (*griddata) fHits.Append(**griddata);
				griddata++;
			}
		}
	}
	
	return fHits;
}	

/* write grid statistics */
template <class sTYPE>
void GridManager2DT<sTYPE>::WriteStatistics(ostream& out) const
{
	/* dimensions */
	out << "\n Search grid statistics:\n";
	out << " Number of grid cells. . . . . . . . . . . . . . = " << fGrid.Length() << '\n';
	out << "     dx = " << fdx << " (" << fnx << ")" <<'\n';
	out << "     dy = " << fdy << " (" << fny << ")" <<'\n';

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
AutoArrayT<sTYPE>** GridManager2DT<sTYPE>::FetchGrid(const double* coords)
{
	/* grid indices */
	int ix = int((coords[0] - fxmin)/fdx);
	int iy = int((coords[1] - fymin)/fdy);
	
	/* range check */
	if (ix < 0 || ix >= fnx || iy < 0 || iy >= fny ) {
		const char caller[] = "GridManager2DT<sTYPE>::FetchGrid";
		cout << "\n " << caller << ": point out of range\n"
		     << "  1: (" << fxmin << " < " << coords[0] << " < " << fxmax << ")\n"
		     << "  2: (" << fymin << " < " << coords[1] << " < " << fymax << ")\n" << endl;
		ExceptionT::GeneralFail(caller);
	}
	
	/* stored column major */
	return fGrid.Pointer(ix*fny + iy);
}

} // namespace Tahoe 
#endif /* _GRIDMANAGER2D_T_H_ */
