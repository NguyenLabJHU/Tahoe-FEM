/* $Id: GridManager3DT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (12/06/1997)                                          */
/* Interface for regular rectangular search and storage grid              */
/* sTYPE requirements:                                                    */
/* Clear() generates a "blank" data card                                  */
/* Coords() returns a double* to the coordinates                          */
/* 	   operator=                                                          */

#ifndef _GRIDMANAGER3D_T_H_
#define _GRIDMANAGER3D_T_H_

#include "Environment.h"
#include "Constants.h"

/* language support */
#include <iostream.h>

/* direct members */
#include "AutoArrayT.h"
#include "pArrayT.h"
#include "dArray2DT.h"

/* forward declarations */
class iArrayT;

template <class sTYPE>
class GridManager3DT
{
public:

	/* constructors */
	GridManager3DT(double xmin, double xmax, int nx,
	             double ymin, double ymax, int ny,
	             double zmin, double zmax, int nz);
	GridManager3DT(int nx, int ny, int nz, const dArray2DT& coords,
		const iArrayT* nodes_used);
	
	/* destructor */
	~GridManager3DT(void);
	
	/* empty grid */
	void Reset(void);
	void Reset(const dArray2DT& coords, const iArrayT* nodes_used);

	/* insert data into the grid */
	void Add(const sTYPE& data);

	/* closest point - error if grid is empty */
	const sTYPE& Closest(double* target);
	
	/* return list of data falling within the defined region */
	const AutoArrayT<sTYPE>& HitsInRegion(double* coords, double distance);
	const AutoArrayT<sTYPE>& HitsInRegion(double* coords, int cellspan);

	/* the distance covered by the given cell span */
	double CellSpan(int cell_span) const;

	/* grid statistics */
	void WriteStatistics(ostream& out) const;

protected:

	/* return pointer to the content list for the given coords */
	AutoArrayT<sTYPE>** FetchGrid(double* coords);

protected:

	double fxmin;
	double fxmax;
	double fdx;
	int	   fnx;
	
	double fymin;
	double fymax;
	double fdy;
	int	   fny;

	double fzmin;
	double fzmax;
	double fdz;
	int	   fnz;
	
	/* indexing offsets */
	int fxjump;
	int fyjump;

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
GridManager3DT<sTYPE>::GridManager3DT(
double xmin, double xmax, int nx,
	double ymin, double ymax, int ny,
	double zmin, double zmax, int nz):
	fxmin(xmin),
	fxmax(xmax),
	fnx(nx),
	fymin(ymin),
	fymax(ymax),
	fny(ny),
	fzmin(zmin),
	fzmax(zmax),
	fnz(nz),
	fxjump(fny*fnz),
	fyjump(fnz),
	fGrid(fnx*fny*fnz)
{
	/* consistency */
	if (fxmax <= fxmin || fymax <= fymin || fzmax <= fzmin)
		throw eGeneralFail;

	/* grid spacings */
	fdx = (fxmax - fxmin)/fnx;
	fdy = (fymax - fymin)/fny;
	fdz = (fzmax - fzmin)/fnz;
}	

template <class sTYPE>
GridManager3DT<sTYPE>::GridManager3DT(int nx, int ny, int nz,
	const dArray2DT& coords, const iArrayT* nodes_used):
	fnx(nx),
	fny(ny),
	fnz(nz)
{
	/* initialize grid */
	Reset(coords, nodes_used);
}	

template <class sTYPE>
GridManager3DT<sTYPE>::~GridManager3DT(void) { }	

/* empty grid */
template <class sTYPE>
void GridManager3DT<sTYPE>::Reset(void)
{
	AutoArrayT<sTYPE>** pgrid = fGrid.Pointer();

	for (int i = 0; i < fGrid.Length(); i++)
	{
		if (*pgrid) (*pgrid)->Allocate(0);
		pgrid++;
	}

}

template <class sTYPE>
void GridManager3DT<sTYPE>::Reset(const dArray2DT& coords,
	const iArrayT* nodes_used)
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
			fzmin = fzmax = coords(0,2);

			for (int i = 0; i < coords.MajorDim(); i++)
			{
				double* p = coords(i);
		
				if (p[0] < fxmin) fxmin = p[0];
				else if (p[0] > fxmax) fxmax = p[0];
	
				if (p[1] < fymin) fymin = p[1];
				else if (p[1] > fymax) fymax = p[1];
	
				if (p[2] < fzmin) fzmin = p[2];
				else if (p[2] > fzmax) fzmax = p[2];
			}
		}
		else
		{
			fxmin = fxmax = 0.0;
			fymin = fymax = 0.0;
			fzmin = fzmax = 0.0;
		}
	}
	else
	{
		n_max = nodes_used->Length();
		if (nodes_used->Length() > 0)
		{
			int* dex = nodes_used->Pointer();

			/* initialize limits */
			fxmin = fxmax = coords(*dex, 0);
			fymin = fymax = coords(*dex, 1);
			fzmin = fzmax = coords(*dex, 2);

			for (int i = 0; i < nodes_used->Length(); i++)
			{
				double* p = coords(*dex++);
		
				if (p[0] < fxmin) fxmin = p[0];
				else if (p[0] > fxmax) fxmax = p[0];
	
				if (p[1] < fymin) fymin = p[1];
				else if (p[1] > fymax) fymax = p[1];
	
				if (p[2] < fzmin) fzmin = p[2];
				else if (p[2] > fzmax) fzmax = p[2];
			}
		}
		else
		{
			fxmin = fxmax = 0.0;
			fymin = fymax = 0.0;
			fzmin = fzmax = 0.0;
		}
	}
	
	/* min n_max */
	n_max = (n_max < 1) ? 1 : n_max;

	/* add a little space */
	fdx = 0.0001*(fxmax - fxmin);
	fdy = 0.0001*(fymax - fymin);
	fdz = 0.0001*(fzmax - fzmin);

	/* check for zero width */
	fdx = (fdx < kSmall) ? 1.0 : fdx;
	fdy = (fdy < kSmall) ? 1.0 : fdy;
	fdz = (fdz < kSmall) ? 1.0 : fdz;

	fxmax += fdx;
	fymax += fdy;
	fzmax += fdz;
	
	fxmin -= fdx;
	fymin -= fdy;
	fzmin -= fdz;

	/* grid spacings */
	fdx = (fxmax - fxmin)/fnx;
	fdy = (fymax - fymin)/fny;
	fdz = (fzmax - fzmin)/fnz;

	/* limit grid stretch */
	double stretch_limit = 2.0;
	if (fdx < fdy && fdx < fdz)
	{
		if (fdy/fdx > stretch_limit)
		{
			double r = fdy/fdx/stretch_limit;
			fny = int((2.0*r*fny + 1.0)/2);
			fny = (fny > n_max) ? n_max : fny;
			fdy = (fymax - fymin)/fny;
		}
		
		if (fdz/fdx > stretch_limit)
		{
			double r = fdz/fdx/stretch_limit;
			fnz = int((2.0*r*fnz + 1.0)/2);
			fnz = (fnz > n_max) ? n_max : fnz;
			fdz = (fzmax - fzmin)/fnz;
		}
	}
	else if (fdy < fdx && fdy < fdz)
	{
		if (fdx/fdy > stretch_limit)
		{
			double r = fdx/fdy/stretch_limit;
			fnx = int((2.0*r*fnx + 1.0)/2);
			fnx = (fnx > n_max) ? n_max : fnx;
			fdx = (fxmax - fxmin)/fnx;
		}
		
		if (fdz/fdy > stretch_limit)
		{
			double r = fdz/fdy/stretch_limit;
			fnz = int((2.0*r*fnz + 1.0)/2);
			fnz = (fnz > n_max) ? n_max : fnz;
			fdz = (fzmax - fzmin)/fnz;
		}
	}
	else
	{
		if (fdx/fdz > stretch_limit)
		{
			double r = fdx/fdz/stretch_limit;
			fnx = int((2.0*r*fnx + 1.0)/2);
			fnx = (fnx > n_max) ? n_max : fnx;
			fdx = (fxmax - fxmin)/fnx;
		}
		
		if (fdy/fdz > stretch_limit)
		{
			double r = fdy/fdz/stretch_limit;
			fny = int((2.0*r*fny + 1.0)/2);
			fny = (fny > n_max) ? n_max : fny;
			fdy = (fymax - fymin)/fny;
		}
	}
	
	/* set grid parameters */
	fxjump = fny*fnz;
	fyjump = fnz;
	fGrid.Allocate(fnx*fny*fnz);
}

/* insert data into the grid */
template <class sTYPE>
void GridManager3DT<sTYPE>::Add(const sTYPE& data)
{
	/* get content list */
	AutoArrayT<sTYPE>** griddata = FetchGrid(data.Coords());

	/* initialize new list */
	if (!(*griddata))
	{
		*griddata = new AutoArrayT<sTYPE>;
		if (!*griddata) throw(eOutOfMemory);
	}

	/* append value */
	(*griddata)->Append(data);	
}

/* closest point */
template <class sTYPE>
const sTYPE& GridManager3DT<sTYPE>::Closest(double* target)
{
	/* broaden search */
	double distance = (fdx + fdy + fdz)/12.0; //so initial is 1/2 of average
	do {
	
		distance *= 2;
		HitsInRegion(target, distance);
	
	} while (fHits.Length() == 0 && (target[0] - distance > fxmin ||
	                                 target[0] + distance < fxmax ||
	                                 target[1] - distance > fymin ||
	                                 target[1] + distance < fymax ||
	                                 target[2] - distance > fzmin ||
	                                 target[2] + distance < fzmax) );

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
	double dz = coords[2] - target[2];

	double minsqr = dx*dx + dy*dy + dz*dz;
	int    mindex = 0;
	for (int i = 1; i < fHits.Length(); i++)
	{
		coords = (++pdata)->Coords();
		dx = coords[0] - target[0];
		dy = coords[1] - target[1];
		dz = coords[2] - target[2];
		double dsqr = dx*dx + dy*dy + dz*dz;
	
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
inline double GridManager3DT<sTYPE>::CellSpan(int cell_span) const
{
	/* max cell dimension */
	double dmax = (fdx > fdy) ?
		((fdx > fdz) ? fdx : fdz) :
		((fdy > fdz) ? fdy : fdz);
	return (cell_span*dmax)/2.0;
}

/* return list of data falling within the defined region */
template <class sTYPE>
inline const AutoArrayT<sTYPE>& GridManager3DT<sTYPE>::
	HitsInRegion(double* coords, int cellspan)
{
	return HitsInRegion(coords, CellSpan(cellspan));
}	

/* return list of data falling within the defined region */
template <class sTYPE>
const AutoArrayT<sTYPE>& GridManager3DT<sTYPE>::
	HitsInRegion(double* coords, double distance)
{
	/* empty hit list */
	fHits.Allocate(0);

	/* grid indices */
	int ixstart = int((coords[0] - fxmin - distance)/fdx);
	int iystart = int((coords[1] - fymin - distance)/fdy);
	int izstart = int((coords[2] - fzmin - distance)/fdz);
	int ixstop  = int((coords[0] - fxmin + distance)/fdx);
	int iystop  = int((coords[1] - fymin + distance)/fdy);
	int izstop  = int((coords[2] - fzmin + distance)/fdz);

	/* keep within grid */
	ixstart = (ixstart < 0) ? 0 : ixstart;
	iystart = (iystart < 0) ? 0 : iystart;
	izstart = (izstart < 0) ? 0 : izstart;
	
	ixstop = (ixstop >= fnx) ? fnx - 1 : ixstop;
	iystop = (iystop >= fny) ? fny - 1 : iystop;
	izstop = (izstop >= fnz) ? fnz - 1 : izstop;

	bool out_of_range = ixstart > ixstop ||
	                    iystart > iystop ||
	                    izstart > izstop;
	if (!out_of_range)
	{
		/* scan section of grid */
		for (int ix = ixstart; ix <= ixstop; ix++)
			for (int iy = iystart; iy <= iystop; iy++)
			{
				/* column top */
				AutoArrayT<sTYPE>** griddata = fGrid.Pointer(ix*fxjump + iy*fyjump + izstart);
			
				/* copy contents from cells */
				for (int iz = izstart; iz <= izstop; iz++)
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
void GridManager3DT<sTYPE>::WriteStatistics(ostream& out) const
{
	/* dimensions */
	out << "\n Search grid statistics:\n";
	out << " Number of grid cells. . . . . . . . . . . . . . = " << fGrid.Length() << '\n';
	out << "     dx = " << fdx << " (" << fnx << ")" <<'\n';
	out << "     dy = " << fdy << " (" << fny << ")" <<'\n';
	out << "     dz = " << fdz << " (" << fnz << ")" <<'\n';

	/* occupancy */
	int min_count =-1;
	int max_count = 0;
	int null_count = 0;
	int tot_count  = 0;
	AutoArrayT<sTYPE>** ppgrid = fGrid.Pointer();
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
	out << " Minimum number of occupants per cell. . . . . . = " << min_count << '\n';
	out << " Maximum number of occupants per cell. . . . . . = " << max_count << '\n';
	out << " Average number of occupants per cell. . . . . . = ";
	out << tot_count/fGrid.Length() << '\n';
}

/***********************************************************************
* Protected
***********************************************************************/

/* return pointer to the content list for the given coords */
template <class sTYPE>
AutoArrayT<sTYPE>** GridManager3DT<sTYPE>::FetchGrid(double* coords)
{
	/* grid indices */
	int ix = int((coords[0] - fxmin)/fdx);
	int iy = int((coords[1] - fymin)/fdy);
	int iz = int((coords[2] - fzmin)/fdz);
	
	/* range check */
	if (ix < 0 || ix >= fnx ||
	    iy < 0 || iy >= fny ||
	    iz < 0 || iz >= fnz)
	{
		int d_width = OutputWidth(cout, coords);
		cout << "\n GridManager3DT<sTYPE>::FetchGrid: coordinates are out of range:\n";
		cout <<   "     x: " << setw(d_width) << coords[0] << '\n';
		cout <<   "     y: " << setw(d_width) << coords[1] << '\n';
		cout <<   "     z: " << setw(d_width) << coords[2] << "\n\n";
		cout <<   " Grid limits:\n";
		cout <<   "     {x_min, x_max} = {"
		     << setw(d_width) << fxmin << ","
		     << setw(d_width) << fxmax << "}\n";	
		cout <<   "     {y_min, y_max} = {"
		     << setw(d_width) << fymin << ","
		     << setw(d_width) << fymax << "}\n";	
		cout <<   "     {z_min, z_max} = {"
		     << setw(d_width) << fzmin << ","
		     << setw(d_width) << fzmax <<  "}" << endl;	

		throw eGeneralFail;
	}
	
	/* stored column major */
	return fGrid.Pointer(ix*fxjump + iy*fyjump + iz);
}

#endif /* _GRIDMANAGER3D_T_H_ */
