/* $Id: Contact3DT.cpp,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (07/17/1999)                                          */

#include "Contact3DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "FEManagerT.h"
#include "eControllerT.h"
#include "NodeManagerT.h"
#include "iGridManager3DT.h"
#include "Vector3T.h"

/* parameters */
const int kNumFacetNodes = 3;
const int kMaxNumGrid    = 50;

/* constructor */
Contact3DT::Contact3DT(FEManagerT& fe_manager):
	ContactT(fe_manager, kNumFacetNodes),
	fGrid3D(NULL)
{
	/* check base class initializations */
	if (fNumSD != 3) throw eGeneralFail;
}

/* destructor */
Contact3DT::~Contact3DT(void) {	delete fGrid3D; }

/***********************************************************************
* Protected
***********************************************************************/

//TEMP - convert all contact surfaces to triangles
void Contact3DT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	/* inherited */
	ContactT::EchoConnectivityData(in, out);

	/* can only handle tri facets */
	for (int i = 0; i < fSurfaces.Length(); i++)
	{	
		/* subdivide quad faces */
		if (fNumFacetNodes == 3 && fSurfaces[i].MinorDim() == 4)
		{
			/* message */
			out << "\n Contact3DT::EchoConnectivityData: subdividing 4-nodes facets on\n";
			out <<   "     surface " << i+1 << " into (2) triangular facets" << endl;

			ConvertQuadToTri(fSurfaces[i]);
		}
	}
}

/* convert quad facets to tri's */
void Contact3DT::ConvertQuadToTri(iArray2DT& surface) const
{
	/* fix empty sets */
	if (surface.MinorDim() == 4)
	{		
		//TEMP - could do this smarter, i.e., choose shorter diagonal. Also,
		//       could cross-cut into 4 facets, but "ghost" node needs
		//       special treatment
		
		/* generate tri's */
		int num_surfaces = surface.MajorDim();
		iArray2DT surf_tmp(2*num_surfaces, 3);
		for (int j = 0; j < num_surfaces; j++)
		{
			int* quad = surface(j);
			int* tri1 = surf_tmp(2*j);
			int* tri2 = surf_tmp(2*j + 1);
		
			*tri1++ = quad[0];
			*tri1++ = quad[1];
			*tri1   = quad[2];
			*tri2++ = quad[0];
			*tri2++ = quad[2];
			*tri2   = quad[3];
		}
		
		/* exchange facet data */
		surf_tmp.Swap(surface);
	}
	else if (surface.MinorDim() != kNumFacetNodes)
	{
		cout << "\n Contact3DT::EchoConnectivityData: only 3- or 4-noded facets\n";
		cout <<   "     are supported" << endl;
		throw eGeneralFail;	
	}
}

/* generate contact element data */
bool Contact3DT::SetActiveInteractions(void)
{
//NOTE - this is very similar to Contact2DT::SetActiveInteractions(), could
//       make search grid the default behavior for all contact

	int last_num_active = fActiveStrikers.Length();

	/* collect current striker node coords */
	if (fStrikerTags.Length() > 0)
		fStrikerCoords.RowCollect(fStrikerTags, fNodes->CurrentCoordinates());
		
	/* construct search grid if needed */
	if (!fGrid3D)
	{
		/* try to get roughly least 10 per grid */
		int ngrid = int(pow(fStrikerCoords.MajorDim()/10.0,
		                    1.0/fStrikerCoords.MinorDim())) + 1;

		ngrid = (ngrid < 2) ? 2 : ngrid;
		ngrid = (ngrid > kMaxNumGrid) ? kMaxNumGrid : ngrid;

		fGrid3D = new iGridManager3DT(ngrid, ngrid, ngrid, fStrikerCoords, 0);
		if (!fGrid3D) throw eOutOfMemory;

		/* search grid statistics */
		ostream& out = fFEManager.Output();
		out << "\n Search grid: group " << fFEManager.ElementGroupNumber(this) + 1 << '\n';
		fGrid3D->WriteStatistics(out);
	}
	
	/* (re-)set grid boundaries */
	fGrid3D->Reset();
		
	/* set striker/facet data */
	SetActiveStrikers();

	/* assume changed unless last and current step have no active */
	if (last_num_active == 0 && fActiveStrikers.Length() == 0)
		return false;
	else
		return true;
}	

/* generate element data (based on current striker/body data) */
void Contact3DT::SetConnectivities(void)
{
	/* check */
	if (fConnectivities.MajorDim() != fActiveStrikers.Length())
	{
		cout << "\n Contact3DT::SetConnectivities: expecting the number of contact\n"
		     <<   "    connectivities " << fConnectivities.MajorDim()
		     << " to equal the number of active strikers "
		     << fActiveStrikers.Length() << endl;
		throw eGeneralFail;
	}

	for (int i = 0; i < fConnectivities.MajorDim(); i++)
	{
		const iArray2DT& surface = fSurfaces[fHitSurface[i]];
		
		int   facet = fHitFacets[i];
		int* pfacet = surface(facet);
		int*  pelem = fConnectivities(i);

		/* all element tags */
		pelem[0] = pfacet[0]; // 1st facet node
		pelem[1] = pfacet[1]; // 2nd facet node
		pelem[2] = pfacet[2]; // 3rd facet node
		pelem[3] = fActiveStrikers[i]; // striker node
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* sets active striker data (based on current bodies data) */
void Contact3DT::SetActiveStrikers(void)
{
	/* clear previous contact config */
	fActiveMap = -1;
	fActiveStrikers.Allocate(0);
	fHitSurface.Allocate(0);
	fHitFacets.Allocate(0);

	/* reference to current coordinates */
	const dArray2DT& allcoords = fNodes->CurrentCoordinates(); //EFFECTIVE_DVA
	
	/* by-striker data */
	int numstrikers = fStrikerTags.Length();
	if (numstrikers == 0) numstrikers = allcoords.MajorDim();
	iArrayT strikerfacet(numstrikers);
	dArrayT strikerdists(numstrikers);

	/* loop over surfaces */
	//dArrayT normal(3);
	AutoArrayT<double> dists;
	Vector3T<double> mid_x;
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		const iArray2DT& surface = fSurfaces[i];
		int numfacets = surface.MajorDim();
		for (int j = 0; j < numfacets; j++)
		{
			/* facet node positions */
			int* pfacet = surface(j);
			allcoords.RowAlias(pfacet[0], fx1);	
			allcoords.RowAlias(pfacet[1], fx2);	
			allcoords.RowAlias(pfacet[2], fx3);
	
			/* facet midpoint coordinates */
			mid_x.Average(fx1.Pointer(), fx2.Pointer(), fx3.Pointer());

			double radius = 1.1*Vector3T<double>::Norm(mid_x, fx1.Pointer())*3.0/2.0;
			const AutoArrayT<iNodeT>& hits = fGrid3D->HitsInRegion(mid_x, radius);
			
			/* set closest, closest point projection over all bodies */	
			for (int k = 0; k < hits.Length(); k++)
			{
				/* possible striker */
				int tag = hits[k].Tag();
				int strikertag = (fStrikerTags.Length() == 0) ?
					tag : fStrikerTags[tag];
				
				/* no self contact (per surface) */
				if (!surface.HasValue(strikertag))
				{
					/* possible striker */
					fStriker.Set(fNumSD, hits[k].Coords());
				
					double h;
					if (Intersect(fx1, fx2, fx3, fStriker, h))
					{
						/* first time to facet */
						if (fActiveMap[tag] == -1)
						{
							fActiveMap[tag] = fHitSurface.Length();

							fActiveStrikers.Append(strikertag);
							fHitSurface.Append(i);
							fHitFacets.Append(j);
							dists.Append(h);
						}
						else /* closer point projection */
						{
							int map = fActiveMap[tag];
							if (fabs(h) < fabs(dists[map]))
							{
								fHitSurface[map] = i;
								fHitFacets[map]  = j;
								dists[map] = h;
							}
						}
					}
				}
			}	
		}
	}
}

bool Contact3DT::Intersect(const dArrayT& x1, const dArrayT& x2,
	const dArrayT& x3, const dArrayT& xs, double& h) const
{
	/* workspace vectors */
	Vector3T<double> a, b, c, n;
	a.Diff(x2.Pointer(), x1.Pointer());
	b.Diff(x3.Pointer(), x1.Pointer());
	c.Diff(xs.Pointer(), x1.Pointer());

	/* facet normal (direction) = a x b */
	n.Cross(a, b);
	double mag = n.Norm();
	n /= mag;
	
	/* height */
	h = Vector3T<double>::Dot(n, c);

	/* no eminent contact */
	double close_by = sqrt(mag)/10.0;
	if (fabs(h) > close_by)
		return false;
	else /* check projection onto facet */
	{
		Vector3T<double> ni;
	
		/* edge 1 */
		ni.Cross(a, c);
		if (Vector3T<double>::Dot(n, ni) < 0.0)
			return false;
		else
		{
			Vector3T<double> xis, edge;
		
			/* edge 2 */
			edge.Diff(x3.Pointer(), x2.Pointer());
			 xis.Diff(xs.Pointer(), x2.Pointer());
			ni.Cross(edge, xis);
			if (Vector3T<double>::Dot(n, ni) < 0.0)
				return false;
			else
			{
				/* edge 3 */
				edge.Diff(x1.Pointer(), x3.Pointer());
				 xis.Diff(xs.Pointer(), x3.Pointer());
				ni.Cross(edge, xis);
				if (Vector3T<double>::Dot(n, ni) < 0.0)
					return false;
				else
					return true;
			}
		}	
	}
}
