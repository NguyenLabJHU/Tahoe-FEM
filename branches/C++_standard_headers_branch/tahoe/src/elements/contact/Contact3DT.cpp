/* $Id: Contact3DT.cpp,v 1.11 2004-07-15 08:26:08 paklein Exp $ */
/* created: paklein (07/17/1999) */
#include "Contact3DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ofstreamT.h"
#include "eIntegratorT.h"
#include "iGridManager3DT.h"
#include "Vector3T.h"
#include "ElementSupportT.h"

using namespace Tahoe;

/* parameters */
const int kNumFacetNodes = 3;
const int kMaxNumGrid    = 50;

/* debugging */
#undef DEBUG

/* constructor */
Contact3DT::Contact3DT(const ElementSupportT& support):
	ContactT(support, kNumFacetNodes),
	fGrid3D(NULL)
{
	SetName("contact_3D");
}

/* destructor */
Contact3DT::~Contact3DT(void) {	delete fGrid3D; }

/***********************************************************************
 * Protected
 ***********************************************************************/

/* Converts quadrilateral faces to triangular faces */
void Contact3DT::ExtractContactGeometry(const ParameterListT& list)
{
	/* inherited */
	ContactT::ExtractContactGeometry(list);

	/* can only handle tri facets */
	for (int i = 0; i < fSurfaces.Length(); i++)
	{	
		/* subdivide quad faces */
		if (fNumFacetNodes == 3 && fSurfaces[i].MinorDim() == 4)
		{
			/* message */
			if (ElementSupport().PrintInput()) {
				ofstreamT& out = ElementSupport().Output();
				out << "\n Contact3DT::EchoConnectivityData: subdividing 4-nodes facets on\n";
				out <<   "     surface " << i+1 << " into (2) triangular facets" << endl;
			}

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
		throw ExceptionT::kGeneralFail;	
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
		fStrikerCoords.RowCollect(fStrikerTags, ElementSupport().CurrentCoordinates());
		
	/* construct search grid if needed */
	if (!fGrid3D)
	{
		/* try to get roughly least 10 per grid */
		int ngrid = int(pow(fStrikerCoords.MajorDim()/10.0,
		                    1.0/fStrikerCoords.MinorDim())) + 1;

		ngrid = (ngrid < 2) ? 2 : ngrid;
		ngrid = (ngrid > kMaxNumGrid) ? kMaxNumGrid : ngrid;

		fGrid3D = new iGridManager3DT(ngrid, ngrid, ngrid, fStrikerCoords, 0);
		if (!fGrid3D) throw ExceptionT::kOutOfMemory;

		/* search grid statistics */
		ostream& out = ElementSupport().Output();
		out << "\n Search grid: group " << ElementSupport().ElementGroupNumber(this) + 1 << '\n';
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
	const char caller[] = "Contact3DT::SetConnectivities";

	/* check */
	if (fConnectivities[0]->MajorDim() != fActiveStrikers.Length())
		ExceptionT::SizeMismatch(caller, "contact interactions %d != active strikers %d",
			fConnectivities[0]->MajorDim(), fActiveStrikers.Length());

	int* pelem = (int*) fConnectivities[0]->Pointer();
	int rowlength = fConnectivities[0]->MinorDim();
	if (fConnectivities[0]->MajorDim() > 0 && rowlength != 4)
		ExceptionT::SizeMismatch(caller, "expecting connectivites length 4 not %d", rowlength);
	
	for (int i = 0; i < fConnectivities[0]->MajorDim(); i++, pelem += rowlength)
	{
		const iArray2DT& surface = fSurfaces[fHitSurface[i]];
		
		int facet = fHitFacets[i];
		const int* pfacet = surface(facet);

		/* all element tags */
		pelem[0] = pfacet[0]; // 1st facet node
		pelem[1] = pfacet[1]; // 2nd facet node
		pelem[2] = pfacet[2]; // 3rd facet node
		pelem[3] = fActiveStrikers[i]; // striker node
	}
}

/* set surface normal derivative matrix */
void Contact3DT::Set_dn_du(const dArray2DT& curr_coords,
	dMatrixT& dn_du) const
{
	double* p = dn_du.Pointer();
	const double* x1 = curr_coords(0);
	const double* x2 = curr_coords(1);
	const double* x3 = curr_coords(2);

	*p++ = 0;
	*p++ =-x2[2] + x3[2];
	*p++ = x2[1] - x3[1];
	*p++ = x2[2] - x3[2];
	*p++ = 0;
	*p++ =-x2[0] + x3[0];
	*p++ =-x2[1] + x3[1];
	*p++ = x2[0] - x3[0];
	*p++ = 0;
	*p++ = 0;
	*p++ = x1[2] - x3[2];
	*p++ =-x1[1] + x3[1];
	*p++ =-x1[2] + x3[2];
	*p++ = 0;
	*p++ = x1[0] - x3[0];
	*p++ = x1[1] - x3[1];
	*p++ =-x1[0] + x3[0];
	*p++ = 0;
	*p++ = 0;
	*p++ =-x1[2] + x2[2];
	*p++ = x1[1] - x2[1];
	*p++ = x1[2] - x2[2];
	*p++ = 0;
	*p++ =-x1[0] + x2[0];
	*p++ =-x1[1] + x2[1];
	*p++ = x1[0] - x2[0];
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p++ = 0;
	*p   = 0;
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* sets active striker data (based on current bodies data) */
void Contact3DT::SetActiveStrikers(void)
{
	/* clear previous contact config */
	fActiveMap = -1;
	fActiveStrikers.Dimension(0);
	fHitSurface.Dimension(0);
	fHitFacets.Dimension(0);

	/* reference to current coordinates */
	const dArray2DT& allcoords = ElementSupport().CurrentCoordinates(); //EFFECTIVE_DVA
	
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
			const int* pfacet = surface(j);
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

#ifdef DEBUG
/* print information about following striker node */
double follow_tag = 311;
bool print_all = false;
if (print_all || strikertag == follow_tag)
{
	ostream& out = ElementSupport().Output();
	out << "\n stiker: " << strikertag+1 << '\n';
	out << " surf index: " << i+1 << '\n';
	out << " face index: " << j+1 << '\n';
	out << " face nodes: " << pfacet[0]+1 << " " << pfacet[1]+1 << " " << pfacet[2]+1 << '\n';
}
#endif
				/* no self contact (per surface) */
				if (!surface.HasValue(strikertag))
				{
					/* possible striker */
					fStriker.Alias(NumSD(), hits[k].Coords());

					double h;
					if (Intersect(fx1, fx2, fx3, fStriker, h))
					{
#ifdef DEBUG
if (print_all || strikertag == follow_tag)
{
	ostream& out = ElementSupport().Output();
	out << "  x = " << fStriker.no_wrap() << '\n';
	out << " x1 = " << fx1.no_wrap() << '\n';
	out << " x2 = " << fx2.no_wrap() << '\n';
	out << " x3 = " << fx3.no_wrap() << '\n';
	out << " h = " << h << '\n';
}
#endif
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
	Vector3T<double> a, b, c, n, xsp;
	a.Diff(x2.Pointer(), x1.Pointer());
	b.Diff(x3.Pointer(), x1.Pointer());
	c.Diff(xs.Pointer(), x1.Pointer());

	/* facet normal (direction) = a x b */
	n.Cross(a, b);
	double mag = n.Norm();
	n /= mag;
	
	/* height */
	h = Vector3T<double>::Dot(n, c);

	/* projection of striker into the plane of the face */
	xsp.Combine(1.0, xs.Pointer(), -h, n);

	/* tolerances */
	double dist_tol = sqrt(mag)/2.0; /* tolerance on perpendicular distance to the face */
	double area_tol = mag/50.0; /* tolerance on node falling "within" the area of the face */
	if (fabs(h) > dist_tol)
		return false;
	else /* check projection onto facet */
	{
		Vector3T<double> xis, ni;
	
		/* edge 1 */
		xis.Diff(xsp, x1.Pointer());
		ni.Cross(a, xis);
		if (Vector3T<double>::Dot(n, ni) < -area_tol)
			return false;
		else
		{
			Vector3T<double> edge;
		
			/* edge 2 */
			edge.Diff(x3.Pointer(), x2.Pointer());
			 xis.Diff(xsp, x2.Pointer());
			ni.Cross(edge, xis);
			if (Vector3T<double>::Dot(n, ni) < -area_tol)
				return false;
			else
			{
				/* edge 3 */
				edge.Diff(x1.Pointer(), x3.Pointer());
				 xis.Diff(xsp, x3.Pointer());
				ni.Cross(edge, xis);
				if (Vector3T<double>::Dot(n, ni) < -area_tol)
					return false;
				else
					return true;
			}
		}	
	}
}
