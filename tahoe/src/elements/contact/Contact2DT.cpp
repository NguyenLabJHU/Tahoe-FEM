/* $Id: Contact2DT.cpp,v 1.10 2004-07-15 08:26:08 paklein Exp $ */
/* created: paklein (05/26/1999) */
#include "Contact2DT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ofstreamT.h"
#include "eIntegratorT.h"
#include "iGridManager2DT.h"
#include "ElementSupportT.h"

/* parameters */
const int kNumFacetNodes = 2;
const int kMaxNumGrid    = 75;

using namespace Tahoe;

/* constructor */
Contact2DT::Contact2DT(const ElementSupportT& support):
	ContactT(support, kNumFacetNodes),
	fGrid2D(NULL),
	fv1(2),
	fv2(2)
{
	SetName("contact_2D");
}

/* destructor */
Contact2DT::~Contact2DT(void) {	delete fGrid2D; }

/* allocates space and reads connectivity data */
void Contact2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ContactT::TakeParameterList(list);
	
	/* dimension */
	int neq = NumElementNodes()*NumDOF();
	fNEEvec.Dimension(neq);
	fNEEmat.Dimension(neq);
	SetShapeFunctionArrays();	
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* generate contact element data */
bool Contact2DT::SetActiveInteractions(void)
{
//NOTE - this is very similar to Contact3DT::SetActiveInteractions(), could
//       make search grid the default behavior for all contact

	int last_num_active = fActiveStrikers.Length();

	/* collect current striker node coords */
	if (fStrikerTags.Length() > 0)
		fStrikerCoords.RowCollect(fStrikerTags, ElementSupport().CurrentCoordinates());
		
	/* construct search grid if needed */
	if (!fGrid2D)
	{
		/* try to get roughly least 10 per grid */
		int ngrid = int(pow(fStrikerCoords.MajorDim()/10.0,
		                    1.0/fStrikerCoords.MinorDim())) + 1;

		ngrid = (ngrid < 2) ? 2 : ngrid;
		ngrid = (ngrid > kMaxNumGrid) ? kMaxNumGrid : ngrid;

		fGrid2D = new iGridManager2DT(ngrid, ngrid, fStrikerCoords, 0);
		if (!fGrid2D) throw ExceptionT::kOutOfMemory;

		/* search grid statistics */
		ostream& out = ElementSupport().Output();
		out << "\n Search grid: group " << ElementSupport().ElementGroupNumber(this) + 1 << '\n';
		fGrid2D->WriteStatistics(out);
	}
	
	/* (re-)set grid boundaries */
	fGrid2D->Reset();

	/* update by-body stored data */
	SetSurfacesData();
		
	/* set striker/facet data */
	SetActiveStrikers();
	
	/* assume changed unless last and current step have no active */
	if (last_num_active == 0 && fActiveStrikers.Length() == 0)
		return false;
	else
		return true;
}	

/* sets active striker data (based on current bodies data) */
void Contact2DT::SetActiveStrikers(void)
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
	dArrayT tangent;
	AutoArrayT<double> dists;
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		/*  current surface data */
		dArray2DT& tanvecs = fTanVecs[i];
		dArrayT&   tanmags = fTanMags[i];
		const iArray2DT& surface = fSurfaces[i];
		int numfacets = surface.MajorDim();

		for (int j = 0; j < numfacets; j++)
		{
			/* facet node positions */
			const int* pfacet = surface(j);
			allcoords.RowAlias(pfacet[0], fx1);	
			allcoords.RowAlias(pfacet[1], fx2);	
		
			/* facet tangent vector */
			tanvecs.RowAlias(j, tangent);
	
			/* facet midpoint coordinates */
			fv1.SetToCombination(1.0, fx1, 0.5, tangent);

			double radius = 0.55*tanmags[j];
			const AutoArrayT<iNodeT>& hits =
				fGrid2D->HitsInRegion(fv1.Pointer(), radius);
			
			/* set closest closest point projection over all bodies */	
			for (int k = 0; k < hits.Length(); k++)
			{
				/* possible striker */
				int tag = hits[k].Tag();
				int strikertag = (fStrikerTags.Length() == 0) ?
					tag : fStrikerTags[tag];
				
				/* no self contact (per surface) */
				//if (tag != pfacet[0] && tag != pfacet[1])
				if (!surface.HasValue(strikertag))
				{
					/* possible striker */
					fStriker.Alias(NumSD(), hits[k].Coords());
			
					/* penetration vectors */
					fv1.DiffOf(fStriker, fx1);
					fv2.DiffOf(fStriker, fx2);

					/* distance to facet */
					double   magtan = tanmags[j];				
					double        b = dArrayT::Dot(tangent, fv1)/magtan;
					double    depth = fabs(fv2[0]*fv1[1] - (fv1[0]*fv2[1])/magtan);
					double    max_d = magtan/2.0;
					double overhang = magtan/50.0; // facets a little oversized
					
					/* within cut-off box */
//					if (depth < max_d && (b >= 0.0 && b <= magtan))
					if (depth < max_d && (b+overhang > kSmall && b-overhang < magtan))
//					if (b + overhang > kSmall && b - overhang < magtan)
					{
						/* first time to facet */
						if (fActiveMap[tag] == -1)
						{
							fActiveMap[tag] = fHitSurface.Length();

							fActiveStrikers.Append(strikertag);
							fHitSurface.Append(i);
							fHitFacets.Append(j);
							dists.Append(depth);
						}
						else /* closer point projection */
						{
							int map = fActiveMap[tag];
							if (fabs(depth) < fabs(dists[map]))
							//if (depth < fDists[map])
							{
								fHitSurface[map] = i;
								fHitFacets[map]  = j;
								dists[map] = depth;

							}
						}
					}
				}
			}	
		}
	}
}

/* generate element data (based on current striker/body data) */
void Contact2DT::SetConnectivities(void)
{
	const char caller[] = "Contact2DT::SetConnectivities";

	/* check */
	if (fConnectivities[0]->MajorDim() != fActiveStrikers.Length())
		ExceptionT::GeneralFail(caller, "number of contact connectivities %d != number of active strikers %d",
			fConnectivities[0]->MajorDim(), fActiveStrikers.Length());

	/* set interacting nodes */
	int* pelem = (int*) fConnectivities[0]->Pointer();
	int rowlength = fConnectivities[0]->MinorDim();
	if (fConnectivities[0]->MajorDim() > 0 && rowlength != 3)
		ExceptionT::SizeMismatch(caller, "expecting connectivites length 3 not %d", rowlength);

	for (int i = 0; i < fConnectivities[0]->MajorDim(); i++, pelem += rowlength)
	{
		const iArray2DT& surface = fSurfaces[fHitSurface[i]];
		
		int facet = fHitFacets[i];
		const int* pfacet = surface(facet);

		/* all element tags */
		pelem[0] = pfacet[0]; // 1st facet node
		pelem[1] = pfacet[1]; // 2nd facet node
		pelem[2] = fActiveStrikers[i]; // striker node
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* set working arrays */
void Contact2DT::SetShapeFunctionArrays(void)
{
	/* allocate workspace - displacement DOF's only */
	int neq = NumElementNodes()*NumDOF();
	int nsd = NumSD();
	fdv1T.Dimension(neq, nsd);
	fdv2T.Dimension(neq, nsd);
	fdtanT.Dimension(neq, nsd);

	/* derivative arrays */
	fdv1T = 0.0;
	fdv2T = 0.0;
		
	/* facet nodes */	
	fdv1T(0,0) =-1.0;
	fdv1T(1,1) =-1.0;
	fdv2T(2,0) =-1.0;
	fdv2T(3,1) =-1.0;
	
	/* striker node */	
	fdv1T(4,0) = 1.0;
	fdv1T(5,1) = 1.0;
	fdv2T(4,0) = 1.0;
	fdv2T(5,1) = 1.0;
	
	/* tangent gradient */
	fdtanT = 0.0;
	fdtanT(0,0) =-1.0;
	fdtanT(1,1) =-1.0;
	fdtanT(2,0) = 1.0;
	fdtanT(3,1) = 1.0;
}

/* update by-body stored data */
void Contact2DT::SetSurfacesData(void)
{
	/* allocate by-body work space */
	if (fTanVecs.Length() != fSurfaces.Length())
	{
		int num_surfaces = fSurfaces.Length();
		fTanVecs.Dimension(num_surfaces);
		fTanMags.Dimension(num_surfaces);
		for (int i = 0; i < fSurfaces.Length(); i++)
		{
			int num_facets = fSurfaces[i].MajorDim();	
			fTanVecs[i].Dimension(num_facets, NumSD());
			fTanMags[i].Dimension(num_facets);
		}
	}
	
	/* reference to current coordinates */
	const dArray2DT& allcoords = ElementSupport().CurrentCoordinates();

	/* loop over bodies */
	dArrayT tangent;
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		/* current surface */
		const iArray2DT& surface = fSurfaces[i];
		int numfacets = surface.MajorDim();

		/* loop over facets */
		dArray2DT& tanvecs = fTanVecs[i];
		dArrayT&   magtans = fTanMags[i];
		for (int j = 0; j < numfacets; j++)
		{
			/* facet node positions */
			const int* pfacet = surface(j);
			allcoords.RowAlias(pfacet[0], fx1);	
			allcoords.RowAlias(pfacet[1], fx2);	
				
			/* facet tangent and length */
			tanvecs.RowAlias(j, tangent);
			tangent.DiffOf(fx2, fx1);
			magtans[j] = tangent.Magnitude();
		}
	}
}

/*
NOTES:
(0) It doesn't appear as if the contact surfaces behave consistently
with implicit time integration algorithms. We'll look into
this now.
(1) A striker can interact with only 1 facet per surface, but
	may interact with more than 1 surface at the same time.
(2) No history information is stored, i.e., the complete contact
	configuration for the time interval and the previous values of
	the Lagrange multipliers. With no history information, there's
	no way to step back. Therefore, no step cuts are admissible.
(3) PrintKinematic occurs during Step, i.e., at the start of the
	next time interval, but before the time has been incremented.
	Since the contact configuration is set at the end of the previous
	time interval, the final values of the Lagrange multipliers
	are not available for output. So what is the solution at the
	end of the a time interval: relaxed or pre-relaxation? With
	new equation numbers or before? AFTER seems best. However,
	seems like configuration of the global equation system must
	NOT occur before the end of the current time interval, i.e.,
	some time at the beginning of the next step, but after the
	call to PrintKinematic.

PLANNED:
(1) Allow each node to interact with only a single facet over
	a given time interval. Select this one facet over all surfaces
	by a closest, closest-point projection.
(2) Carry the contact "elements" {node_1, node_2, stiker} into
	the next time interval:
	(i) allows better initialization of the multiplier if a match
		is found.
	(ii) allows reset if convergence fails.
(3) Store all other contact configuration things minimally, i.e.,
	need the multiplier values, but all of the geometrically based
	data should reconstruct identically after stepping back.
*/
