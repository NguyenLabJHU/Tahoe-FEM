* $Id: ContactSearchT.cpp,v 1.1 2001-03-22 17:57:42 rjones Exp $ */

#include "ContactSearchT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "FEManagerT.h"
#include "eControllerT.h"
#include "NodeManagerT.h"
#include "iGridManagerT.h"
#include "Vector3T.h"

/* parameters */
const int    kMaxNumGrid    = 50;
const double kFaceTolerance = 1.1;
const int    kNumPerGridCell= 10;

/* constructor */
ContactSearchT::ContactSearchT(FEManagerT& fe_manager):
	fGrid(NULL)
{
}

/* destructor */
ContactSearchT::~ContactSearchT(void) {	delete fGrid; }

/***********************************************************************
* Protected
***********************************************************************/

/* generate contact element data */
bool ContactSearchT::SetInteractions(void)
{
	/* stats on last search */

	/* current coordinates and normals */
	UpdateKinematicData();
	ResetOpposingSurfaceData;
		

	for (i = 0 ; i < fNumSurfaces ; i++) {
//		surf1 = 
		for (i = 0 ; i < fNumSurfaces ; i++) {
//surf2 = 
	/* loop over pairs of surfaces */
	/* construct search grid if needed */
	if (!fGrid)
	{
		/* try to get roughly least 10 per grid */
		int ngrid = int(
		  pow(surface1.fContactNodeCoords.MajorDim()/kNumPerGridCell,
		      1.0/surface1.fContactNodeCoords.MinorDim()) ) + 1;

		ngrid = (ngrid < 2) ? 2 : ngrid;
		ngrid = (ngrid > kMaxNumGrid) ? kMaxNumGrid : ngrid;

		fGrid = new iGridManagerT
                         (ngrid, ngrid, ngrid, fContactNodeCoords, 0);
		if (!fGrid) throw eOutOfMemory;

		/* search grid statistics */
		ostream& out = fFEManager.Output();
		out << "\n ContactSearch grid: group " 
		    << fFEManager.ElementGroupNumber(this) + 1 << '\n';
		fGrid->WriteStatistics(out);
	}
	
	/* (re-)set grid boundaries */
	fGrid->Reset();
		
	/* set node-face data */
	NodeFaceSearch(surface1,surface2,parameters);

	/* assume changed unless last and current step have no active */
	if (last_num_active == 0 && fActiveContactNodes.Length() == 0)
		return false;
	else
		return true;

		}
	}
}	

bool ContactSearchT::UpdateInteractions(void)
{
	UpdateProjection(surf1,surf2,parameters)
}


/***********************************************************************
* Private
***********************************************************************/

/* sets active striker data (based on current bodies data) */
void ContactSearchT::NodeFaceSearch(void)
{
	/* clear previous contact config */
	fActiveMap = -1;
	fActiveContactNodes.Allocate(0);
	fActiveSurface.Allocate(0);
	fActiveFacets.Allocate(0);

	/* reference to current coordinates */
	const dArray2DT& allcoords = fNodes->CurrentCoordinates(); 
	
	int num_contact_nodes = fContactNodeTags.Length();

	/* loop over faces */
	for (int i = 0; i < fSurfaces.Length(); i++)
	{
		const iArray2DT& surface = fSurfaces[i];
		int numfaces = surface.MajorDim();
		for (int j = 0; j < numfaces; j++)
		{
		   int* pface = surface(j);
		   /* tracking of old interactions */
		   if (pface.Has_Interaction() )
		   {
		   } else {
		        /* face centroid*/
			pface.ComputeCentroid(centroid);
			/* face "radius"*/	
			double radius = kFaceTolerance*pface.ComputeRadius();
			const AutoArrayT<iNodeT>& 
			 close_nodes = fGrid->ActivesInRegion(centroid, radius);
			
			for (int k = 0; k < close_nodes.Length(); k++)
			{
				/* possible contact */
				int tag = close_nodes[k].Tag();
				int strikertag = 
				  (fContactNodeTags.Length() == 0) ?
				  tag : fContactNodeTags[tag];
				
				/* no self contact within face */
				if (!surface.HasValue(strikertag))
				{
					/* inside and normals opposed */
					if ( pface.Projection(x0,n,loc_coor,g) )
					{
						/* no previous interaction */
// Active map or Opposing Surface Pointer?????
						if (fActiveMap[tag] == -1)
						{
						   /* create interaction */
						}
						else 
						{
						   int map = fActiveMap[tag];
						   /* keep closest */
						   if (fabs(g) < fabs(gap[map]))
						   {
						   /* create interaction */
						   }
						}
					}
				}
			}	
		   }
		}
	}
}

ContactSearch::UpdateProjection
{
	// if not in old face, look at neighbors
}
