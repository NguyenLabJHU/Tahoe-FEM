/* $Id: ContactSearchT.cpp,v 1.2 2001-04-09 22:28:55 rjones Exp $ */

#include "ContactSearchT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "FEManagerT.h"
#include "ContactSurfaceT.h"
#include "eControllerT.h"
#include "NodeManagerT.h"
#include "iGridManagerT.h"

/* parameters */
const int    kMaxNumGrid    = 50;
const double kFaceTolerance = 1.1;
const int    kNumPerGridCell= 10;

/* constructor */
ContactSearchT::ContactSearchT
(FEManagerT& fe_manager, 
ArrayT<ContactSurfaceT>& surfaces,
nMatrixT<dArrayT>& search_parameters):
	fGrid(NULL),
	fSurfaces(surfaces),
	fSearchParameters(search_parameters)
{
}

/* destructor */
ContactSearchT::~ContactSearchT(void) {	delete fGrid; }

/* generate contact element data */
bool ContactSearchT::SetInteractions(void)
{
  /* current coordinates and normals */
  /* reset opposing node-face pairs and do tracking */
  Initialize();

  /* loop over pairs of surfaces */
  // HACK should be dependent on parameters
  for (int i = 0 ; i < fSurfaces.Length() ; i++) {
	ContactSurfaceT& surface1 = fSurfaces[i]; // "node" surface
	/* construct search grid */
	/* with roughly 10 nodes per grid cell */
	int ngrid = int( pow(surface1.NumNodes()/kNumPerGridCell,
		      1.0/surface1.NumNodesPerFace()) ) + 1;

	ngrid = (ngrid < 2) ? 2 : ngrid;
	ngrid = (ngrid > kMaxNumGrid) ? kMaxNumGrid : ngrid;
	
	fGrid = new iGridManagerT (ngrid, surface1.Coordinates(), 0);
	if (!fGrid) throw eOutOfMemory;
	
	/* search grid statistics */
//ostream& out = fFEManager.Output();
//out << "\n ContactSearch grid: group " 
//<< fFEManager.ElementGroupNumber(this) + 1 << '\n';
//fGrid->WriteStatistics(out);
	
	/* (re-)set grid boundaries */
	fGrid->Reset();
		
	for (int j = 0 ; j < fSurfaces.Length() ; j++) {
		ContactSurfaceT& surface2  = fSurfaces[j]; // "face" surface
		/* set node-face data */
		NodeFaceSearch(surface1,surface2);
 	}
	delete fGrid;
  }
  return 1; // HACK will be a status code
}


bool ContactSearchT::UpdateInteractions(void)
{
  /* current coordinates and normals */
  for (int i = 0 ; i < fSurfaces.Length() ; i++) {
        ContactSurfaceT& surface = fSurfaces[i]; // "node" surface
	surface.UpdateConfiguration();
  }
                
  /* loop over pairs of surfaces */
  for (int i = 0 ; i < fSurfaces.Length() ; i++) {
        ContactSurfaceT& surface1 = fSurfaces[i]; // "node" surface
        for (int j = 0 ; j < fSurfaces.Length(); j++) {
                ContactSurfaceT& surface2  = fSurfaces[j]; // "face" surface
                /* update node-face data */
		UpdateProjection(surface1,surface2);
        }
  }
  return 1; // HACK will be a status code
}


/***********************************************************************
* Private
***********************************************************************/

void ContactSearchT::Initialize(void)
{
  /* loop over surfaces */
  for (int i = 0 ; i < fSurfaces.Length() ; i++) {
	ContactSurfaceT& surface = fSurfaces[i];
	surface.UpdateConfiguration();
	for (int j = 0; j < fSurfaces.Length(); j++) {
#if 0
// HACK fillin later
		if (node.OpposingSurface()) {
			opposing_face = node.OpposingFace() ;
			opposing_face.Projection();
		} else {
			node.OpposingFace() = NULL;
		}
#endif
	}
  }
}

void ContactSearchT::NodeFaceSearch
(ContactSurfaceT& node_surface,ContactSurfaceT& face_surface)
{
  /* loop over faces */
  for (int i = 0; i < face_surface.NumFaces(); i++) {
#if 0
	FaceT* face = Faces[i];
	/* face centroid*/
        face.ComputeCentroid(centroid);
        /* face "radius"*/
        double radius = kFaceTolerance*pface.ComputeRadius();
        const AutoArrayT<iNodeT>&
	/* get nodes in neighborhood */
        close_nodes = fGrid->ActivesInRegion(centroid, radius);
	for (int j = 0; j < close_nodes.Length(); j++) {
		dArrayT& paramters = fSearchParameters(i,j);
    		/* taking first interaction defined */
    		if (!node.OpposingSurface()) {
			/* projection checks : */ 
			/* opposing normals, gap, local coords */
			pface.Projection(x0,n,loc_coor,g, parameters) ;
		}
	}
#endif
  }
}

void ContactSearchT::UpdateProjection
(ContactSurfaceT& node_surface,ContactSurfaceT& face_surface)
{
  /* loop over faces */
  for (int i = 0; i < node_surface.NumNodes(); i++) {
#if 0
        pnode = fNode[i];
	opposing_face = node.OpposingFace() ;
	opposing_face.Projection();

	if (!node.OpposingFace() ){
		while (!node.OpposingFace() && <MORE> ){
			neigbor_face = opposing_face.NeighborFaces() ;
			neigbor_face.Projection();
		}
	}
#endif
  }
}
