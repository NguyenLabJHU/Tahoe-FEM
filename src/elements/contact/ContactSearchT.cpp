/* $Id: ContactSearchT.cpp,v 1.4 2001-04-19 23:47:01 rjones Exp $ */

#include "ContactSearchT.h"

#include "ContactSurfaceT.h"
#include "iGridManagerT.h"
#include "AutoArrayT.h"

/* parameters */
const int    kMaxNumGrid    = 50;
const double kFaceTolerance = 1.1;
const int    kNumPerGridCell= 10;

/* constructor */
ContactSearchT::ContactSearchT
(ArrayT<ContactSurfaceT>& surfaces,
nMatrixT<dArrayT>& search_parameters):
	fSurfaces(surfaces),
	fSearchParameters(search_parameters),
	fGrid(NULL)
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
	
	/* (re-)set grid boundaries */
	fGrid->Reset();
		
	for (int j = 0 ; j < fSurfaces.Length() ; j++) {
		ContactSurfaceT& surface2  = fSurfaces[j]; // "face" surface
		dArrayT& parameters = fSearchParameters(i,j);
		/* set node-face data */
		NodeFaceSearch(surface1,surface2,parameters);
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
		dArrayT& parameters = fSearchParameters(i,j);
                /* update node-face data */
		UpdateProjection(surface1,surface2,parameters);
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
	/* update surface geometry */
	surface.UpdateConfiguration();

	/* update face normals */
	ArrayT<FaceT*>& faces = surface.Faces();
	for (int j = 0; j < faces.Length(); j++) {
		FaceT* face = faces[j];
		face->FaceNormal();
	}

	/* track previous node-face pairs */
	for (int j = 0; j < fSurfaces.Length(); j++) {
#if 0
// HACK fillin later
	  for (int k = 0; k < fSurfaces.Length(); k++) {
		ContactNodeT node = nodes[k];
		if (node.OpposingSurface()) {
			FaceT* opposing_face = node.OpposingFace() ;
			opposing_face->Projection(node,parameters);
		}
	  }
#endif
	}
  }
}

void ContactSearchT::NodeFaceSearch
(ContactSurfaceT& node_surface,ContactSurfaceT& face_surface, 
dArrayT& parameters)
{
  /* loop over faces */
  ArrayT<FaceT*>& faces = face_surface.Faces();
  for (int i = 0; i < face_surface.NumFaces(); i++) {
	FaceT* face = faces[i];
	/* face centroid*/
        face->ComputeCentroid(*centroid);
        /* face "radius"*/
        radius = kFaceTolerance*(face->ComputeRadius());
	/* get nodes in neighborhood */
        const AutoArrayT<iNodeT>&
          close_nodes = fGrid->HitsInRegion(centroid, radius);
	for (int j = 0; j < close_nodes.Length(); j++) {
#if 0
    		if (!node.OpposingSurface()) {
		/* checks : opposing normals, gap, local coords */
			face.Projection(node,parameters) ;
		}
#endif
	}
  }
}

void ContactSearchT::UpdateProjection
(ContactSurfaceT& node_surface,ContactSurfaceT& face_surface,
dArrayT& parameters)
{
  /* loop over faces */
  for (int i = 0; i < node_surface.NumNodes(); i++) {
#if 0
        pnode = fNode[i];
	opposing_face = node.OpposingFace() ;
	opposing_face->Projection(node,parameters);

	if (!node.OpposingFace() ){
		while (!node.OpposingFace() && <MORE> ){
			neigbor_face = opposing_face.NeighborFaces() ;
			neigbor_face.Projection(node,parameters);
		}
	}
#endif
  }
}
