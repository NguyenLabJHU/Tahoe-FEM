/* $Id: ContactSearchT.cpp,v 1.5 2001-04-23 17:50:26 rjones Exp $ */

#include "ContactSearchT.h"

#include "ContactSurfaceT.h"
#include "iGridManagerT.h"
#include "AutoArrayT.h"
#include "iNodeT.h"

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
  bool ok = UpdateProjection();
  return ok; 
}


/***********************************************************************
* Private
***********************************************************************/

void ContactSearchT::Initialize(void)
{
  int i,j;
  /* loop over surfaces */
  for (i = 0 ; i < fSurfaces.Length() ; i++) {
	ContactSurfaceT& surface = fSurfaces[i];
	/* update surface geometry */
	surface.UpdateConfiguration();

	/* update face normals */
	ArrayT<FaceT*>& faces = surface.Faces();
	for (j = 0; j < faces.Length(); j++) {
		FaceT* face = faces[j];
		face->FaceNormal();
	}
  }

  bool found = 0;
  /* track previous node-face pairs and reset others */
  int tag;
  for (i = 0; i < fSurfaces.Length(); i++) {
	ContactSurfaceT& surface = fSurfaces[i];
	ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	for (j = 0; j < nodes.Length(); j++) {
		ContactNodeT* node = nodes[j];
		const SurfaceT* osurface = node->OpposingSurface();
 		if (osurface) {
			tag = osurface->Tag();	
			dArrayT& parameters = fSearchParameters(i,tag);
			const FaceT* face = node->OpposingFace() ;
			found = face->Projection(node,parameters);
			if (!found) node->ClearOpposing();
		}
		else {
			node->ClearOpposing();
		}
    	}
  }
}

void ContactSearchT::NodeFaceSearch
(ContactSurfaceT& node_surface,ContactSurfaceT& face_surface, 
dArrayT& parameters)
{
  bool found = 0;
  /* loop over faces */
  ArrayT<FaceT*>&        faces = face_surface.Faces();
  ArrayT<ContactNodeT*>& nodes = node_surface.ContactNodes();
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
		ContactNodeT* node = nodes[close_nodes[j].Tag()];	
		/* take first one FOR NOW */
    		if (!(node->OpposingSurface()) ) {
		/* checks : opposing normals, gap, local coords */
			found = face->Projection(node,parameters) ;
		}
	}
  }
}

bool ContactSearchT::UpdateProjection (void)
{
  int i,j,k;
  bool found = 0;
  /* track previous node-face pairs and reset others */
  int tag;
  for (i = 0; i < fSurfaces.Length(); i++) {
        ContactSurfaceT& surface = fSurfaces[i];
        ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
        for (j = 0; j < nodes.Length(); j++) {
                ContactNodeT* node = nodes[j];
                const SurfaceT* osurface = node->OpposingSurface();
                if (osurface) {
                        tag = osurface->Tag();
                        dArrayT& parameters = fSearchParameters(i,tag);
                        const FaceT* face = node->OpposingFace() ;
                        found = face->Projection(node,parameters);
                        if (!found) {
#if 0
			  neigbor_faces = opposing_face->NeighborFaces();
			  while (!found && k < neighbor_faces.Length() ) {
				face = neighbor_faces[k] ;
				found = 
				  neigbor_face->Projection(node,parameters);

			  }
#endif
			}
                }
		found = 0;
        }
  }
  return 1;
}
