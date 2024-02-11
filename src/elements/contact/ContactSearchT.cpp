/* $Id: ContactSearchT.cpp,v 1.29 2003-12-20 01:22:14 rjones Exp $ */
#include "ContactSearchT.h"

#include "ContactSurfaceT.h"
#include "ContactNodeT.h"
#include "iGridManagerT.h"
#include "AutoArrayT.h"
#include "iNodeT.h"
#include "ContactElementT.h"

using namespace Tahoe;

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
ContactSearchT::~ContactSearchT(void) {	}

/* generate contact element data */
bool ContactSearchT::SetInteractions(void)
{
  /* current coordinates and normals */
  /* reset opposing node-face pairs and do tracking */
  int i,j;
  /* loop over surfaces */
  for (i = 0 ; i < fSurfaces.Length() ; i++) {
	ContactSurfaceT& surface = fSurfaces[i];

	/* update surface geometry */
	surface.UpdateConfiguration();
	}
	
  bool found = 0;
  /* track previous node-face pairs and reset others */
  int tag;
  for (i = 0; i < fSurfaces.Length(); i++) {
	ContactSurfaceT& surface = fSurfaces[i];
	ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	for (j = 0; j < nodes.Length(); j++) {
		ContactNodeT* node = nodes[j];
		node->AssignOriginal();
		node->AssignLastGap();
		const SurfaceT* osurface = node->OpposingSurface();
 		if (osurface) {
			tag = osurface->Tag();	
			const dArrayT& parameters = fSearchParameters(i,tag);
			found = 
			  node->OpposingFace()->Projection(node,parameters);
			if (!found) node->Initialize();
		}
		else {
			node->Initialize();
		}
    }
  }

  /* loop over pairs of surfaces */
  for (i = 0 ; i < fSurfaces.Length() ; i++) {
	ContactSurfaceT& surface1 = fSurfaces[i]; // "node" surface
	/* construct search grid */
	/* with roughly 10 nodes per grid cell */
	int ngrid = int(pow(double(surface1.NumNodes())/kNumPerGridCell,
		      1.0/surface1.NumNodesPerFace()) ) + 1;

	ngrid = (ngrid < 2) ? 2 : ngrid;
	ngrid = (ngrid > kMaxNumGrid) ? kMaxNumGrid : ngrid;
	
	/* construct a search grid */
	const dArray2DT& coordinates = surface1.Coordinates();
	grid_nodes.Dimension(coordinates.MajorDim());
	grid_nodes.SetValueToPosition();
	iArrayT n_grid(surface1.NumSD());
	n_grid = ngrid; 
	n_grid = 1; //HACK<<<<<
	/* not optimized for thin bodies */
	fGrid = new iGridManagerT (n_grid, coordinates, &grid_nodes);
	if (!fGrid) throw ExceptionT::kOutOfMemory;
	
	/* (re-)set grid boundaries */
	fGrid->Reset();
	/*
	cout << "\nContact Search Grid Statistics for surface " << i <<'\n';
	fGrid->WriteStatistics(cout);
	*/
		
	for (j = 0 ; j < fSurfaces.Length() ; j++) {
		ContactSurfaceT& surface2  = fSurfaces[j]; // "face" surface
		const dArrayT& parameters = fSearchParameters(i,j);
		/* set node-face data */
		if (parameters.Length() != 0) 
			NodeFaceSearch(surface1,surface2,parameters);
 	}
	delete fGrid;
  }

  // history
  for (i = 0; i < fSurfaces.Length(); i++) {
	ContactSurfaceT& surface = fSurfaces[i];
	ArrayT<ContactNodeT*>& nodes = surface.ContactNodes();
	for (j = 0; j < nodes.Length(); j++) {
		ContactNodeT* node = nodes[j];
		if (node->OpposingSurface()) { 
			node->AssignOriginal();
			node->AssignLastGap();
		}
    }
  }

  return 1; // implies contact has changed, HACK will be a status code
}


bool ContactSearchT::UpdateInteractions(void)
{
  /* current coordinates and normals */
  for (int i = 0 ; i < fSurfaces.Length() ; i++)
	fSurfaces[i].UpdateConfiguration(); // "node" surface
                
  /* loop over pairs of surfaces */
  bool ok = UpdateProjection();
  return ok; 
}


/***********************************************************************
* Private
***********************************************************************/


void ContactSearchT::NodeFaceSearch
(ContactSurfaceT& node_surface,ContactSurfaceT& face_surface, 
const dArrayT& parameters)
{
  bool found = 0;
  /* loop over faces */
  const ArrayT<FaceT*>&  faces = face_surface.Faces();
  ArrayT<ContactNodeT*>& nodes = node_surface.ContactNodes();
  /* search tolerance */
  double tol_g = parameters[ContactElementT::kGapTol];
  for (int i = 0; i < face_surface.NumFaces(); i++) {
	const FaceT* face = faces[i];
	/* face centroid*/
	face->ComputeCentroid(centroid);
  	/* face "radius"*/
   	radius = kFaceTolerance * (face->ComputeRadius());
	/* get nodes in neighborhood */
	radius = (radius > tol_g) ? radius : tol_g;
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
			found = 0;
			if (osurface) {
			tag = osurface->Tag();
			const dArrayT& parameters = fSearchParameters(i,tag);
			found = 
			  node->OpposingFace()->Projection(node,parameters);
#if 0
			if (!found) cout << "Warning : "
				<< surface.GlobalNodes()[j]
				<< " node does not project in orginal face\n";
#endif
			if (!found) {
				  /* neighborhood face patch */
				  const ArrayT<FaceT*>& neighbor_faces 
					= node->OriginalOpposingFace()->Neighbors();
				  k = 0;
				  /* look in neighborhood of original face */
				  while (!found && k < neighbor_faces.Length() ) {
					found = 
					 neighbor_faces[k]->Projection(node,parameters);
					k++;
				  }
			}
			if (!found) node->ResetStatus();
			if (!found) cout << "Warning: "
				<< surface.GlobalNodes()[j]
				<< " node lost projection during iteration" << "\n";
//				<< " (a smaller step may be req'd) \n";

            }
        }
  }
  return 1;
}
