/*  $Id: ContactSurfaceT.cpp,v 1.5 2001-04-30 19:30:19 rjones Exp $ */
#include "ContactSurfaceT.h"

#include "SurfaceT.h"
#include "ContactNodeT.h"

/* parameters */

ContactSurfaceT::ContactSurfaceT(void)
{
}

ContactSurfaceT::~ContactSurfaceT(void)
{
}

void
ContactSurfaceT::AllocateContactNodes(void)
{
	fContactNodes.Allocate(fGlobalNodes.Length());
	for(int i = 0; i < fContactNodes.Length(); i++){
		fContactNodes[i] = new ContactNodeT(*this,i);
	}
#if 0
	if (friction)
		fPreviousContactPoints.Allocate(fGlobalNodes.Length());
#endif
}

void
ContactSurfaceT::CopyCurrentToPrevious(void)
{
	for (int i = 0 ; i < fContactNodes.Length() ; i++) {
#if 0
		fPreviousContactPoints[i] = fContactPoints[i];
		fContactPoints[i].OpposingSurface() = NULL;
#endif
	}
}

void 
ContactSurfaceT::SetPotentialConnectivity(void)
{
	int i,j,k,count;
	ContactNodeT* node;
	const FaceT* face = NULL;
	iArrayT node_face_counts;
	node_face_counts.Allocate(fContactNodes.Length());
	node_face_counts = 0;

	/* count connectivity */
        for (i = 0; i < fContactNodes.Length(); i++){
          node = fContactNodes[i];
          face = node->OpposingFace();
          /* connectivities for potential interactions, based on search tol */
          if (face) {
	    node_face_counts[i] = 1;
            /* inclusive of opposing face */
            const ArrayT<FaceT*>&  faces 
		= node->OpposingFace()->Neighbors();
            for (j = 0; j < faces.Length() ; j++) {
		face = faces[j]; // this is a cast
		node_face_counts[i] += face->Connectivity().Length();
            }
          }
        }

	/* configure connectivity and equation numbers */
	fConnectivities.Configure(node_face_counts);
	fEqNums.Configure(node_face_counts,fNumSD);

	/* fill connectivity */
        for (i = 0; i < fContactNodes.Length(); i++){
          node = fContactNodes[i];
	  face = node->OpposingFace();
	  count = 0;
	  /* connectivities for potential interactions, based on search tol */
	  if (face) {
	    int* node_face_connectivity = fConnectivities(i);
            node_face_connectivity[count] 
		= fGlobalNodes[i]; // node
	    count++;
	    /* inclusive of opposing face */
	    const iArrayT& global_nodes=node->OpposingSurface()->GlobalNodes();
            const ArrayT<FaceT*>&  faces 
		= node->OpposingFace()->Neighbors();
	    for (j = 0; j < faces.Length() ; j++) {
		face = faces[j] ; // this is cast
                const iArrayT& face_connectvity = face->Connectivity();
                for (int k = 0; k < face_connectvity.Length(); k++ ) {
                  node_face_connectivity[count] 
		    = global_nodes[face_connectvity[k]];// face nodes
		  count++;
                }
	    }
	  }
        }
}
