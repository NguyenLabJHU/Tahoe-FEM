/*  $Id: ContactSurfaceT.cpp,v 1.10 2001-06-27 18:16:21 rjones Exp $ */
#include "ContactSurfaceT.h"

#include "SurfaceT.h"
#include "ContactNodeT.h"
#include <iostream.h>
#include "ofstreamT.h"

/* parameters */

ContactSurfaceT::ContactSurfaceT(void)
{
}

ContactSurfaceT::~ContactSurfaceT(void)
{
        for (int i=0 ; i < fContactNodes.Length() ; i++) {
                delete fContactNodes[i];
        }

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
//	    node_face_counts[i] = 1;
            /* all nodes in associated primary faces */
            for (j = 0; j <  fNodeNeighbors.MinorDim(i) ; j++) {
                face = fNodeNeighbors(i)[j]; 
                node_face_counts[i] += face->Connectivity().Length();
            }
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
            /* all nodes in associated primary faces */
	    const iArrayT& global_nodes = this->GlobalNodes();
            for (j = 0; j < fNodeNeighbors.MinorDim(i) ; j++) {
                face = fNodeNeighbors(i)[j]; 
                const iArrayT& face_connectivity = face->Connectivity();
                for (k = 0; k < face_connectivity.Length(); k++ ) {
                  node_face_connectivity[count] 
		    = global_nodes[face_connectivity[k]];// face nodes
		  count++;
		}
            }
//            node_face_connectivity[count] 
//		= fGlobalNodes[i]; // node
//	    count++;
	    /* inclusive of opposing face */
	    const iArrayT& opp_global_nodes
		= node->OpposingSurface()->GlobalNodes();
            const ArrayT<FaceT*>&  faces 
		= node->OpposingFace()->Neighbors();
	    for (j = 0; j < faces.Length() ; j++) {
		face = faces[j] ; // this is cast
                const iArrayT& face_connectivity = face->Connectivity();
                for (k = 0; k < face_connectivity.Length(); k++ ) {
                  node_face_connectivity[count] 
		    = opp_global_nodes[face_connectivity[k]];// face nodes
		  count++;
                }
	    }
	  }
        }
}

void
ContactSurfaceT::PrintContactArea(ostream& out) const
{
  	dArrayT weights(this->NumNodesPerFace());	
	dArray2DT points(this->NumNodesPerFace(),fNumSD);

	double total_area = 0.0, contact_area = 0.0;
        for (int f = 0;  f < fFaces.Length(); f++) {
          const FaceT* face = fFaces[f];
          face->Quadrature(points,weights);
          for (int i = 0 ; i < weights.Length() ; i++) {
             total_area += weights[i];
	     // should be toleranced
             if (fContactNodes[face->Node(i)]->Gap() < 0) { 
                contact_area += weights[i];
             }
          } 
	}
	out << "Surface " << this->Tag() << ":" ;
	out << " total area "<< total_area 
	    << " contact area " << contact_area << '\n';
}

void
ContactSurfaceT::PrintGap(ostream& out) const
{
	out << "#Surface " << this->Tag() << '\n';

	for (int n = 0 ; n < fContactNodes.Length(); n++) {
	    // HACK this tolerance should agree with the one in Projection
	    // for 2d these should be sorted
	    if (fContactNodes[n]->Gap() < 1.e7) {
		for (int i = 0; i < fNumSD; i++) {
			out << fContactNodes[n]->Position()[i] << " ";
		}
		out << fContactNodes[n]->Gap() << '\n';
	    }
	}
}

void
ContactSurfaceT::PrintGap(ofstream& out) const
{
        out << "#Surface " << this->Tag() << '\n';

        for (int n = 0 ; n < fContactNodes.Length(); n++) {
            // HACK this tolerance should agree with the one in Projection
            // for 2d these should be sorted
            if (fContactNodes[n]->Gap() < 1.e7) {
                for (int i = 0; i < fNumSD; i++) {
                        out << fContactNodes[n]->Position()[i] << " ";
                }
                out << fContactNodes[n]->Gap() << '\n';
	    }
        }
}


void
ContactSurfaceT::PrintNormals(ofstream& out) const
{
        out << "#Surface " << this->Tag() << '\n';

        for (int n = 0 ; n < fContactNodes.Length(); n++) {
                for (int i = 0; i < fNumSD; i++) {
                        out << fContactNodes[n]->Position()[i] << " ";
                }
                for (int i = 0; i < fNumSD; i++) {
                        out << fContactNodes[n]->Normal()[i] << " ";
                }
		out << '\n';
        }
}

