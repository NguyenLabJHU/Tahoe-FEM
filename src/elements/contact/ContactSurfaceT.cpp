/*  $Id: ContactSurfaceT.cpp,v 1.18 2001-09-24 20:37:24 rjones Exp $ */
#include "ContactSurfaceT.h"

#include <iostream.h>
#include "ofstreamT.h"

#include "ContactNodeT.h"
#include "ContactFaceT.h"
#include "ContactLineL2FaceT.h"
#include "ContactLineQ3FaceT.h"
#include "ContactQuadL4FaceT.h"


/* parameters */

ContactSurfaceT::ContactSurfaceT(void):
	fNumPotentialContactNodes(0)
{
}

ContactSurfaceT::~ContactSurfaceT(void)
{
        for (int i=0 ; i < fContactNodes.Length() ; i++) {
                delete fContactNodes[i];
        }

		for (int i=0 ; i < fContactFaces.Length() ; i++) {
                delete fContactFaces[i];
        }
}

void
ContactSurfaceT::Initialize
(const NodeManagerT* node_manager, int num_multipliers)
{

	/* inherited */
	SurfaceT::Initialize(node_manager);

	fNumMultipliers = num_multipliers;

	/* allocate contact nodes */
	fContactNodes.Allocate(fGlobalNodes.Length());
	for(int i = 0; i < fContactNodes.Length(); i++){
		fContactNodes[i] = new ContactNodeT(*this,i);
	}

	/* allocate contact faces */
	fContactFaces.Allocate(fFaces.Length());
	/* attach them to the geometry faces */
    for (int i = 0 ; i < fContactFaces.Length() ; i++)
    {
      switch (fGeometryType)
      {
        case GeometryT::kLine :
          switch (fNumNodesPerFace)
          {
            case 2: fContactFaces[i] = new ContactLineL2FaceT(fFaces[i]);
            break;
            case 3: fContactFaces[i] = new ContactLineQ3FaceT(fFaces[i]);
            break;
            default: throw eGeneralFail;
          }
          break;
        case GeometryT::kTriangle :
          switch (fNumNodesPerFace)
          {
#if 0
			case 3: fContactFaces[i] = new ContactTriaL3FaceT(fFaces[i]);
			break;
#endif
			default: throw eGeneralFail;
          }
          break;
        case GeometryT::kQuadrilateral :
          switch (fNumNodesPerFace)
          {
            case 4: fContactFaces[i] = new ContactQuadL4FaceT(fFaces[i]);
            break;
            default: throw eGeneralFail;
          }
          break;
        default: throw eGeneralFail;
      }
	}

	if (fNumMultipliers) {
		fMultiplierMap.Allocate(fGlobalNodes.Length());
		fLastMultiplierMap.Allocate(fGlobalNodes.Length());
		fRealGhostNodePairs.Allocate(fGlobalNodes.Length(),2);
		/* fill real node column */
		for(int i = 0; i < fContactNodes.Length(); i++){
			fRealGhostNodePairs(i,0) = fGlobalNodes[i];
		}
	}
}

void 
ContactSurfaceT::SetMultiplierConnectivity(void)
{
	ContactFaceT* face = NULL;
	for (int i = 0; i < fContactFaces.Length(); i++){
		face = fContactFaces[i];
		face->SetMultiplierConnectivity();
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
            /* (1) all nodes in associated primary faces */
            for (j = 0; j <  fNodeNeighbors.MinorDim(i) ; j++) {
                face = fNodeNeighbors(i)[j]; 
                node_face_counts[i] += face->Connectivity().Length();
            }
            /* (2) all nodes in opposing neighbor faces */
			/* inclusive of opposing face */
			const ArrayT<FaceT*>&  faces 
			= node->OpposingFace()->Neighbors();
            for (j = 0; j < faces.Length() ; j++) {
				face = faces[j]; // this is a cast
				node_face_counts[i] += face->Connectivity().Length();
			}
		}
	}

	if (fNumMultipliers) {
		fConnectivities.Configure(node_face_counts,2);
	 	fEqNums.Configure(node_face_counts,fNumSD+fNumMultipliers);
	} else {
		/* configure connectivity and equation numbers */
		fConnectivities.Configure(node_face_counts);
	 	fEqNums.Configure(node_face_counts,fNumSD);
	}

	/* fill connectivity */
	for (i = 0; i < fContactNodes.Length(); i++){
		node = fContactNodes[i];
		face = node->OpposingFace();
		count = 0;
		/* connectivities for potential interactions, based on search tol */
		if (face) {
			/* if node has opposing face it is potentially in contact */
			int* node_face_connectivity = fConnectivities(i);
            /* (1) all nodes in associated primary faces */
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
            /* (2) all nodes in opposing neighbor faces */
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
	    	if (count != node_face_counts[i]) {
			cout <<"Error in ContactSurface::SetPotentialConnectivities\n";
			cout <<" count " << count <<" "<<  node_face_counts[i] <<'\n';
			throw eGeneralFail;
		    }
			if (fNumMultipliers) {
				int* node_face_connectivity = fConnectivities(i);
				int local_multiplier_node;
                /* all nodes in associated primary faces */
                for (j = 0; j < fNodeNeighbors.MinorDim(i) ; j++) {
                    face = fNodeNeighbors(i)[j];
                    const iArrayT& face_connectivity = face->Connectivity();
                    for (k = 0; k < face_connectivity.Length(); k++ ) {
                        node_face_connectivity[count]
                        = fMultiplierTags[fMultiplierMap[face_connectivity[k]]];
                      	count++;
                    }
                }
                /* all nodes in opposing neighbor faces */
                /* inclusive of opposing face */
                const iArrayT& opp_multiplier_tags
                    = node->OpposingSurface()->MultiplierTags();
                const iArrayT& opp_multiplier_map
                    = node->OpposingSurface()->MultiplierMap();
                const ArrayT<FaceT*>&  faces
                    = node->OpposingFace()->Neighbors();
                for (j = 0; j < faces.Length() ; j++) {
                    face = faces[j] ; // this is cast
                    const iArrayT& face_connectivity = face->Connectivity();
                    for (k = 0; k < face_connectivity.Length(); k++ ) {
                        node_face_connectivity[count]
                        = opp_multiplier_tags[opp_multiplier_map[face_connectivity[k]]];
                      	count++;
                    }
                }
	    	}
		}
	}
}

/* this is for debugging */
bool 
ContactSurfaceT::IsInConnectivity
(int primary_local_node, int secondary_global_node) const
{
	int ln = primary_local_node;
	for (int i = 0 ; i < fConnectivities.MinorDim(ln); i++){
		if (fConnectivities(ln)[i] == secondary_global_node)
			return 1;
	}
	return 0;
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

void
ContactSurfaceT::PrintStatus(ostream& out) const
{
        out << "#Surface " << this->Tag() << '\n';

        for (int n = 0 ; n < fContactNodes.Length(); n++) {
                out << fContactNodes[n]->Tag()<< " ";
                out << " status " << fContactNodes[n]->Status()  << "\n";
        }
}


void
ContactSurfaceT::PrintMultipliers(ostream& out) const
{
        out << "#Surface " << this->Tag() << '\n';

        for (int n = 0 ; n < fMultiplierMap.Length(); n++) {
			int tag = fMultiplierMap[n];
			if (tag > -1) {
                out << n << " ";
                out << " pressure " << fMultiplierValues(tag,0) << "\n";
			}
			else {
                out << fContactNodes[n]->Tag()<< " ";
                out << " no multipliers " << "\n";
			}
        }
}


void
ContactSurfaceT::DetermineMultiplierExtent(void)
{
    ContactNodeT* node = NULL;
    const FaceT* face = NULL;

    for (int i = 0; i < fContactNodes.Length(); i++){
        node = fContactNodes[i];
        /* if node has opposing face it is potentially in contact */
        if (node->OpposingFace()) {
            /* (1) all nodes in associated primary faces */
            for (int j = 0; j < fNodeNeighbors.MinorDim(i) ; j++) {
                face = fNodeNeighbors(i)[j];
                const iArrayT& face_connectivity = face->Connectivity();
                for (int k = 0; k < face_connectivity.Length(); k++ ) {
                    fMultiplierMap[face_connectivity[k]] = 1;
                }
            }
            /* (2) all nodes in opposing neighbor faces */
            /* inclusive of opposing face */
            const ArrayT<FaceT*>&  faces
                = node->OpposingFace()->Neighbors();
			//node->OpposingSurface()->TagMultiplierMap(faces);
			//casting away const-ness
			ContactSurfaceT* opposing_surface = (ContactSurfaceT*) node->OpposingSurface();
			opposing_surface->TagMultiplierMap(faces);
		}
	}
}

void
ContactSurfaceT::TagMultiplierMap(const ArrayT<FaceT*>&  faces)
{
	for (int j = 0; j < faces.Length() ; j++) {
		const iArrayT& face_connectivity = faces[j]->Connectivity();
		for (int k = 0; k < face_connectivity.Length(); k++ ) {
			fMultiplierMap[face_connectivity[k]] = 1;
		}
	}
}

void
ContactSurfaceT::AllocateMultiplierTags(void) 
{
	/* set last muliplier array to local node map and store values */
	fLastMultiplierMap = fMultiplierMap; // these should be "copy"s
	fLastMultiplierValues = fMultiplierValues;

	/* assign active numbering and total */
	int count = 0;
	for (int n = 0 ; n < fContactNodes.Length(); n++) {
		if (fMultiplierMap[n] > -1) { fMultiplierMap[n] = count++;}
	}
	fNumPotentialContactNodes = count;

	/* Allocate space for ghost node tags for multipliers */
	fMultiplierTags.Allocate(fNumPotentialContactNodes);
	
}

void
ContactSurfaceT::ResetMultipliers(dArray2DT& multiplier_values)
{
        /* set last muliplier array to local node map and store values */
        multiplier_values = 0.0;
        for (int i = 0; i < fMultiplierMap.Length(); i++)
        {
                int old_map = fLastMultiplierMap[i];
                int new_map = fMultiplierMap[i];
                if (old_map > -1 && new_map > -1)
                        multiplier_values[new_map] 
			   = fLastMultiplierValues[old_map];
        }

}

void
ContactSurfaceT::MultiplierTags
(const iArrayT& local_nodes, iArrayT& multiplier_tags)
{
	for (int i = 0; i < local_nodes.Length(); i++)
	{
		multiplier_tags[i] 
			= fMultiplierTags[fMultiplierMap[local_nodes[i]]];
	}
}

iArray2DT& 
ContactSurfaceT::RealGhostNodePairs(void)
{ // for ConnectsDOF
	int ghostnode;
	for (int i = 0; i < fRealGhostNodePairs.MajorDim(); i++) {
	     if(fMultiplierMap[i] > -1) {
		fRealGhostNodePairs(i,1) = fMultiplierTags[fMultiplierMap[i]]; 
	     } else {
		fRealGhostNodePairs(i,1) = -1;
	     }
	}
	return  fRealGhostNodePairs;
}
