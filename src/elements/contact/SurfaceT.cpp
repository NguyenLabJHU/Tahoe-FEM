/*  $Id: SurfaceT.cpp,v 1.23 2002-03-18 19:24:23 rjones Exp $ */
#include "SurfaceT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "ModelManagerT.h"
#include "fstreamT.h"
#include "IOBaseT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "ContinuumElementT.h" // For conversion of side sets to facets.
#include "FaceT.h"
#include "LineL2FaceT.h"
#include "LineQ3FaceT.h"
#include "QuadL4FaceT.h"
#include "TriaL3FaceT.h"
#include "AutoFill2DT.h"

/* parameters */
const int kHeadroom = 10; // this is a percentage
const int kMaxNumFacesAtNode = 8;
const int kMaxNumFaceNeighbors = 10;
const int kMaxNumNodesAround = kMaxNumFacesAtNode;

/* vector functions */
inline static void Zero(double* v, int dim)
{
	v[0] = 0.0; v[1] = 0.0; if (dim == 3) v[2] = 0.0;
	
}

inline static void Add(double* v1, double* v2, double* v, int dim)
{
	v[0] = v1[0] + v2[0];
	v[1] = v1[1] + v2[1];
	if (dim == 3) v[2] = v1[2] + v2[2];
}


inline static void Normalize(double* v, int dim)
{
	double scale;
        if (dim == 3) {
                scale = 1.0/sqrt(v[0]*v[0] +v[1]*v[1] +v[2]*v[2]);
                v[0] *= scale; v[1] *= scale; v[2] *= scale;
        } 
        else {
                scale = 1.0/sqrt(v[0]*v[0] +v[1]*v[1]);
                v[0] *= scale; v[1] *= scale;
        }
}

SurfaceT::SurfaceT(void)
{
}

SurfaceT::~SurfaceT(void)
{
	for (int i=0 ; i < fFaces.Length() ; i++) {
		delete fFaces[i];
	}
}

void SurfaceT::PrintKinematicData(ostream& out)
{
	out << "\n Surface " << fTag
            << " nodes:" << setw(kIntWidth) << fGlobalNodes.Length() << '\n' ;

	for (int i = 0 ; i < fNormals.MajorDim() ; i++) {
		out << i << " position: " ;
		for (int j = 0 ; j < fNumSD ; j++) {
			out << fCoordinates(i)[j] << ", ";
		}
		out << '\n';
                out << i << " normal  : " ;
                for (int j = 0 ; j < fNumSD ; j++) {
                        out << fNormals(i)[j] << ", ";
                }
                out << '\n';

	}

}

void SurfaceT::PrintConnectivityData(ostream& out)
{ /* surface data */
	/* nodes */
	out << " nodes:" << setw(kIntWidth) << fGlobalNodes.Length() << '\n' ;
	fGlobalNodes++;
	out << fGlobalNodes.wrap(8) << '\n';
	fGlobalNodes--;

	/* face connectivities */
	out << " faces:" << setw(kIntWidth) << fFaces.Length() 
	    << ", geometry type : " << GeometryT::ToString(fGeometryType) <<  '\n' ;
	for (int i = 0 ; i < fFaces.Length() ; i++) {
		iArrayT connectivity = fFaces[i]->Connectivity();
		connectivity++;
        	out << connectivity.wrap(8) << '\n';
        	connectivity--;
	}
	out << '\n';

}



/* surface input functions */
void SurfaceT::InputSideSets 
(const FEManagerT& fe_manager, ifstreamT& in, ostream& out)
{
// parent element determines number of face nodes per face
#ifdef __NO_RTTI__
	cout << "\n ContactT::InputSideSets: RTTI required,"
             << " but not available.\n";
	cout <<   "     Use different surface specification mode." << endl;
	throw;
#endif

	int elem_group;
	in >> elem_group;
	elem_group--;
	ContinuumElementT* pelem_group =
		dynamic_cast<ContinuumElementT*>
                (fe_manager.ElementGroup(elem_group));

	/* checks */
	if (!pelem_group)
	{
		cout << "\n SurfaceT::InputSideSets: element group " 
                     << elem_group;
		cout << " must be of type\n" <<   "     ContinuumElementT" 
		     << endl;
		throw eBadInputValue;
	}

	/* read side set: element, local face pair */
	iArray2DT side_set;
	out <<" Surface: "<< fTag << "\n" ;

	/* read data from parameter file */
	ArrayT<StringT> ss_ID;
	bool multidatabasesets = false; /* change to positive and the parameter file format changes */
	ModelManagerT* model = fe_manager.ModelManager();
	model->SideSetList (in, ss_ID, multidatabasesets);

	if (ss_ID.Length () != 1) 
	  {
	    cout << "\n\nContactT::InputSideSets: Model Manager read more than one side set, not programmed for this.\n\n";
	    throw eBadInputValue;
	  }

	/* read side set */
	StringT block_ID;
	side_set = model->SideSet(ss_ID[0]);
	if (side_set.MajorDim() > 0 && model->IsSideSetLocal(ss_ID[0]))
	    block_ID = model->SideSetGroupID(ss_ID[0]);
	else {
		iArray2DT temp = side_set;
		model->SideSetGlobalToLocal(temp, side_set, block_ID);
	}
	
	/* global node numbers of faces from element group */
	/* allocates to number of nodes per face */
	iArray2DT faces_tmp;
	pelem_group->SideSetToFacets(block_ID, side_set, faces_tmp);
	int num_faces = faces_tmp.MajorDim();

	/* make node list and convert connectivities to local numbering */
        int num_nodes = (fe_manager.NodeManager())->NumNodes();
        iArrayT counts(num_nodes);
        iArrayT global2local(num_nodes);
        counts = 0;

        /* tally occurrences */
        int* pnode  = faces_tmp.Pointer();
        int  length = faces_tmp.Length();
        for (int j = 0; j < length; j++)
        {
                counts[*pnode]++;
                pnode++;
        }

        /* count surface nodes */
        int  node_count = 0;
        int* pcount = counts.Pointer();
        for (int j = 0; j < num_nodes; j++)
	{
                if (*pcount++ > 0) node_count++;
	}

        /* collect */
        fGlobalNodes.Allocate(node_count);
        pcount = counts.Pointer();
        int nsurf_nodes = 0;
        for (int k = 0; k < num_nodes; k++)
	{
                if (*pcount++ > 0)
                {
                        fGlobalNodes[nsurf_nodes] = k;
                        global2local[k] = nsurf_nodes;
                        nsurf_nodes++;
                }
	}
        /* convert connectvities to local numbering */
        for (int k = 0; k < length; k++) 
	{
                faces_tmp[k] = global2local[faces_tmp[k]];
        }

	/* create faces */
	ArrayT <GeometryT::CodeT> geometry_code;
	iArrayT num_face_nodes;
	pelem_group->FacetGeometry(geometry_code, num_face_nodes);
        fFaces.Allocate(num_faces);
	/* assuming all faces have same code */
	fNumNodesPerFace = num_face_nodes[0];
	fGeometryType = geometry_code[0];
	for (int i = 0 ; i < num_faces ; i++) 
	{
	  switch (fGeometryType)
	  {
		case GeometryT::kLine :
		  switch (fNumNodesPerFace)
	 	  { 
			case 2:
			fFaces[i] = 
			  new LineL2FaceT(*this,fCoordinates, 
			  fNumNodesPerFace,faces_tmp(i));
			break;
			case 3:
			fFaces[i] = 
			  new LineQ3FaceT(*this,fCoordinates, 
			  fNumNodesPerFace,faces_tmp(i) );
			break;
			default:               
			cout << "\n SurfaceT::InputSideSets:" 
			     << " no LineFace " << fGeometryType 
			     << " with " << fNumNodesPerFace 
			     << " face nodes \n" ;
			throw eGeneralFail;

		  }
		  break;
		case GeometryT::kTriangle :
		  switch (fNumNodesPerFace)
	 	  { 
			case 3:
			fFaces[i] =
			  new TriaL3FaceT(*this,fCoordinates, 
			  fNumNodesPerFace,faces_tmp(i) );
			break;

			default:
			cout << "\n SurfaceT::InputSideSets:"
		   	     << " no TriangleFace " << fGeometryType
			     << " with " << fNumNodesPerFace 
    			 << " face nodes \n" ;
			
			throw eGeneralFail;

		  }
		  break;
		case GeometryT::kQuadrilateral :
		  switch (fNumNodesPerFace)
	 	  { 

			case 4:
			fFaces[i] = 
			  new QuadL4FaceT(*this,fCoordinates, 
			  fNumNodesPerFace,faces_tmp(i));
			break;

			default:
			cout << "\n SurfaceT::InputSideSets:"
			     << " no QuadFace " << fGeometryType
                             << " with " << fNumNodesPerFace 
                             << " face nodes \n" ;

			throw eGeneralFail;                    
		  }
		  break;
		default:
		   cout << "\n SurfaceT::InputSideSets:"
			<< " unknown face type " << fGeometryType <<"\n";
		   throw eGeneralFail;                    
	  }
	}

}

void SurfaceT::Initialize (const NodeManagerT* node_manager) 
{

	kNodeManager = node_manager;
	fNumSD = kNodeManager->NumSD();

	int num_nodes = fGlobalNodes.Length();
        fCoordinates.Allocate(num_nodes,fNumSD);
        fNormals.Allocate(num_nodes,fNumSD);
        fTangent1s.Allocate(num_nodes,fNumSD);
        if (fNumSD == 3) fTangent2s.Allocate(num_nodes,fNumSD);

	/* initialize faces */
	for (int i=0 ; i < fFaces.Length() ; i++) {
		FaceT*  pface = fFaces[i];
 		pface->Initialize();
	}

	ComputeNeighbors();

	UpdateConfiguration();
}

void SurfaceT::UpdateConfiguration ()
{
  /* update current coordinates */ 
  fCoordinates.RowCollect
	(fGlobalNodes,kNodeManager->CurrentCoordinates());

  /* update averaged outward normals (and tangents) */
  ComputeSurfaceBasis();

  /* update face normals */
  for (int j = 0; j < fFaces.Length(); j++) { fFaces[j]->CalcFaceNormal(); }


#if 0
  PrintKinematicData(cout);
#endif
}

void SurfaceT::ComputeNeighbors(void)
{// 3D: assume vertex nodes are ordered CCW and first in connectivity lists
 // 2D: assume vertex nodes are ordered LR  and first in connectivity lists
	int i,j,k;
	int num_surf_nodes =  fGlobalNodes.Length();
	int num_faces      =  fFaces.Length();
	/* find all faces a node is connected to */
        AutoFill2DT<FaceT*> 
		faces_at_node(num_surf_nodes,kHeadroom,kMaxNumFacesAtNode);
        for (i = 0; i < fFaces.Length() ; i++) {
                FaceT* face = fFaces[i];
                const iArrayT& conn = face->Connectivity();
                for (j = 0; j < conn.Length(); j++) { 
                        faces_at_node.Append(conn[j],face);
                }

        }

	fNodeNeighbors.CopyCompressed(faces_at_node);

	iArrayT face_counts;
	face_counts.Allocate(fNodeNeighbors.MajorDim());
	for (i = 0; i < fNodeNeighbors.MajorDim() ; i++) {
		face_counts[i] = fNodeNeighbors.MinorDim(i);
	}

	/* find the local node number in each neighbor face */
	fLocalNodeInNeighbors.Configure(face_counts);
	for (i = 0; i < fLocalNodeInNeighbors.MajorDim() ; i++) {
		for (j = 0; j < fLocalNodeInNeighbors.MinorDim(i) ; j++) {
		  FaceT* face = fNodeNeighbors(i)[j];
		  fLocalNodeInNeighbors(i)[j] = face->LocalNodeNumber(i);
		}
	}


	/* find all neighbor faces inclusive of orig. face currently */
	AutoFill2DT<FaceT*>
		faces_next_to_face(num_faces,kHeadroom,kMaxNumFaceNeighbors);
	int node_num;
        for (i = 0; i < fFaces.Length() ; i++) {
                FaceT* face = fFaces[i];
                const iArrayT& conn = face->Connectivity();
                for (j = 0; j < conn.Length(); j++) {
		  node_num = conn[j];
		  for (k = 0; k < faces_at_node.MinorDim(node_num); k++) {
			FaceT* neighbor_face = faces_at_node(node_num,k);
                        faces_next_to_face.AppendUnique(i,neighbor_face);
		  }
                }
        }
	
//fFaceNeighbors.CopyCompressed(faces_next_to_face);
	/* copy neighbors faces to face data */
	int num_neighbors;
	for (i = 0; i < fFaces.Length() ; i++) {
		ArrayT<FaceT*>& neighbor_faces = fFaces[i]->AssignNeighbors();
		num_neighbors = faces_next_to_face.MinorDim(i);
		neighbor_faces.Allocate(num_neighbors);
		for (j = 0; j < num_neighbors ; j++) {
		  neighbor_faces[j] = faces_next_to_face(i,j);
		}
	}

}

void SurfaceT::ComputeSurfaceBasis(void)
{ 
	double normal_i[3];
	for (int i = 0; i < fNodeNeighbors.MajorDim(); i++) {
		double* normal = fNormals(i);
		Zero(normal,fNumSD);
  		for (int j = 0; j < fNodeNeighbors.MinorDim(i); j++) {
			FaceT* face = fNodeNeighbors(i)[j];
			int lnn = fLocalNodeInNeighbors(i)[j];
			face->NodeNormal(lnn,normal_i);
			Add(normal,normal_i,normal,fNumSD);
		}
		Normalize(normal,fNumSD);
		/* compute tangents */
		double* tangent1 = fTangent1s(i);
		double* tangent2 = (fNumSD ==3) ? fTangent2s(i) : NULL; 
		FaceT* face = fNodeNeighbors(i)[0];
		face->LocalBasis(normal,tangent1,tangent2);
  	}
}
