/*  $Id: SurfaceT.cpp,v 1.10 2001-04-24 00:33:22 rjones Exp $ */
#include "SurfaceT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "IOBaseT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "ExodusT.h"
#include "ModelFileT.h"
#include "ContinuumElementT.h" // For conversion of side sets to facets.
#include "FaceT.h"
#include "LineL2FaceT.h"
#include "QuadL4FaceT.h"
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


void SurfaceT::PrintConnectivityData(ostream& out)
{
	/* surface data */
        /* echo data and correct numbering offset */
	/* nodes */
	out << "\n Surface " << fTag
            << " nodes:" << setw(kIntWidth) << fGlobalNodes.Length() << '\n' ;
	fGlobalNodes++;
	out << fGlobalNodes.wrap(8) << '\n';
	fGlobalNodes--;

	/* face connectivities */
	out << "\n Faces :" << setw(kIntWidth) 
	    << fFaces.Length() 
	    << ", geometry type :" << GeometryType() <<  '\n' ;
	for (int i = 0 ; i < fFaces.Length() ; i++) {
		iArrayT connectivity = fFaces[i]->Connectivity();
		connectivity++;
        	out << connectivity.wrap(8) << '\n';
        	connectivity--;
	}
	out << '\n';

}

void SurfaceT::PrintKinematicData(ostream& out)
{
	out << "\n Surface " << fTag
            << " nodes:" << setw(kIntWidth) << fGlobalNodes.Length() << '\n' ;

	for (int i = 0 ; i < fNormals.MajorDim() ; i++) {
		double* position = fCoordinates(i);
		out << i << " position: " << position[0] << ", "
		                          << position[1] << ", "
		                          << position[2] << '\n';
		double* normal = fNormals(i);
		out << i << " normal  : " << normal[0] << ", "
		                          << normal[1] << ", "
		                          << normal[2] << '\n';
		out << '\n';
	}

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
	int block_ID;
	int input_format = fe_manager.InputFormat();
	switch (input_format)
	{
		case IOBaseT::kExodusII:
		{
			/* ExodusII database info */
			const StringT& file = fe_manager.ModelFile();
			ostream& out = fe_manager.Output();
			ExodusT database(out);
			if (!database.OpenRead(file))
			{
				cout << "\n ContactT::InputSideSets:"
				     << " error opening file: "
		     		     << file << endl;
				throw eGeneralFail;
			}		

			/* read side set info */
			int set_ID;
			in >> set_ID;
			int num_sides = database.NumSidesInSet(set_ID);
			side_set.Allocate(num_sides, 2);
			database.ReadSideSet(set_ID, block_ID, side_set);

			/* correct offset */
			side_set--;

			/* echo dimensions */
			out << " side set ID: " << set_ID << '\n';
			out << "  element ID: " << block_ID << '\n';
			out << "       sides: " << side_set.MajorDim() << '\n';
			break;
		}
		case IOBaseT::kTahoe:
		{	
			/* read */
			int num_faces;
			in >> block_ID >> num_faces;
			if (num_faces < 0) throw eBadInputValue;
			side_set.Allocate(num_faces, 2);
			in >> side_set;
			
			/* correct offset */
			block_ID--;
			side_set--;
			break;
		}
		case IOBaseT::kTahoeII:
		{
			/* open database */
			ModelFileT database;
			if (database.OpenRead(fe_manager.ModelFile()) 
			    != ModelFileT::kOK)
			{
				cout << "\n ContactT::InputSideSets:"
				     << " error opening file: "
		     		     << fe_manager.ModelFile() << endl;
				throw eGeneralFail;
			}		

			/* read side set info */
			int set_ID;
			in >> set_ID;			
			database.GetSideSet(set_ID, block_ID, side_set);
			
			/* correct offset */
			side_set--;

			/* echo dimensions */
			out << " side set ID: " << set_ID << '\n';
			out << "  element ID: " << block_ID << '\n';
			out << "       sides: " << side_set.MajorDim() << '\n';
			break;
		}
		default:		
			cout << "\n ContactT::InputSideSets:"
			     << " input format not supported: ";
			cout << input_format << endl;
			throw eGeneralFail;
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
	int number_of_face_nodes = num_face_nodes[0];
	GeometryT::CodeT face_geometry_code = geometry_code[0];
	for (int i = 0 ; i < num_faces ; i++) 
	{
	  switch (face_geometry_code) //
	  {
		case GeometryT::kLine :
		  switch (number_of_face_nodes)
	 	  { 
			case 2:
			fFaces[i] = 
			  new LineL2FaceT(*this,fCoordinates, 
			  number_of_face_nodes,faces_tmp(i));
			break;
#if 0
			case 3:
			fFaces[i] = 
			  new LineQ3FaceT(*this,fCoordinates, 
			  number_of_face_nodes,faces_tmp(i) );
#endif
			default:               
			cout << "\n SurfaceT::InputSideSets:" 
			     << " no " << face_geometry_code 
			     << number_of_face_nodes << endl;
			throw eGeneralFail;

		  }
#if 0
		case GeometryT::kTriangle :
		  switch (number_of_face_nodes)
	 	  { 
                        case 3:
                        fFaces[i] =
                          new TriaL3FaceT(*this,fCoordinates, 
			  number_of_face_nodes,faces_tmp(i) );
                        break;

                        default:
                        cout << "\n SurfaceT::InputSideSets:"
                             << " no " << face_geometry_code
                             << number_of_face_nodes << endl;
                        throw eGeneralFail;

		  }
#endif
		case GeometryT::kQuadrilateral :
		  switch (number_of_face_nodes)
	 	  { 

			case 4:
			fFaces[i] = 
			  new QuadL4FaceT(*this,fCoordinates, 
			  number_of_face_nodes,faces_tmp(i));
			break;

			default:
			cout << "\n SurfaceT::InputSideSets:"
			     << " no " << face_geometry_code
			     << number_of_face_nodes << endl;
			throw eGeneralFail;                    
		  }
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

  PrintKinematicData(cout);
}

void SurfaceT::ComputeNeighbors (void)
{ 

  switch(fNumSD) {
     case 2:
	ComputeNeighbors2D();
	break;
     case 3:
	ComputeNeighbors3D();
	break;
     default:
	cout << "\n SurfaceT::ComputeNeighbors, not 2D nor 3D geometry \n";
  }

}

void SurfaceT::ComputeNeighbors2D()
{// assume vertex nodes are ordered LR and first in connectivity lists
#if 0
        //ComputeNeighbors2D();
        // need left/right sense
        for (i = 0; i < fFaces.Length() ; i++) {
                face = fFace[i];
                conn = face.Connectivity();
                for (j = 0; j < conn.Length(); j++) {
                        Inface
                        NeighborNode[conn[1],0] = conn[0]; /left -> right
                        NeighborNode[conn[0],1] = conn[1]; /right <- left
                }
        }
        for (i = 0; i < fFaces.Length() ; i++) {
                face = fFace[i];
                conn = face.Connectivity();
                for (j = 0; j < conn.Length(); j++) {
                       FaceNeighbor.Insert() = InFace[conn[j]];
                }
        }


#endif
}

void SurfaceT::ComputeNeighbors3D()
{// assume vertex nodes are ordered CCW and first in connectivity lists
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
	
	fFaceNeighbors.CopyCompressed(faces_next_to_face);

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
		face->NodeNormal(lnn,*normal_i);
		Add(normal,normal_i,normal,fNumSD);
          }
	  Normalize(normal,fNumSD);
	  /* compute tangents */
	  double* tangent1 = fTangent1s(i);
	  double* tangent2 = fTangent2s(i); // 2D ??
	  FaceT* face = fNodeNeighbors(i)[0];
	  face->LocalBasis(normal,tangent1,tangent2);
        }
}
