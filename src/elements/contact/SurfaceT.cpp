/*  $Id: SurfaceT.cpp,v 1.2 2001-04-09 22:28:55 rjones Exp $ */
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

/* parameters */

SurfaceT::SurfaceT(void)
{
}

SurfaceT::~SurfaceT(void)
{
}


void SurfaceT::PrintData(ostream& out)
{
        /* echo data and correct numbering offset */
        out << setw(kIntWidth) << "surface"
            << setw(kIntWidth) << "faces"
            << setw(kIntWidth) << "size" << '\n';
 
	out << setw(kIntWidth) << "X"
	<< setw(kIntWidth) << fFaces.Length() << "\n\n";

	/* set offset for output */
	{

#if 0
// need iArray connectivities
                        fFaces++;
                        fFaces.WriteNumbered(out);
                        fFaces--;
                        out << '\n';
# endif
	}
	out << "\n Surface nodes:\n";
	fNodes++;
	out << fNodes.wrap(8) << '\n';
	fNodes--;

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
	fNumFaces = faces_tmp.MajorDim();

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
        fNodes.Allocate(node_count);
        pcount = counts.Pointer();
        int nsurf_nodes = 0;
        for (int k = 0; k < num_nodes; k++)
	{
                if (*pcount++ > 0)
                {
                        fNodes[nsurf_nodes] = k;
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
        fFaces.Allocate(fNumFaces);
	/* assuming all faces have same code */
	int number_of_face_nodes = num_face_nodes[0];
	GeometryT::CodeT face_geometry_code = geometry_code[0];
	for (int i = 0 ; i < fNumFaces ; i++) 
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

        // allocate space for add'l data members
 	int fNumNodes = fNodes.Length();


}

void SurfaceT::Initialize (const NodeManagerT* node_manager) 
{
	ComputeNeighbors();
	kNodeManager = node_manager;
//fNumSD = fe_manager.NodeManager()->NumSD();
	fNumSD = kNodeManager->NumSD();
        fNormals.Allocate(fNumNodes,fNumSD);
	UpdateConfiguration();
}

void SurfaceT::UpdateConfiguration ()
{
 	/* update current coordinates */ 
	fCoordinates.RowCollect(fNodes,kNodeManager->CurrentCoordinates());
	/* update averaged outward normals */
	ComputeSurfaceNormals();
}

void SurfaceT::ComputeNeighbors (void)
{ // assume vertex nodes are ordered CCW and first in connectivity lists

# if (0)
  switch(NumSD) {
     case 3:
	RaggedArray2DT<int> set_data;
       	set_data.Configure(count);

	/* determine node neighbors CCW ordered */
	/*       and face neighbors */
        //HACK : need type and 2d storage for each node
	iArray2DT ahead;
	iArray2DT behind;
	for (i = 0; i < fNumFaces ; i++) {
		face = fFace[i];
		conn = face.Connectivity();
		for (j = 0; j < nVNodes; j++) {
			curr = conn(j);
			next = conn(mod(j + 1,nVNodes) + 1);
			prev = conn(mod(j - 1,nVNodes) + 1);
			ahead(curr)  = next;
			behind(curr) = prev;
			faces = face
		}
	}
	/* search next and previous to get ring */
	/* condense associated faces for member nodes */
	break;
     case 2:
	break;
     default:
	cout << "\n SurfaceT::ComputeNeighbors, not 2D or 3D geometry \n";
  }
#endif

}

void SurfaceT::ComputeSurfaceNormals(void)
{
	/*compute each normal, add, then normalize*/
}
