/* 
 * File : SurfaceT.cpp
 */
#include "SurfaceT.h"
#include "ContactT.h"

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

/* parameters */

SurfaceT::SurfaceT(void)
{
}


void SurfaceT::PrintData(ostream& out)
{
        /* echo data and correct numbering offset */
        out << setw(kIntWidth) << "surface"
            << setw(kIntWidth) << "faces"
            << setw(kIntWidth) << "size" << '\n';
 
	out << setw(kIntWidth) << "X"
	<< setw(kIntWidth) << fFaces.MajorDim()
	<< setw(kIntWidth) << fFaces.MinorDim() << "\n\n";

	/* set offset for output */
	{
                        fFaces++;
                        fFaces.WriteNumbered(out);
                        fFaces--;
                        out << '\n';
	}
	out << "\n Surface nodes:\n";
	fNodes++;
	out << fNodes.wrap(8) << '\n';
	fNodes--;

}


/* surface input functions */
void SurfaceT::InputSideSets 
(const FEManagerT& kFEManager,ifstreamT& in, ostream& out)
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
                (kFEManager.ElementGroup(elem_group));

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
	int input_format = kFEManager.InputFormat();
	switch (input_format)
	{
		case IOBaseT::kExodusII:
		{
			/* ExodusII database info */
			const StringT& file = kFEManager.ModelFile();
			ostream& out = kFEManager.Output();
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
			if (database.OpenRead(kFEManager.ModelFile()) 
			    != ModelFileT::kOK)
			{
				cout << "\n ContactT::InputSideSets:"
				     << " error opening file: "
		     		     << kFEManager.ModelFile() << endl;
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
	
	////???????????????????????
	pelem_group->SurfaceFacets(geometry_code, 
			surface_facets, surface_nodes); 
	/* global node numbers of faces from element group */
	/* allocates to number of nodes per face */
	pelem_group->SideSetToFacets(block_ID, side_set, faces_tmp);

	/* make node list and convert connectivities to local numbering */
        int num_nodes = (kFEManager.NodeManager())->NumNodes();
        iArrayT counts(num_nodes);
        iArrayT global2local(num_nodes);
        counts = 0;

        /* tally occurrences */
        int* pnode  = fFaces.Pointer();
        int  length = fFaces.Length();
        for (int j = 0; j < length; j++)
        {
                counts[*pnode]++;
                pnode++;
        }

        /* count surface nodes */
        int  node_count = 0;
        int* pcount = counts.Pointer();
        for (int j = 0; j < num_nodes; j++)
                if (*pcount++ > 0)
                        node_count++;

        /* collect */
        fNodes.Allocate(node_count);
        pcount = counts.Pointer();
        int nsurf_nodes = 0;
        for (int k = 0; k < num_nodes; k++)
                if (*pcount++ > 0)
                {
                        fNodes[nsurf_nodes] = k;
                        global2local[k] = nsurf_nodes;
                        nsurf_nodes++;
                }
        /* convert connectvities to local numbering */
        for (int k = 0; k < length; k++) {
                fFaces[k] = global2local[fFaces[k]];
        }

	/* create faces */
	switch (geometry_code)
	{
		case GeometryT::kLine :
		  switch (num_nodes)
	 	  { 
			case 2:
			fFace[i] = new LineL2Face (XXX )
		  }
		case GeometryT::kTriangle :
		  switch (num_nodes)
	 	  { 
			case 3:
			fFace[i] = new TriaL3Face (XXX )
		  }
		case GeometryT::kQuadrilateral :
		  switch (num_nodes)
	 	  { 
			case 4:
			fFace[i] = new QuadL4Face (XXX )
		  }
	}

        // allocate space for add'l data members
 	int NumNodes = fNodes.Length();

        fJacobians.Allocate(NumNodes);
	fNumSD = kFEManager.NodeManager()->NumSD();
        fNormals.Allocate(NumNodes,NumSD);

}

void SurfaceT::Initialize (void) 
{
	ComputeNeighbors();
}

void SurfaceT::UpdateConfiguration (void)
{
   // use NodeManager to get current coordinates
   // and the RowCollect function
   fCoordinates.RowCollect(fNodes,CurrentCoordinates);
   ComputeSurfaceNormals();
}

void SurfaceT::ComputeNeighbors (void)
{
	// need to allocate
	// HOW TO DESTINGUISH BETWEEN VERTEX, EDGE, & INTERIOR NODES
	for (i = 0; i < fNumFaces ; i++) {
		nVNodes =  face.NumVertexNodes();
		for (j = 0; j < nVNodes; j++) {
		// treat connectivities as a CCW ring
			curr = conn(j);
			next = conn(mod(j + 1,nVNodes) + 1)
			prev = mod(j - 1,nVNodes) + 1
			fNodeNeighbor(curr) = next;
			fNodeNeighbor(curr) = prev;

		}
	}
	// node neigbors need to be ordered CCW for 3D
	// face neigbor ++ common node

}

