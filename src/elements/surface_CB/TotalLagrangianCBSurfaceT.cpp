/* $Id: TotalLagrangianCBSurfaceT.cpp,v 1.3 2005-06-30 16:52:58 hspark Exp $ */
#include "TotalLagrangianCBSurfaceT.h"

#include "ModelManagerT.h"
#include "ShapeFunctionT.h"
#include "FCC3D_Surf.h"

using namespace Tahoe;

/* constructor */
TotalLagrangianCBSurfaceT::TotalLagrangianCBSurfaceT(const ElementSupportT& support):
	TotalLagrangianT(support),
	fSurfaceCBSupport(NULL)
{
	SetName("total_lagrangian_CBsurface");
}

/* destructor */
TotalLagrangianCBSurfaceT::~TotalLagrangianCBSurfaceT(void)
{
	/* free surface models */
	for (int i = 0; i < fSurfaceCB.Length(); i++)
		delete fSurfaceCB[i];
	delete fSurfaceCBSupport;
}

/* accept parameter list */
void TotalLagrangianCBSurfaceT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "TotalLagrangianCBSurfaceT::TakeParameterList";
	
	/* inherited */
	TotalLagrangianT::TakeParameterList(list);

	/* the shape functions */
	const ShapeFunctionT& shape = ShapeFunction();
	int nsd = shape.NumSD();	// # of spatial dimensions in problem
	int nfs = shape.NumFacets();	// # of total possible surface facets?
	int nsi = shape.FacetShapeFunction(0).NumIP();		// # IPs per what?
	int nfn = shape.FacetShapeFunction(0).NumNodes();	// # nodes on each surface face?

	/* support for the surface model */
	fF_Surf_List.Dimension(nsi);
	for (int i = 0; i < fF_Surf_List.Length(); i++)
		fF_Surf_List[i].Dimension(nsd);
	fSurfaceCBSupport = new FSMatSupportT(nsd, nsi);
	fSurfaceCBSupport->SetContinuumElement(this);
	fSurfaceCBSupport->SetDeformationGradient(&fF_Surf_List);

	/* hard coded for hex's with faces parallel to the coordinate axes */
	if (GeometryCode() != GeometryT::kHexahedron)
		ExceptionT::GeneralFail(caller, "only implemented for hex elements");
	double normals_dat[6*3] = {
        1.0, 0.0, 0.0,
       -1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0,-1.0, 0.0,
        0.0, 0.0, 1.0,
        0.0, 0.0,-1.0,
	};
	dArray2DT normals(6, 3, normals_dat);

	/* initialize surface information create all possible (12) surface clusters */
	fNormal.Dimension(nfs);
	fSurfaceCB.Dimension(nfs);
	fSurfaceCB = NULL;
	for (int i = 0; i < nfs; i++)
	{
		/* face normal */
		fNormal[i].Dimension(nsd);
		fNormal[i] = normals(i);
	
		/* face C-B model */
		fSurfaceCB[i] = new FCC3D_Surf;
		fSurfaceCB[i]->SetFSMatSupport(fSurfaceCBSupport);
		
		/* pass parameters to the surface model */ // initialize surface material here
		//fSurfaceCB[i]->TakeParameterList
	}

	/* collect surface element information */
	ArrayT<StringT> block_ID;
	ElementBlockIDs(block_ID);
	ModelManagerT& model_manager = ElementSupport().ModelManager();
	model_manager.BoundingElements(block_ID, fSurfaceElements, fSurfaceElementNeighbors);
	// One element (0), and all faces on the surface
	
	/* determine normal type of each face */
	dMatrixT Q(nsd);
	dMatrixT jacobian(nsd, nsd-1);
	LocalArrayT face_coords(LocalArrayT::kInitCoords, nfn, nsd);
	iArrayT face_nodes(nfn), face_nodes_index(nfn);
	ElementSupport().RegisterCoordinates(face_coords);
	fSurfaceElementFacesType = fSurfaceElementNeighbors;
	fSurfaceElementFacesType = -1;
	for (int i = 0; i < fSurfaceElements.Length(); i++)
	{
		for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) /* loop over faces */
			if (fSurfaceElementNeighbors(i,j) == -1) /* no neighbor => surface */
			{
				/* face parent domain */
				const ParentDomainT& surf_shape = shape.FacetShapeFunction(j);

				/* collect coordinates of face nodes */
				ElementCardT& element_card = ElementCard(fSurfaceElements[i]);
				shape.NodesOnFacet(j, face_nodes_index);	// fni = 4 nodes of surface face
				face_nodes.Collect(face_nodes_index, element_card.NodesX());
				face_coords.SetLocal(face_nodes);

				/* face normal (using 1st integration point) */
				surf_shape.DomainJacobian(face_coords, 0, jacobian);
				surf_shape.SurfaceJacobian(jacobian, Q);	// Q is matrix of normals and tangents
				
				/* match to face normal types - match to normals_dat.  normal_type between 0 and 5 */
				int normal_type = -1;
				for (int k = 0; normal_type == -1 && k < fNormal.Length(); k++)
				{
					if ((Q.DotCol(nsd-1, fNormal[k]) - 1.0) < -kSmall) 
						normal_type = -1;
					else
						normal_type = k;	
				}
				/* no match */
				//if (normal_type == -1)
				//	ExceptionT::GeneralFail(caller, "could not classify normal on face %d of element %d",
				// 		j+1, fSurfaceElements[i]+1);

				/* store */
				fSurfaceElementFacesType(i,j) = normal_type;
			}
	}
}
