/* $Id: TotalLagrangianCBSurfaceT.cpp,v 1.12 2005-07-13 05:26:42 hspark Exp $ */
#include "TotalLagrangianCBSurfaceT.h"

#include "ModelManagerT.h"
#include "ShapeFunctionT.h"
#include "FCC3D_Surf.h"
#include "FCC3D.h"
#include "MaterialListT.h"

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
	int nsi = shape.FacetShapeFunction(0).NumIP();		// # IPs per surface face (2x2=4 for 2D surface)
	int nfn = shape.FacetShapeFunction(0).NumNodes();	// # nodes on each surface face?

	/* support for the surface model */
	fF_Surf_List.Dimension(nsi);
	
	/* Need to actually place values into fF_Surf_List when testing (identity) */
	for (int i = 0; i < fF_Surf_List.Length(); i++)
		fF_Surf_List[i].Dimension(nsd);
		
	/* DUMMY INITIALIZE fF_Surf_List - SPECIFY DEFORMATION GRADIENT */
	fF_Surf_List[0].Identity();

	/* Back to normal flow */
	fSurfaceCBSupport = new FSMatSupportT(nsd, nsi);
	fSurfaceCBSupport->SetContinuumElement(this);
	fSurfaceCBSupport->SetDeformationGradient(&fF_Surf_List);

	/* hard coded for hex's with faces parallel to the coordinate axes */
	if (GeometryCode() != GeometryT::kHexahedron)
		ExceptionT::GeneralFail(caller, "only implemented for hex elements");
	
	/* Do we need to redefine this in "canonical" normal order? */
	double normals_dat[6*3] = {
        1.0, 0.0, 0.0,
       -1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0,-1.0, 0.0,
        0.0, 0.0, 1.0,
        0.0, 0.0,-1.0
	};
	dArray2DT normals(6, 3, normals_dat);

	/* find parameter list for the bulk material */
	int num_blocks = list.NumLists("large_strain_element_block");
	if (num_blocks > 1)
		ExceptionT::GeneralFail(caller, "expecting only 1 not %d element blocks", num_blocks);
	const ParameterListT& block = list.GetList("large_strain_element_block");
	const ParameterListT& mat_list = block.GetListChoice(*this, "large_strain_material_choice");
	const ArrayT<ParameterListT>& mat = mat_list.Lists();
	const ParameterListT& bulk_params = mat[0];
	if (bulk_params.Name() != "FCC_3D")
		ExceptionT::GeneralFail(caller, "expecting \"FCC_3D\" not \"%s\"", bulk_params.Name().Pointer());
	
	/* initialize surface information & create all possible (6) surface clusters */
	fNormal.Dimension(nfs);
	fSurfaceCB.Dimension(nfs);
	fSurfaceCB = NULL;

	/* get pointer to the bulk model */
	FCC3D* fcc_3D = NULL;
	if (fMaterialList->Length() == 1) {
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[0];
		fcc_3D = TB_DYNAMIC_CAST(FCC3D*, pcont_mat);
		if (!fcc_3D) ExceptionT::GeneralFail(caller, "could not resolve FCC3D material");
	} else ExceptionT::GeneralFail(caller, "expecting 1 not %d materials", fMaterialList->Length());

	/* Update parameter list for FCC3D_Surf to include the surface normal codes */
	for (int i = 0; i < nfs; i++)
	{
		/* face normal */
		fNormal[i].Dimension(nsd);
		fNormal[i] = normals(i);
	
		/* face C-B model */
		fSurfaceCB[i] = new FCC3D_Surf;
		fSurfaceCB[i]->SetFSMatSupport(fSurfaceCBSupport);
		
		/* pass parameters to the surface model, including surface normal code */
		ParameterListT surf_params = bulk_params;
		surf_params.SetName("FCC_3D_Surf");
		surf_params.AddParameter(i, "normal_code");
		surf_params.AddParameter(fcc_3D->NearestNeighbor(), "bulk_nearest_neighbor");

		/* Initialize a different FCC3D_Surf for each different surface normal type (6 total) */
		fSurfaceCB[i]->TakeParameterList(surf_params);
	}

	/* collect surface element information */
	ArrayT<StringT> block_ID;
	ElementBlockIDs(block_ID);
	ModelManagerT& model_manager = ElementSupport().ModelManager();
	model_manager.BoundingElements(block_ID, fSurfaceElements, fSurfaceElementNeighbors);
	
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
				surf_shape.SurfaceJacobian(jacobian, Q);	// Last column of Q is normal vector to surface face
				
				/* match to face normal types, i.e. normals_dat */
				int normal_type = -1;
				for (int k = 0; normal_type == -1 && k < fNormal.Length(); k++)
				{
					if ((Q.DotCol(nsd-1, fNormal[k]) - 1.0) < -kSmall) 
						normal_type = -1;
					else
						normal_type = k;	
				}
				/* no match */
				if (normal_type == -1)
					ExceptionT::GeneralFail(caller, "could not classify normal on face %d of element %d",
				 		j+1, fSurfaceElements[i]+1);

				/* store */
				fSurfaceElementFacesType(i,j) = normal_type;
			}
	}
}

/* TO DO LIST */
// (1) Add reference areas for surface clusters into FCC3D_Surf.cpp - need to do for modulus?
// (2) Check stress/modulus calculations in FCC3D_Surf.cpp, see what modifications need
// to be made for surface cluster calculations
// (3) Add Stress/Modulus function calls in TotalLagrangianCBSurfaceT.cpp