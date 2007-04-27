/* $Id: TotalLagrangianCBSurfaceT.cpp,v 1.30 2007-04-27 00:56:08 hspark Exp $ */
#include "TotalLagrangianCBSurfaceT.h"

#include "ModelManagerT.h"
#include "ShapeFunctionT.h"
#include "FCC3D_Surf.h"
#include "FCC3D.h"
#include "EAMFCC3D_surf.h"
#include "EAMFCC3D.h"
#include "EAMFCC3DMatT.h"
#include "EAMFCC3DMatT_surf.h"
#include "MaterialListT.h"
#include "eIntegratorT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "InverseMapT.h"

using namespace Tahoe;

/* vector functions */
inline static void CrossProduct(const double* A, const double* B, double* AxB) {   
	AxB[0] = A[1]*B[2] - A[2]*B[1];
	AxB[1] = A[2]*B[0] - A[0]*B[2];
	AxB[2] = A[0]*B[1] - A[1]*B[0];
};

inline static double Dot(const double* A, const double* B){ 
	return A[0]*B[0] + A[1]*B[1] + A[2]*B[2]; 
};

inline static void Vector(const double* start, const double* end, double* v) {
	v[0] = end[0] - start[0];
	v[1] = end[1] - start[1];
	v[2] = end[2] - start[2];
};

inline static void Scale(double* v, double scale) {
	v[0] *= scale;
	v[1] *= scale;
	v[2] *= scale;
};

inline static void Sum(const double* A, const double* B, double* AB) {
	AB[0] = A[0] + B[0];
	AB[1] = A[1] + B[1];
	AB[2] = A[2] + B[2];
};

/* constructor */
TotalLagrangianCBSurfaceT::TotalLagrangianCBSurfaceT(const ElementSupportT& support):
	TotalLagrangianT(support),
	fSurfaceCBSupport(NULL),
	fSplitInitCoords(LocalArrayT::kInitCoords),
	fIndicator(0),
	fSplitShapes(NULL)
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
	delete fSplitShapes;
}

/* accumulate the residual force on the specified node */
void TotalLagrangianCBSurfaceT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
	/* fRHS != 0 AT BEGINNING - WHY? */
	const char caller[] = "TotalLagrangianCBSurfaceT::AddNodalForce";

	/* bulk forces from inherited function */
	TotalLagrangianT::AddNodalForce(field, node, force);

	/* not my field */
	if (&field != &(Field())) return;

	/* quick exit */
	bool hasnode = false;
	for (int i=0; i < fBlockData.Length() && !hasnode; i++)
		if (fConnectivities[i]->HasValue(node)) hasnode = true;
	if (!hasnode) return;

	/*************** surface force contribution ***************/

	/* temp for nodal force */
	dArrayT nodalforce;

	/* time integration parameters */
	double constKd = 0.0;
	int formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* dimensions */
	const ShapeFunctionT& shape = ShapeFunction();
	int nsd = shape.NumSD();                          // # of spatial dimensions in problem
	int nfs = shape.NumFacets();                      // # of total possible element faces
	int nsi = shape.FacetShapeFunction(0).NumIP();    // # IPs per surface face
	int nfn = shape.FacetShapeFunction(0).NumNodes(); // # nodes on each surface face
	int nen = NumElementNodes();                      // # nodes in bulk element

	/* matrix alias to fNEEvec */
	dMatrixT WP(nsd, fStressStiff.Rows(), fNEEvec.Pointer());
	
	/* loop over surface elements */
	dMatrixT jacobian(nsd, nsd-1);
	LocalArrayT face_coords(LocalArrayT::kInitCoords, nfn, nsd);
	iArrayT face_nodes(nfn), face_nodes_index(nfn);
	ElementSupport().RegisterCoordinates(face_coords);
	dArrayT ip_coords_X(nsd);
	dArrayT ip_coords_Xi(nsd);
	dArrayT Na(nen);
	dArray2DT DNa_X(nsd,nen), DNa_Xi(nsd,nen);
	dMatrixT DXi_DX(nsd);
	dMatrixT F_inv(nsd);
	dMatrixT PK1(nsd), cauchy(nsd);
	double t_surface;
	
	/* TEMPORARY */
	dArrayT force1, force2, reactionforce, force3, totalforce;
	force1.Dimension(24);
	force2.Dimension(24);
	reactionforce.Dimension(3);	
	force3.Dimension(3);
	totalforce.Dimension(24);
	force1 = 0.0;
	force2 = 0.0;
	reactionforce = 0.0;
	force3 = 0.0;
	totalforce = 0.0;
	
	int count = 0;
	
	for (int i = 0; i < fSurfaceElements.Length(); i++)
	{
		/* bulk element information */
		int element = fSurfaceElements[i];
		const ElementCardT& element_card = ElementCard(element);
		fLocInitCoords.SetLocal(element_card.NodesX()); /* reference coordinates over bulk element */
		fLocDisp.SetLocal(element_card.NodesU()); /* displacements over bulk element */
	
		/* integrate surface contribution to nodal forces */
		fRHS = 0.0;
		for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) /* loop over faces */
			if (fSurfaceElementNeighbors(i,j) == -1) /* no neighbor => surface */
			{
				/* face parent domain */
				const ParentDomainT& surf_shape = shape.FacetShapeFunction(j);
			
				/* collect coordinates of face nodes */
				ElementCardT& element_card = ElementCard(fSurfaceElements[i]);
				shape.NodesOnFacet(j, face_nodes_index);	// fni = 4 nodes of surface face
				face_nodes.Collect(face_nodes_index, element_card.NodesX());

				const iArrayT& nodes_u = element_card.NodesU();
				int nodeposition = 0;
//				if (nodes_u.HasValue(node, nodeposition))
//				{
				
				/* get coordinates */
				face_coords.SetLocal(face_nodes);

				/* set up split integration */
				int normal_type = fSurfaceElementFacesType(i,j);
				
				if (fIndicator == "FCC_3D")
					t_surface = fSurfaceCB[normal_type]->SurfaceThickness();
				else if (fIndicator == "FCC_EAM")
					t_surface = fEAMSurfaceCB[normal_type]->SurfaceThickness();
				else
					ExceptionT::GeneralFail(caller, "unrecognized SCB \"%s\"", fIndicator.Pointer());
	
				fSplitInitCoords = fLocInitCoords;
				SurfaceLayer(fSplitInitCoords, j, t_surface);

				/* remove bulk contribution to surface layer (see TotalLagrangianT::FormKd) */
				const double* Det    = fSplitShapes->IPDets();
				const double* Weight = fSplitShapes->IPWeights();
				fSplitShapes->SetDerivatives(); /* set coordinate mapping over the split domain */
				fSplitShapes->TopIP();
				fShapes->TopIP(); /* synch bulk shape functions */
//				cout << "fNEEvec before = " << fNEEvec << endl;
				while (fSplitShapes->NextIP())
				{
					/* synch bulk shape functions */
					fShapes->NextIP();
				
					/* ip coordinates in the split domain */
					fSplitShapes->IPCoords(ip_coords_X);
					
					/* map ip coordinates to bulk parent domain */
					shape.ParentDomain().MapToParentDomain(fLocInitCoords, ip_coords_X, ip_coords_Xi);

					/* bulk shape functions/derivatives */
					shape.GradU(fLocInitCoords, DXi_DX, ip_coords_Xi, Na, DNa_Xi);
					DXi_DX.Inverse();
					shape.TransformDerivatives(DXi_DX, DNa_Xi, DNa_X);

					/* deformation gradient/shape functions/derivatives at the surface ip */
					dMatrixT& F = fF_List[fSplitShapes->CurrIP()];
					shape.ParentDomain().Jacobian(fLocDisp, DNa_X, F);
					F.PlusIdentity();

					/* F^(-1) */
					double J = F.Det();
					if (J <= 0.0)
						ExceptionT::BadJacobianDet(caller);
					else
						F_inv.Inverse(F);

					/* bulk material model */
					ContinuumMaterialT* pcont_mat = (*fMaterialList)[element_card.MaterialNumber()];
					fCurrMaterial = (SolidMaterialT*) pcont_mat;

					/* get Cauchy stress */
					(fCurrMaterial->s_ij()).ToMatrix(cauchy);

					/* compute PK1/J */
					PK1.MultABT(cauchy, F_inv);

					/* Wi,J PiJ */
					shape.GradNa(DNa_X, fGradNa);
					WP.MultAB(PK1, fGradNa);

					/* accumulate */
					fRHS.AddScaled(J*constKd*(*Weight++)*(*Det++), fNEEvec);
					
					/* If normal code = 1, check what happens */
					//force1.AddScaled(J*constKd*(*Weight++)*(*Det++), fNEEvec);
				}

//				cout << "fNEEvec = " << fNEEvec << endl;

				/* integrate over the face */
				int face_ip;
				fSurfaceCBSupport->SetCurrIP(face_ip);
				const double* w = surf_shape.Weight();				
				for (face_ip = 0; face_ip < nsi; face_ip++) {

					/* coordinate mapping on face */
					surf_shape.DomainJacobian(face_coords, face_ip, jacobian);
					double detj = surf_shape.SurfaceJacobian(jacobian);
				
					/* ip coordinates on face */
					surf_shape.Interpolate(face_coords, ip_coords_X, face_ip);
					
					/* ip coordinates in bulk parent domain */
					shape.ParentDomain().MapToParentDomain(fLocInitCoords, ip_coords_X, ip_coords_Xi);

					/* bulk shape functions/derivatives */
					shape.GradU(fLocInitCoords, DXi_DX, ip_coords_Xi, Na, DNa_Xi);
					DXi_DX.Inverse();
					shape.TransformDerivatives(DXi_DX, DNa_Xi, DNa_X);

					/* deformation gradient/shape functions/derivatives at the surface ip */
					dMatrixT& F = fF_Surf_List[face_ip];
					shape.ParentDomain().Jacobian(fLocDisp, DNa_X, F);
					F.PlusIdentity();
					
					/* F^-1 */
					double J = F.Det();
					if (J <= 0.0)
						ExceptionT::BadJacobianDet(caller);
					else
						F_inv.Inverse(F);
					
					/* stress at the surface */
					if (fIndicator == "FCC_3D")
						(fSurfaceCB[normal_type]->s_ij()).ToMatrix(cauchy);
					else if (fIndicator == "FCC_EAM")
					{
						(fEAMSurfaceCB[normal_type]->s_ij()).ToMatrix(cauchy);
					}
					else
						int blah = 0;
				
					/* compute PK1/J */
					PK1.MultABT(cauchy, F_inv);
					
					/* Wi,J PiJ */
					shape.GradNa(DNa_X, fGradNa);
					WP.MultAB(PK1, fGradNa);

					/* accumulate */
					fRHS.AddScaled(-J*constKd*w[face_ip]*detj, fNEEvec);
					
					/* If normal code = 1, check what happens */
					/* USING ADD SCALED SCREWS STUFF UP FOR SOME REASON */			
					//force2.AddScaled(-J*constKd*w[face_ip]*detj, fNEEvec);
					count ++;

				}
				
				/* MOVED HERE TEMPORARILIY */
				//totalforce -= force1;
				//totalforce += force2;
				
				/* fRHS different from totalforce - WHY? */
				//cout << "fRHS = " << fRHS << endl;
				//cout << "totalforce = " << totalforce << endl;
			
				/* assemble force */
				int dex = 0;
				dArrayT asdf(3);
				asdf = 0.0;
			
				if (nodes_u.HasValue(node, nodeposition))
				{
					for (int i = 0; i < nodes_u.Length(); i++)
					{
					//	cout << "nodes_u = " << nodes_u << endl;
						if (nodes_u[i] == node)
						{
							/* components for node */
							//nodalforce.Set(NumDOF(), fRHS.Pointer(dex));
							//force3.Set(NumDOF(), fRHS.Pointer(dex));
							//force3 += nodalforce;
							/* accumulate */
							//force += force3;
							asdf.CopyPart(0,fRHS,dex,3);
						//cout << "asdf = " << asdf << endl;
							reactionforce += asdf;
							force += asdf;
							//cout << "reactionforce = " << reactionforce << endl;
							//cout << "nodalforce = " << force << endl;
						}
						dex += NumDOF();
					}			
				
				} /* found node */
			} /* surface face */
	}	
//	cout << "nodalforce = " << force << endl;
//	cout << "reactionforce = " << reactionforce << endl;
}

/* information about subordinate parameter lists */
void TotalLagrangianCBSurfaceT::DefineSubs(SubListT& sub_list) const
{
	const char caller[] = "TotalLagrangianCBSurfaceT::TakeParameterList";
	
	/* inherited */
	TotalLagrangianT::DefineSubs(sub_list);

	/* list of passivated surfaces - side set list */
	sub_list.AddSub("passivated_surface_ID_list", ParameterListT::ZeroOrOnce);	
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
	int nen = NumElementNodes();

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

	/* START DISTINGUISHING HERE BETWEEN FCCEAMCB AND FCC3D */
	/* find parameter list for the bulk material */
	int num_blocks = list.NumLists("large_strain_element_block");
	if (num_blocks > 1)
		ExceptionT::GeneralFail(caller, "expecting only 1 not %d element blocks", num_blocks);
	const ParameterListT& block = list.GetList("large_strain_element_block");
	const ParameterListT& mat_list = block.GetListChoice(*this, "large_strain_material_choice");
	const ArrayT<ParameterListT>& mat = mat_list.Lists();
	const ParameterListT& bulk_params = mat[0];
	const char* blah = bulk_params.Name();
	fIndicator = blah;
	//if (bulk_params.Name() != "FCC_3D")
	//if (fIndicator != "FCC_3D")
	//	ExceptionT::GeneralFail(caller, "expecting \"FCC_3D or FCC_EAM\" not \"%s\"", bulk_params.Name().Pointer());
	
	/* Initialize either fSurfaceCB or fEAMSurfaceCB depending on fIndicator */
	if (fIndicator == "FCC_3D")
	{
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
	}
	else if (fIndicator == "FCC_EAM")
	{
		/* initialize surface information & create all possible (6) surface clusters */
		fNormal.Dimension(nfs);
		fEAMSurfaceCB.Dimension(nfs);
		fEAMSurfaceCB = NULL;

		/* get pointer to the bulk model */
		EAMFCC3DMatT* EAMfcc_3D = NULL;
		if (fMaterialList->Length() == 1) {
			ContinuumMaterialT* pcont_mat = (*fMaterialList)[0];
			EAMfcc_3D = TB_DYNAMIC_CAST(EAMFCC3DMatT*, pcont_mat);
			if (!EAMfcc_3D) ExceptionT::GeneralFail(caller, "could not resolve EAMFCC3D material");
		} else ExceptionT::GeneralFail(caller, "expecting 1 not %d materials", fMaterialList->Length());

		/* Update parameter list for EAMFCC3D_Surf to include the surface normal codes */
		for (int i = 0; i < nfs; i++)
		{
			/* face normal */
			fNormal[i].Dimension(nsd);
			fNormal[i] = normals(i);
	
			/* face C-B model */
			fEAMSurfaceCB[i] = new EAMFCC3DMatT_surf;
			fEAMSurfaceCB[i]->SetFSMatSupport(fSurfaceCBSupport);

			/* pass parameters to the surface model, including surface normal code */
			ParameterListT surf_params = bulk_params;
			surf_params.SetName("FCC_EAM_Surf");
			surf_params.AddParameter(i, "normal_code");
			
			/* ADD DUMMY VALUE FOR SHELLS SINCE SHELLS NOT READ IN IN INPUT FILE */
			surf_params.AddParameter(4, "shells");

			/* Initialize a different EAMFCC3D_Surf for each different surface normal type (6 total) */
			fEAMSurfaceCB[i]->TakeParameterList(surf_params);
		}
	}
	else ExceptionT::GeneralFail(caller, "expecting \"FCC_3D or FCC_EAM\" not \"%s\"", bulk_params.Name().Pointer());

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
	
	/* process passivated surfaces */
	const ParameterListT* passivated_surfaces = list.List("passivated_surface_ID_list");
	if (passivated_surfaces) {
		ArrayT<StringT> ss_ID;
		StringListT::Extract(*passivated_surfaces, ss_ID);
		if (ss_ID.Length() > 0) {

			/* need map into the surface element list */
			InverseMapT surf_elem_map;
			surf_elem_map.SetOutOfRange(InverseMapT::Throw);
			surf_elem_map.SetMap(fSurfaceElements);

			/* model manager */
			ModelManagerT& model = ElementSupport().ModelManager();

			/* loop over side set ID's */
			for (int i = 0; i < ss_ID.Length(); i++) {
				const StringT& id = ss_ID[i];
				
				/* side set parameters */
				iArray2DT sides = model.SideSet(id);
				const StringT& block_id = model.SideSetGroupID(id);

				if (sides.MajorDim() > 0) {

					/* convert to element numbering within the group */
					iArrayT elems(sides.MajorDim());
					sides.ColumnCopy(0, elems);
					BlockToGroupElementNumbers(elems, block_id);
					sides.SetColumn(0, elems);
			
					/* mark passivated faces */
					for (int j = 0; j < sides.MajorDim(); j++) {
						int elem = sides(j,0);
						int s_elem = surf_elem_map.Map(elem);
						int side = sides(j,1);
						
						/* mark as non-surface (by giving negative neighbor id) */
						fSurfaceElementNeighbors(s_elem,side) = -2;
					}
				}
			}
		}
	}

//TEMP
if (0) {
fSurfaceElements++;
cout << "surface elements:\n" << fSurfaceElements.wrap(10) << "\n";
fSurfaceElementNeighbors.WriteNumbered(cout);
fSurfaceElements--;
cout << endl;
ExceptionT::Stop();
}

	/* initialize data for split integration */
	fSplitInitCoords.Dimension(nen, nsd);
	fSplitShapes = new ShapeFunctionT(GeometryCode(), NumIP(), fSplitInitCoords);
	fSplitShapes->Initialize();
}

/*************************************************************************
 * Protected
 *************************************************************************/

/* form group contribution to the stiffness matrix */
void TotalLagrangianCBSurfaceT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
	const char caller[] = "TotalLagrangianCBSurfaceT::LHSDriver";
	
	/* inherited - bulk contribution */
	TotalLagrangianT::LHSDriver(sys_type);
	
	/* time integration parameters */
	double constK = 0.0;
	int formK = fIntegrator->FormK(constK);
	if (!formK) return;
	
	/* dimensions */
	const ShapeFunctionT& shape = ShapeFunction();
	int nsd = shape.NumSD();                          // # of spatial dimensions in problem
	int nfs = shape.NumFacets();                      // # of total possible element faces
	int nsi = shape.FacetShapeFunction(0).NumIP();    // # IPs per surface face
	int nfn = shape.FacetShapeFunction(0).NumNodes(); // # nodes on each surface face
	int nen = NumElementNodes();                      // # nodes in bulk element
	
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;
	
	/* loop over surface elements */
	dMatrixT jacobian(nsd, nsd-1);
	iArrayT face_nodes(nfn), face_nodes_index(nfn);
	LocalArrayT face_coords(LocalArrayT::kInitCoords, nfn, nsd);
	ElementSupport().RegisterCoordinates(face_coords);
	dArrayT ip_coords_X(nsd);
	dArrayT ip_coords_Xi(nsd);
	dArrayT Na(nen);
	dArray2DT DNa_X(nsd,nen), DNa_Xi(nsd,nen), DNa_x(nsd,nen);
	dMatrixT DXi_DX(nsd);
	dMatrixT F_inv(nsd);
	dMatrixT PK1(nsd), cauchy(nsd);
	double t_surface;
	for (int i = 0; i < fSurfaceElements.Length(); i++)
	{
		/* bulk element information */
		int element = fSurfaceElements[i];
		const ElementCardT& element_card = ElementCard(element);
		fLocInitCoords.SetLocal(element_card.NodesX()); /* reference coordinates over bulk element */
		fLocDisp.SetLocal(element_card.NodesU()); /* displacements over bulk element */
	
		/* initialize */
		fStressStiff = 0.0;
		fLHS = 0.0;

		/* integrate surface contribution to nodal forces */
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

				/* set up split integration */
				int normal_type = fSurfaceElementFacesType(i,j);
				//double t_surface = fSurfaceCB[normal_type]->SurfaceThickness();
				
				if (fIndicator == "FCC_3D")
					t_surface = fSurfaceCB[normal_type]->SurfaceThickness();
				else if (fIndicator == "FCC_EAM")
					t_surface = fEAMSurfaceCB[normal_type]->SurfaceThickness();
				else
					ExceptionT::GeneralFail(caller, "unrecognized SCB \"%s\"", fIndicator.Pointer());
				
				fSplitInitCoords = fLocInitCoords;
				SurfaceLayer(fSplitInitCoords, j, t_surface);

				/* remove bulk contribution to surface layer (see TotalLagrangianT::FormKd) */
				const double* Det    = fSplitShapes->IPDets();
				const double* Weight = fSplitShapes->IPWeights();
				fSplitShapes->SetDerivatives(); /* set coordinate mapping over the split domain */
				fSplitShapes->TopIP();
				fShapes->TopIP(); /* synch bulk shape functions */				
				while (fSplitShapes->NextIP())
				{
					/* synch bulk shape functions */
					fShapes->NextIP();

				/* MAPPING/DEFORMATION */

					/* ip coordinates in the split domain */
					fSplitShapes->IPCoords(ip_coords_X);
					
					/* map ip coordinates to bulk parent domain */
					shape.ParentDomain().MapToParentDomain(fLocInitCoords, ip_coords_X, ip_coords_Xi);

					/* bulk shape functions/derivatives */
					shape.GradU(fLocInitCoords, DXi_DX, ip_coords_Xi, Na, DNa_Xi);
					DXi_DX.Inverse();
					shape.TransformDerivatives(DXi_DX, DNa_Xi, DNa_X);

					/* deformation gradient/shape functions/derivatives at the surface ip */
					dMatrixT& F = fF_List[fSplitShapes->CurrIP()];
					shape.ParentDomain().Jacobian(fLocDisp, DNa_X, F);
					F.PlusIdentity();

					/* F^(-1) */
					double J = F.Det();
					if (J <= 0.0)
						ExceptionT::BadJacobianDet(caller);
					else
						F_inv.Inverse(F);

					/* bulk material model */
					ContinuumMaterialT* pcont_mat = (*fMaterialList)[element_card.MaterialNumber()];
					fCurrMaterial = (SolidMaterialT*) pcont_mat;

				/* STRESS STIFFNESS */

					/* shape function gradient wrt current configuration */
					shape.TransformDerivatives(F_inv, DNa_X, DNa_x);

					/* get Cauchy stress */
					(fCurrMaterial->s_ij()).ToMatrix(cauchy);

					/* integration weight */
					double scale = -constK*(*Det++)*(*Weight++)*J;

					/* integration constants */
					cauchy *= scale;

					/* using the stress symmetry - watch big X vs. little x */
					shape.GradNa(DNa_x, fGradNa);
					fStressStiff.MultQTBQ(fGradNa, cauchy, format, dMatrixT::kAccumulate);

				/* MATERIAL STIFFNESS */

					/* strain displacement matrix */
					Set_B(DNa_x, fB);

					/* Get D Matrix */
					fD.SetToScaled(scale, fCurrMaterial->c_ijkl());

					/* accumulate */
					fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);
				}

				/* integrate over the face */
				int face_ip;
				fSurfaceCBSupport->SetCurrIP(face_ip);
				const double* w = surf_shape.Weight();
				for (face_ip = 0; face_ip < nsi; face_ip++)
				{
				/* MAPPING/DEFORMATION */
				
					/* reference coordinate mapping on face */
					surf_shape.DomainJacobian(face_coords, face_ip, jacobian);
					double detj = surf_shape.SurfaceJacobian(jacobian);
				
					/* ip coordinates on face */
					surf_shape.Interpolate(face_coords, ip_coords_X, face_ip);
					
					/* ip coordinates in bulk parent domain */
					shape.ParentDomain().MapToParentDomain(fLocInitCoords, ip_coords_X, ip_coords_Xi);

					/* shape functions/derivatives */
					shape.GradU(fLocInitCoords, DXi_DX, ip_coords_Xi, Na, DNa_Xi);
					DXi_DX.Inverse();
					shape.TransformDerivatives(DXi_DX, DNa_Xi, DNa_X);

					/* deformation gradient/shape functions/derivatives at the surface ip */
					dMatrixT& F = fF_Surf_List[face_ip];
					shape.ParentDomain().Jacobian(fLocDisp, DNa_X, F);
					F.PlusIdentity();
					
					/* F^-1 */
					double J = F.Det();
					F_inv.Inverse(F);

				/* STRESS STIFFNESS */
					
					/* shape function gradient wrt current configuration */
					shape.TransformDerivatives(F_inv, DNa_X, DNa_x);
					
					/* stress at the surface - USE INDICATOR HERE */
					if (fIndicator == "FCC_3D")
						(fSurfaceCB[normal_type]->s_ij()).ToMatrix(cauchy);
					else if (fIndicator == "FCC_EAM")
						(fEAMSurfaceCB[normal_type]->s_ij()).ToMatrix(cauchy);
					else
						int blah = 0;
					
					/* integration weight */
					double scale = constK*detj*w[face_ip]*J;
					
					/* integration constants */
					cauchy *= scale;
					
					/* using the stress symmetry - watch big X vs. little x */
					shape.GradNa(DNa_x, fGradNa);
					fStressStiff.MultQTBQ(fGradNa, cauchy, format, dMatrixT::kAccumulate);

				/* MATERIAL STIFFNESS */
				
					/* strain displacement matrix */
					Set_B(DNa_x, fB);
					
					/* Get D Matrix */
					if (fIndicator == "FCC_3D")
						fD.SetToScaled(scale, fSurfaceCB[normal_type]->c_ijkl());
					else if (fIndicator == "FCC_EAM")
						fD.SetToScaled(scale, fEAMSurfaceCB[normal_type]->c_ijkl());
					else
						int blah = 0;
					
					/* accumulate */
					fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);
				}
			}
			
		/* add/expand stress stiffness contribution into fLHS */
		fLHS.Expand(fStressStiff, NumDOF(), dMatrixT::kAccumulate);
		
		/* assemble stiffness */
		ElementSupport().AssembleLHS(Group(), fLHS, element_card.Equations());		
	}
}

/* form group contribution to the residual */
void TotalLagrangianCBSurfaceT::RHSDriver(void)
{
	const char caller[] = "TotalLagrangianCBSurfaceT::RHSDriver";
	
	/* check wave speed directions */
	dArrayT normal(3), speeds(3);
	normal[0] = 0.0;
	normal[1] = 0.0;
	normal[2] = 1.0;

	/* inherited - bulk contribution */
	TotalLagrangianT::RHSDriver();

	/* time integration parameters */
	double constKd = 0.0;
	int formKd = fIntegrator->FormKd(constKd);
	if (!formKd) return;

	/* dimensions */
	const ShapeFunctionT& shape = ShapeFunction();
	int nsd = shape.NumSD();                          // # of spatial dimensions in problem
	int nfs = shape.NumFacets();                      // # of total possible element faces
	int nsi = shape.FacetShapeFunction(0).NumIP();    // # IPs per surface face
	int nfn = shape.FacetShapeFunction(0).NumNodes(); // # nodes on each surface face
	int nen = NumElementNodes();                      // # nodes in bulk element

	/* matrix alias to fNEEvec */
	dMatrixT WP(nsd, fStressStiff.Rows(), fNEEvec.Pointer());

	/* loop over surface elements */
	dMatrixT jacobian(nsd, nsd-1);
	LocalArrayT face_coords(LocalArrayT::kInitCoords, nfn, nsd);
	iArrayT face_nodes(nfn), face_nodes_index(nfn);
	ElementSupport().RegisterCoordinates(face_coords);
	dArrayT ip_coords_X(nsd);
	dArrayT ip_coords_Xi(nsd);
	dArrayT Na(nen);
	dArray2DT DNa_X(nsd,nen), DNa_Xi(nsd,nen);
	dMatrixT DXi_DX(nsd);
	dMatrixT F_inv(nsd);
	dMatrixT PK1(nsd), cauchy(nsd);
	double t_surface;
	
	/* TEMPORARY */
	dArrayT force1,force2,reactionforce,totalforce;
	force1.Dimension(24);
	force2.Dimension(24);
	reactionforce.Dimension(3);
	totalforce.Dimension(24);
	force2 = 0.0;
	force1 = 0.0;
	reactionforce = 0.0;
	totalforce = 0.0;
	
	int count = 0;
	
	for (int i = 0; i < fSurfaceElements.Length(); i++)
	{
		/* bulk element information */
		int element = fSurfaceElements[i];
		const ElementCardT& element_card = ElementCard(element);
		fLocInitCoords.SetLocal(element_card.NodesX()); /* reference coordinates over bulk element */
		fLocDisp.SetLocal(element_card.NodesU()); /* displacements over bulk element */
	
		/* integrate surface contribution to nodal forces */
		fRHS = 0.0;
		for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) /* loop over faces */
			if (fSurfaceElementNeighbors(i,j) == -1) /* no neighbor => surface */
			{
				/* face parent domain */
				const ParentDomainT& surf_shape = shape.FacetShapeFunction(j);
				
				/* TEMPORARY BY HSP TO TEST REACTION FORCE OUTPUT */
				const iArrayT& nodes_u = element_card.NodesU();
			
				/* collect coordinates of face nodes */
				ElementCardT& element_card = ElementCard(fSurfaceElements[i]);
				shape.NodesOnFacet(j, face_nodes_index);	// fni = 4 nodes of surface face
				face_nodes.Collect(face_nodes_index, element_card.NodesX());
				face_coords.SetLocal(face_nodes);

				/* set up split integration */
				int normal_type = fSurfaceElementFacesType(i,j);
				//cout << "normal_type = " << normal_type << endl;
				
				if (fIndicator == "FCC_3D")
					t_surface = fSurfaceCB[normal_type]->SurfaceThickness();
				else if (fIndicator == "FCC_EAM")
					t_surface = fEAMSurfaceCB[normal_type]->SurfaceThickness();
				else
					int blah = 0;
	
				fSplitInitCoords = fLocInitCoords;
				SurfaceLayer(fSplitInitCoords, j, t_surface);

				/* remove bulk contribution to surface layer (see TotalLagrangianT::FormKd) */
				const double* Det    = fSplitShapes->IPDets();
				const double* Weight = fSplitShapes->IPWeights();
				fSplitShapes->SetDerivatives(); /* set coordinate mapping over the split domain */
				fSplitShapes->TopIP();
				fShapes->TopIP(); /* synch bulk shape functions */
				while (fSplitShapes->NextIP())
				{
					/* synch bulk shape functions */
					fShapes->NextIP();
				
					/* ip coordinates in the split domain */
					fSplitShapes->IPCoords(ip_coords_X);
					
					/* map ip coordinates to bulk parent domain */
					shape.ParentDomain().MapToParentDomain(fLocInitCoords, ip_coords_X, ip_coords_Xi);

					/* bulk shape functions/derivatives */
					shape.GradU(fLocInitCoords, DXi_DX, ip_coords_Xi, Na, DNa_Xi);
					DXi_DX.Inverse();
					shape.TransformDerivatives(DXi_DX, DNa_Xi, DNa_X);

					/* deformation gradient/shape functions/derivatives at the surface ip */
					dMatrixT& F = fF_List[fSplitShapes->CurrIP()];
					shape.ParentDomain().Jacobian(fLocDisp, DNa_X, F);
					F.PlusIdentity();

					/* F^(-1) */
					double J = F.Det();
					if (J <= 0.0)
						ExceptionT::BadJacobianDet(caller);
					else
						F_inv.Inverse(F);

					/* bulk material model */
					ContinuumMaterialT* pcont_mat = (*fMaterialList)[element_card.MaterialNumber()];
					fCurrMaterial = (SolidMaterialT*) pcont_mat;

					/* get Cauchy stress */
					(fCurrMaterial->s_ij()).ToMatrix(cauchy);

					/* compute PK1/J */
					PK1.MultABT(cauchy, F_inv);

					/* Wi,J PiJ */
					shape.GradNa(DNa_X, fGradNa);
					WP.MultAB(PK1, fGradNa);

					/* accumulate */
					fRHS.AddScaled(J*constKd*(*Weight++)*(*Det++), fNEEvec);
					
					/* If normal code = 1, check what happens */

					//force1.AddScaled(J*constKd*(*Weight++)*(*Det++), fNEEvec);
				}

				/* integrate over the face */
				int face_ip;
				fSurfaceCBSupport->SetCurrIP(face_ip);
				const double* w = surf_shape.Weight();				
				for (face_ip = 0; face_ip < nsi; face_ip++) {

					/* coordinate mapping on face */
					surf_shape.DomainJacobian(face_coords, face_ip, jacobian);
					double detj = surf_shape.SurfaceJacobian(jacobian);
				
					/* ip coordinates on face */
					surf_shape.Interpolate(face_coords, ip_coords_X, face_ip);
					
					/* ip coordinates in bulk parent domain */
					shape.ParentDomain().MapToParentDomain(fLocInitCoords, ip_coords_X, ip_coords_Xi);

					/* bulk shape functions/derivatives */
					shape.GradU(fLocInitCoords, DXi_DX, ip_coords_Xi, Na, DNa_Xi);
					DXi_DX.Inverse();
					shape.TransformDerivatives(DXi_DX, DNa_Xi, DNa_X);

					/* deformation gradient/shape functions/derivatives at the surface ip */
					dMatrixT& F = fF_Surf_List[face_ip];
					shape.ParentDomain().Jacobian(fLocDisp, DNa_X, F);
					F.PlusIdentity();
					
					/* F^-1 */
					double J = F.Det();
					if (J <= 0.0)
						ExceptionT::BadJacobianDet(caller);
					else
						F_inv.Inverse(F);
					
					/* stress at the surface */
					if (fIndicator == "FCC_3D")
						(fSurfaceCB[normal_type]->s_ij()).ToMatrix(cauchy);
					else if (fIndicator == "FCC_EAM")
					{
						(fEAMSurfaceCB[normal_type]->s_ij()).ToMatrix(cauchy);
						//fEAMSurfaceCB[normal_type]->WaveSpeeds(normal, speeds);
						//cout << "wave speeds are " << speeds << endl;
					}
					else
						int blah = 0;
				
					/* compute PK1/J */
					PK1.MultABT(cauchy, F_inv);
					
					/* Wi,J PiJ */
					shape.GradNa(DNa_X, fGradNa);
					WP.MultAB(PK1, fGradNa);

					/* accumulate */
					fRHS.AddScaled(-J*constKd*w[face_ip]*detj, fNEEvec);
					
					/* If normal code = 1, check what happens */

					//force2.AddScaled(J*constKd*w[face_ip]*detj, fNEEvec);
					count++;

				}	
//				cout << "force2 = " << force2 << endl;
				/* assemble force */
				//totalforce += force1;
				//totalforce += force2;
				
				int dex = 0;
				dArrayT asdf(3);
				asdf = 0.0;
				for (int i = 0; i < nodes_u.Length(); i++)
				{
					/* JUST TRY A SINGLE NODE */
					if (nodes_u[i] == 8)
					{
						//cout << "nodes_u = " << nodes_u << endl;
						/* components for node */
						//reactionforce.Set(NumDOF(), fRHS.Pointer(dex));
						//cout << "dex = " << dex << endl;
						asdf.CopyPart(0,fRHS,dex,3);
						//cout << "asdf = " << asdf << endl;
						reactionforce += asdf;
						
						/* accumulate */
						//force += nodalforce;
					}
					dex += NumDOF();
				}	
//			cout << "reactionforce_RHDS = " << reactionforce << endl;
			}
			//cout << "reactionforce_RHDS= " << reactionforce << endl;
		/* assemble forces */
		ElementSupport().AssembleRHS(Group(), fRHS, element_card.Equations());
//	cout << "reactionforce_RHDS = " << reactionforce << endl;
	}
}
	

/***********************************************************************
 * Private
 ***********************************************************************/

/* reduce the coordinates to a surface layer on the given face */
void TotalLagrangianCBSurfaceT::SurfaceLayer(LocalArrayT& coords, int face, double thickness) const
{
	const char caller[] = "TotalLagrangianCBSurfaceT::SurfaceLayer";
	if (GeometryCode() != GeometryT::kHexahedron)
		ExceptionT::GeneralFail(caller, "implemented for hex geometry only");
	int nen = NumElementNodes();
	if (nen != 8 && nen != 20)
		ExceptionT::GeneralFail(caller, "implemented for 8 or 20 node hexes only: %d", nen);

	/* transpose coordinate data */
	double d60[3*20]; /* oversize */
	dArray2DT coords_tmp(nen, 3, d60);
	coords.ReturnTranspose(coords_tmp);

	int i8[8]; /* oversize */
	iArrayT face_nodes(nen, i8);
	ShapeFunction().NodesOnFacet(face, face_nodes);
	double v1[3], v2[3];
	Vector(coords_tmp(face_nodes[1]), coords_tmp(face_nodes[0]), v1);
	Vector(coords_tmp(face_nodes[1]), coords_tmp(face_nodes[2]), v2);
	double normal[3];
	CrossProduct(v2, v1, normal);
	Scale(normal, 1.0/sqrt(Dot(normal,normal)));

	/* opposite face information - vertex nodes only */
	int opp_face[6] = {1,0,4,5,2,3};
	int opp_face_nodes_dat[6*4] = {
		4,7,6,5,
		0,1,2,3,
		3,2,6,7,
		0,3,7,4,
		1,0,4,5,
		2,1,5,6};
	iArray2DT opp_face_nodes(6, 4, opp_face_nodes_dat);

	/* "shorten" vectors to back face */
	for (int i = 0; i < 4; i++)
	{
		/* nodes */
		int n_front = face_nodes[i];
		int n_back  = opp_face_nodes(face,i); 
	
		/* vector to back face */
		Vector(coords_tmp(n_front), coords_tmp(n_back), v1);
		double h = -Dot(v1,normal);
		if (h < 0.0) ExceptionT::GeneralFail(caller, "geometry error");
		if (h < 2.0*thickness)
			ExceptionT::GeneralFail(caller, "layer thickness %g exceeds half element depth %g",
				thickness, h/2.0);
		
		/* compute point */
		Scale(v1, thickness/h);
		Sum(coords_tmp(n_front), v1, v2);
		
		/* store */
		coords_tmp.SetRow(n_back, v2);
	}

	/* move mid-side nodes */
	if (nen == 20) {
		int edges_dat[4*2] = {0,1,1,2,2,3,3,0};
		iArray2DT edges(4, 2, edges_dat);
		for (int i = 0; i < 6; i++)
			if (i != face) {
				ShapeFunction().NodesOnFacet(i, face_nodes);
				for (int j = 0; j < 4; j++) /* loop over edges */ {
					Sum(coords_tmp(face_nodes[edges(j,0)]), 
					    coords_tmp(face_nodes[edges(j,1)]), 
					    coords_tmp(face_nodes[j+4]));
					coords_tmp.ScaleRow(face_nodes[j+4], 0.5);
				}
			}
	}
	
	/* write back */
	coords.FromTranspose(coords_tmp);	
}
