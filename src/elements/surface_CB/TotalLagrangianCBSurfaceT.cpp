/* $Id: TotalLagrangianCBSurfaceT.cpp,v 1.21 2005-07-21 16:21:23 paklein Exp $ */
#include "TotalLagrangianCBSurfaceT.h"

#include "ModelManagerT.h"
#include "ShapeFunctionT.h"
#include "FCC3D_Surf.h"
#include "FCC3D.h"
#include "MaterialListT.h"
#include "eIntegratorT.h"

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
				double t_surface = fSurfaceCB[normal_type]->SurfaceThickness();
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
					
					/* stress at the surface */
					(fSurfaceCB[normal_type]->s_ij()).ToMatrix(cauchy);

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
					fD.SetToScaled(scale, fSurfaceCB[normal_type]->c_ijkl());
					
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
				face_coords.SetLocal(face_nodes);

				/* set up split integration */
				int normal_type = fSurfaceElementFacesType(i,j);
				double t_surface = fSurfaceCB[normal_type]->SurfaceThickness();
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
					(fSurfaceCB[normal_type]->s_ij()).ToMatrix(cauchy);

					/* compute PK1/J */
					PK1.MultABT(cauchy, F_inv);
					
					/* Wi,J PiJ */
					shape.GradNa(DNa_X, fGradNa);
					WP.MultAB(PK1, fGradNa);

					/* accumulate */
					fRHS.AddScaled(-J*constKd*w[face_ip]*detj, fNEEvec);
				}				
			}
			
		/* assemble forces */
		ElementSupport().AssembleRHS(Group(), fRHS, element_card.Equations());
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
		ExceptionT::GeneralFail(caller, "implemented for 4 or 20 node hexes only: %d", nen);

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
