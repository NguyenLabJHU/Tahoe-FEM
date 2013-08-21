#include "FSDielectricElastomerQ1P0SurfaceT.h"
#include "FSDEMatQ1P0T.h"
#include "FSDEMatSupportQ1P0T.h"
#include "ShapeFunctionT.h"

/* From TLCBSurfaceT */
#include "ModelManagerT.h"
#include "MaterialListT.h"
#include "eIntegratorT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "InverseMapT.h"
#include "RowAutoFill2DT.h"
#include "OutputSetT.h"

namespace Tahoe {

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

  FSDielectricElastomerQ1P0SurfaceT::FSDielectricElastomerQ1P0SurfaceT(
      const ElementSupportT& support) :
    FSDielectricElastomerQ1P0T(support),
    fLocCurrCoords(LocalArrayT::kCurrCoords),
    fSurfTension(0.0),
    fSurfaceCBSupport(NULL)
  {
    SetName("dielectric_elastomer_Q1P0Elastocapillary");
  }


/* Destructor */
  FSDielectricElastomerQ1P0SurfaceT::~FSDielectricElastomerQ1P0SurfaceT()
  {
	  /* TLCBSurface Stuff */
	  /* free surface models */
	  delete fSurfaceCBSupport;
  }

  // specify parameters needed by the interface  
  void FSDielectricElastomerQ1P0SurfaceT::DefineParameters(ParameterListT& list) const
  {
    // inherited
    FSDielectricElastomerQ1P0T::DefineParameters(list);
    
	/* TLCBSurface Stuff */
	/* associated fields */
	ParameterT output_surface(ParameterT::Boolean, "output_surface");
	output_surface.SetDefault(false);
	list.AddParameter(output_surface);	 
  }

/* From TLCBSurfaceT */
/* information about subordinate parameter lists */
void FSDielectricElastomerQ1P0SurfaceT::DefineSubs(SubListT& sub_list) const
{
	const char caller[] = "FSDielectricElastomerQ1P0SurfaceT::DefineSubs";
	
	/* inherited */
	FSDielectricElastomerQ1P0T::DefineSubs(sub_list);

	/* list of passivated surfaces - side set list */
	sub_list.AddSub("passivated_surface_ID_list", ParameterListT::ZeroOrOnce);	
}

  // accept parameter list
  void FSDielectricElastomerQ1P0SurfaceT::TakeParameterList(const ParameterListT& list)
  {
  	const char caller[] = "FSDielectricElastomerQ1P0SurfaceT::TakeParameterList";
  	
    // inherited
    FSDielectricElastomerQ1P0T::TakeParameterList(list);

    /* Define matrix sizes */
    int nen = NumElementNodes();
    int nsd = NumSD();
    int nel = nen;	// # electrical DOFs per element
    int nme = nen * nsd;	// # of mechanical DOFs per element
    int dof = nsd + 1;	// total # of DOFs per node (mech + elec)
    int neq = nen * dof;	// total # of DOFs per element (mech + elec)

	/* Dimension matrices */
    fAmm_mat2.Dimension(nme, nme);
    fAmm_geo2.Dimension(nen, nen);

	/* allocate workspace - from UpdatedLagrangianT.cpp */
	fGradNa.Dimension(nsd, nen);
	
	/* TLCBSurface Stuff */
	/* the shape functions */
	const ShapeFunctionT& shape = ShapeFunction();
	int nfs = shape.NumFacets();	// # of total possible surface facets?
	int nsi = shape.FacetShapeFunction(0).NumIP();		// # IPs per surface face (2x2=4 for 2D surface)
	int nfn = shape.FacetShapeFunction(0).NumNodes();	// # nodes on each surface face?

	/* Dimension matrices */
	fB2.Dimension(dSymMatrixT::NumValues(NumSD()), NumSD()*NumElementNodes());
	fD2.Dimension(dSymMatrixT::NumValues(NumSD()));

	/* Define LHS type based upon analysis type, i.e. static vs. dynamic */
// 	int order = fIntegrator->Order();
// 	if (order == 2)
// 		fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
// 	else
// 		fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
// 		
// 	fLHS.Dimension(neq);
// 	fRHS.Dimension(neq);
	
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

	fNormal.Dimension(nfs);
	for (int i = 0; i < nfs; i++)
	{
		/* face normal */
		fNormal[i].Dimension(nsd);
		fNormal[i] = normals(i);
	}

	/* find parameter list for the bulk material */
	int num_blocks = list.NumLists("large_strain_element_block");
	if (num_blocks > 1)
		ExceptionT::GeneralFail(caller, "expecting only 1 not %d element blocks", num_blocks);
	const ParameterListT& block = list.GetList("large_strain_element_block");
	const ParameterListT& mat_list = block.GetListChoice(*this, "large_strain_material_choice");
	const ArrayT<ParameterListT>& mat = mat_list.Lists();
	const ParameterListT& bulk_params = mat[0];

	/* ACCESS LAMBDA USING BULK_PARAMS? */
	fSurfTension = bulk_params.GetParameter("gamma");

	/* output surface parameters */
	bool output_surface = list.GetParameter("output_surface");
//	if (fIndicator != "Tersoff_CB") { /* no surface output for the other surface types */
//		output_surface = false;
//	}
	
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
	
	/* nodes per surface type */
	RowAutoFill2DT<int> nodes_on_surfaces(nfs, 25);

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
				
				/* collect face nodes */
				if (output_surface) nodes_on_surfaces.AppendUnique(normal_type, face_nodes);
			}
	}
	
	/* copy nodes on face types data */
	if (output_surface) {
		fSurfaceNodes.Dimension(nodes_on_surfaces.MajorDim());
		for (int i = 0; i < fSurfaceNodes.Length(); i++) {
			fSurfaceNodes[i].Dimension(nodes_on_surfaces.MinorDim(i));
			fSurfaceNodes[i].Copy(nodes_on_surfaces(i));
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

  }

/***********************************************************************
 * Private
 ***********************************************************************/

/* calculate the LHS of residual, or element stiffness matrix */
  void FSDielectricElastomerQ1P0SurfaceT::FormStiffness(double constK)
  {
  	/* inherited - bulk contribution */
  	FSDielectricElastomerQ1P0T::FormStiffness(constK);
  	
  	/* matrix format */
    dMatrixT::SymmetryFlagT format = (fLHS.Format()
        == ElementMatrixT::kNonSymmetric)
        ? dMatrixT::kWhole
        : dMatrixT::kUpperOnly;
  	
	/* dimensions */
	const ShapeFunctionT& shape = ShapeFunction();
	int nsd = shape.NumSD();                          // # of spatial dimensions in problem
	int nfs = shape.NumFacets();                      // # of total possible element faces
	int nsi = shape.FacetShapeFunction(0).NumIP();    // # IPs per surface face
	int nfn = shape.FacetShapeFunction(0).NumNodes(); // # nodes on each surface face
	int nen = NumElementNodes();                      // # nodes in bulk element
	
	/* loop over surface elements */
	dMatrixT Q(nsd);
	dMatrixT jacobian(nsd, nsd-1);
	dMatrixT surfGradNa(nsd, nen);
	dArray2DT surfGradNa1(nsd, nen);
	iArrayT face_nodes(nfn), face_nodes_index(nfn);
	LocalArrayT face_coords(LocalArrayT::kInitCoords, nfn, nsd);
	LocalArrayT face_curr_coords(LocalArrayT::kCurrCoords, nfn, nsd);
	ElementSupport().RegisterCoordinates(face_coords);
	ElementSupport().RegisterCoordinates(face_curr_coords);
	dArrayT ip_coords_X(nsd);
	dArrayT ip_coords_Xi(nsd);
	dArrayT Na(nen);
	dArray2DT DNa_X(nsd,nen), DNa_Xi(nsd,nen), DNa_x(nsd,nen);
	dMatrixT DXi_DX(nsd);
	dMatrixT F_inv(nsd), cauchy2(nsd);
	
	for (int i = 0; i < fSurfaceElements.Length(); i++)
	{
		/* bulk element information */
		int element = fSurfaceElements[i];
		const ElementCardT& element_card = ElementCard(element);
		fLocInitCoords.SetLocal(element_card.NodesX()); /* reference coordinates over bulk element */
		fLocDisp.SetLocal(element_card.NodesU()); /* displacements over bulk element */

		/* initialize */
    	fAmm_mat2 = 0.0;
    	fAmm_geo2 = 0.0;

		/* integrate surface contribution to nodal forces */
		for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) /* loop over faces */
			if (fSurfaceElementNeighbors(i,j) == -1) /* no neighbor => surface */
//			if (fSurfaceElementNeighbors(i,j) == -1 && (fSurfaceElementFacesType(i,j) == 2 || fSurfaceElementFacesType(i,j) == 3)) /* no neighbor => surface */
			{			
				/* face parent domain */
				const ParentDomainT& surf_shape = shape.FacetShapeFunction(j);
			
				/* collect coordinates of face nodes */
				ElementCardT& element_card = ElementCard(fSurfaceElements[i]);
				shape.NodesOnFacet(j, face_nodes_index);	// fni = 4 nodes of surface face
				face_nodes.Collect(face_nodes_index, element_card.NodesX());
				face_coords.SetLocal(face_nodes);
				face_curr_coords.SetLocal(face_nodes);

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
				
					/* face normal (for all face IPs) */
					surf_shape.DomainJacobian(face_curr_coords, face_ip, jacobian);
					surf_shape.SurfaceJacobian(jacobian, Q);	// Last column of Q is normal vector to surface
				
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
					cauchy2 = Surf_Tension_Stress();
					
					/* integration weight */
					double scale = constK*detj*w[face_ip]*J;
					
					/* integration constants */
					cauchy2 *= scale;
					
					/* using the stress symmetry - watch big X vs. little x */
					shape.GradNa(DNa_x, fGradNa);

					/* MODIFY SHAPE FUNCTIONS USING SURFACE GRADIENT */
					surfGradNa = Surf_BMatrix(fGradNa, Q);
					surfGradNa1 = MatrixTodArray2DT(surfGradNa);
					
					/* using the stress symmetry */
//					fAmm_geo2.MultQTBQ(fGradNa, cauchy2, format, dMatrixT::kAccumulate);					

				/* MATERIAL STIFFNESS */
				
					/* strain displacement matrix */
//					Set_B(DNa_x, fB2);
					Set_B(surfGradNa1, fB2);
					
					/* Get D Matrix */
//					fD2.SetToScaled(scale, fSurfaceCB->c_ijkl());
					
					/* accumulate */
//					fAmm_mat2.MultQTBQ(fB2, fD2, format, dMatrixT::kAccumulate);
				}
			}
		/* stress stiffness into fLHS (i.e. fAmm_mat) */
//		fAmm_mat2.Expand(fAmm_geo2, NumDOF(), dMatrixT::kAccumulate);			

		/* Assemble into fLHS, or element stiffness matrix */
//		fLHS.AddBlock(0, 0, fAmm_mat2);
		
		/* assemble stiffness */
//		ElementSupport().AssembleLHS(Group(), fLHS, element_card.Equations());			
	}		
  }

/* Compute RHS, or residual of element equations */
  void FSDielectricElastomerQ1P0SurfaceT::FormKd(double constK)
  {  	
	const char caller[] = "FSDielectricElastomerQ1P0SurfaceT::FormKd";
	
	/* inherited - bulk contribution */
	FSDielectricElastomerQ1P0T::FormKd(constK);  
			  
	/* TLCBSurfaceT.cpp STUFF */
	/* dimensions */
	const ShapeFunctionT& shape = ShapeFunction();
	int nsd = shape.NumSD();                          // # of spatial dimensions in problem
	int nfs = shape.NumFacets();                      // # of total possible element faces
	int nsi = shape.FacetShapeFunction(0).NumIP();    // # IPs per surface face
	int nfn = shape.FacetShapeFunction(0).NumNodes(); // # nodes on each surface face
	int nen = NumElementNodes();                      // # nodes in bulk element
    
    /* Define mechanical and electrical residuals */
	dArrayT Rtotal2((nsd+1)*nen);
	Rtotal2 = 0.0;
	dArrayT Rmech2(nen*nsd);
	
	/* loop over surface elements */
	dMatrixT Q(nsd);	
	dMatrixT jacobian(nsd, nsd-1);
	iArrayT face_nodes(nfn), face_nodes_index(nfn);
	LocalArrayT face_coords(LocalArrayT::kInitCoords, nfn, nsd);
	LocalArrayT face_curr_coords(LocalArrayT::kCurrCoords, nfn, nsd);
	ElementSupport().RegisterCoordinates(face_coords);
	ElementSupport().RegisterCoordinates(face_curr_coords);
	dArrayT ip_coords_X(nsd);
	dArrayT ip_coords_Xi(nsd);
	dArrayT Na(nen);
	dArray2DT DNa_X(nsd,nen), DNa_Xi(nsd,nen), DNa_x(nsd,nen);
	dMatrixT DXi_DX(nsd);
	dMatrixT F_inv(nsd), cauchy(nsd), PK1(nsd);

	/* matrix alias to NEEvec */
	dArrayT NEEvec(NumSD()*NumElementNodes());
	NEEvec = 0.0;
	dMatrixT WP(nsd, fAmm_geo2.Rows(), NEEvec.Pointer());

	for (int i = 0; i < fSurfaceElements.Length(); i++)
	{
		/* bulk element information */
		int element = fSurfaceElements[i];
		const ElementCardT& element_card = ElementCard(element);
		fLocInitCoords.SetLocal(element_card.NodesX()); /* reference coordinates over bulk element */
		fLocDisp.SetLocal(element_card.NodesU()); /* displacements over bulk element */
	
		/* integrate surface contribution to nodal forces */
		Rmech2 = 0.0;
		for (int j = 0; j < fSurfaceElementNeighbors.MinorDim(); j++) /* loop over faces */
			if (fSurfaceElementNeighbors(i,j) == -1) /* no neighbor => surface */
//			if (fSurfaceElementNeighbors(i,j) == -1 && (fSurfaceElementFacesType(i,j) == 2 || fSurfaceElementFacesType(i,j) == 3)) /* no neighbor => surface */
			{
				/* face parent domain */
				const ParentDomainT& surf_shape = shape.FacetShapeFunction(j);
			
				/* collect coordinates of face nodes */
				ElementCardT& element_card = ElementCard(fSurfaceElements[i]);
				shape.NodesOnFacet(j, face_nodes_index);	// fni = 4 nodes of surface face
				face_nodes.Collect(face_nodes_index, element_card.NodesX());
				face_coords.SetLocal(face_nodes);
				face_curr_coords.SetLocal(face_nodes);

				/* integrate over the face */
				int face_ip;
				fSurfaceCBSupport->SetCurrIP(face_ip);
				const double* w = surf_shape.Weight();	
				int count = 0;		
				for (face_ip = 0; face_ip < nsi; face_ip++) {

				/* MAPPING/DEFORMATION */
						
					/* reference coordinate mapping on face */
					surf_shape.DomainJacobian(face_coords, face_ip, jacobian);
					double detj = surf_shape.SurfaceJacobian(jacobian);
				
					/* face normal (for all face IPs) */
					surf_shape.DomainJacobian(face_curr_coords, face_ip, jacobian);
					surf_shape.SurfaceJacobian(jacobian, Q);	// Last column of Q is normal vector to surface
					
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
					if (J <= 0.0)
						ExceptionT::BadJacobianDet(caller);
					else
						F_inv.Inverse(F);
					
					/* stress at the surface */
//					fSurfaceCB->s_ij().ToMatrix(cauchy);
					cauchy = Surf_Tension_Stress();

					/* compute PK1/J */
//					PK1.MultABT(cauchy, F_inv);
					
					/* Wi,J PiJ */
					shape.GradNa(DNa_X, fGradNa);
//					WP.MultAB(PK1, fGradNa);

					/* accumulate */
//					Rmech2.AddScaled(-J*constK*w[face_ip]*detj, NEEvec);
				}			
			}
// 		Rtotal2.CopyIn(0, Rmech2);

		fRHS += Rtotal2;
		
		/* assemble forces */
//		ElementSupport().AssembleRHS(Group(), fRHS, element_card.Equations());		

	}		
  }
/***********************************************************************
 * Protected
 ***********************************************************************/
dMatrixT FSDielectricElastomerQ1P0SurfaceT::Surf_BMatrix(const dMatrixT& B_original, const dMatrixT& Q)
{
	dMatrixT surfb(NumSD(), NumElementNodes());
	dMatrixT surfb1(NumSD(), NumElementNodes());
	dMatrixT surfb2(NumSD());
	dArrayT normal(NumSD());
	
	/* extract normal from Q */
	Q.CopyColumn(2,normal);

	/* modify B matrix using surface gradient */
	surfb = B_original;
	surfb2.Outer(normal, normal);
	surfb1.MultAB(surfb2, B_original);
	surfb-=surfb1;
	
	return surfb;
}

dMatrixT FSDielectricElastomerQ1P0SurfaceT::Surf_Tension_Stress()
{
	// MAY NOT NEED THIS - JUST MULTIPLY SURF B-MATRIX BY GAMMA
	dMatrixT tempstress(NumSD());
	tempstress.Identity();
	tempstress*=fSurfTension;

	return tempstress;
}

dMatrixT FSDielectricElastomerQ1P0SurfaceT::Surf_Tension_Stiffness()
{
	dMatrixT tempstiff(dSymMatrixT::NumValues(NumSD()));
	
	return tempstiff;
}

dArray2DT FSDielectricElastomerQ1P0SurfaceT::MatrixTodArray2DT(const dMatrixT& source)
{
	dArray2DT B(NumSD(), NumElementNodes());
	
	for (int i = 0; i < NumSD(); i++)
	{
		for (int j = 0; j < NumElementNodes(); j++)
		{
			B(i,j) = source(i,j);
		}	
	}

	return B;
}

/***********************************************************************
 * Private
 ***********************************************************************/

  // extrapolate from integration points and compute output nodal/element values
  void FSDielectricElastomerQ1P0SurfaceT::ComputeOutput(const iArrayT& n_codes,
      dArray2DT& n_values, const iArrayT& e_codes, dArray2DT& e_values)
  {
  	/* inherited - bulk contribution */
  	FSDielectricElastomerQ1P0T::ComputeOutput(n_codes, n_values, e_codes, e_values);

   }

} // namespace Tahoe
