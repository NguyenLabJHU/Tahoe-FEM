#include "FSDielectricElastomerT.h"
#include "FSDEMatSupportT.h"
#include "FSDEMatT.h"
#include "ParameterContainerT.h"
#include "OutputSetT.h"
#include "ShapeFunctionT.h"
#include "eIntegratorT.h"

/* REMAINING ISSUES (NOVEMBER 7, 2010):
2.  Shape function gradients in LHS
3.  Shape function gradients in RHS
4.  Check stress, stiffness, electric displacement and electromechanical coupling 
5.  SetGlobalShape:  calculating Efield from gradient of Psi?
6.  Related to 6:  why can't call fShapes->Derivatives_X in SetGlobalShape?
7.  Some reasonable boundary value problems (2D/3D)?
8.  External "forces" (i.e. charges) in electrical system:
	does this need to be applied here or at an element level?
9.  IP weight has different signs in FormStiffness and FormKd 
*/

//
// materials lists (3D only)
//
#include "FSSolidMatList2DT.h"
#include "FSSolidMatList3DT.h"

namespace Tahoe {


/* Destructor */
  FSDielectricElastomerT::~FSDielectricElastomerT()
  {
	  if (0 != fFSDEMatSupport) delete fFSDEMatSupport;
  }

  //
  // specify parameters needed by the interface
  //
  void FSDielectricElastomerT::DefineParameters(ParameterListT& list) const
  {
//  	cout << "FSDielectricElastomerT::DefineParameters" << endl;
    // inherited
    FiniteStrainT::DefineParameters(list);

    // additional fields
    list.AddParameter(ParameterT::Word, "electric_field_name");
  }

  //
  // accept parameter list
  //
  void FSDielectricElastomerT::TakeParameterList(const ParameterListT& list)
  {
//  	cout << "FSDielectricElastomerT::TakeParameterList" << endl;
    //
    // inherited
    //
    FiniteStrainT::TakeParameterList(list);

    //
    // get electric scalar potential field
    //
    // for now use same integration and interpolation schemes as primary field
    //
    const StringT& electric_field_name = list.GetParameter(
        "electric_field_name");

    fElectricScalarPotentialField = ElementSupport().Field(electric_field_name);
    if (!fElectricScalarPotentialField) {
      ExceptionT::GeneralFail("FSDielectricElastomerT::TakeParameterList",
          "could not resolve \"%s\" field", electric_field_name.Pointer());
    }

    /* Define matrix sizes */
    int nen = NumElementNodes();
    int nsd = NumSD();
    int nel = nen;	// # electrical DOFs per element
    int nme = nen * nsd;	// # of mechanical DOFs per element
    int dof = nsd + 1;	// total # of DOFs per node (mech + elec)
    int neq = nen * dof;	// total # of DOFs per element (mech + elec)

	/* Dimension electric field arrays */
	const int nip = NumIP();
	fE_all.Dimension(nip*nsd);
	fE_all = 0.0;	// testing HSP
	fE_List.Dimension(nip);
	
	/* what does this do? */
    for (int i = 0; i < nip; ++i) {
      fE_List[i].Alias(nsd, fE_all.Pointer(i * nsd));
    }	

    fAmm_mat.Dimension(nme, nme);
    fAmm_geo.Dimension(nme, nme);
//	fAmm_geo.Dimension(nen, nen);
    fAme.Dimension(nme, nel);
    fAem.Dimension(nel, nme);
    fAee.Dimension(nel, nel);

    fLHS.Dimension(neq);
    fRHS.Dimension(neq);
    
// 	/* dimension workspace */
// 	fStressMat.Dimension(NumSD());
// 	fTempMat1.Dimension(NumSD());
// 	fTempMat2.Dimension(NumSD());
// 
// 	fGradNa.Dimension(NumSD(), NumElementNodes());
// 	fStressStiff.Dimension(NumElementNodes());
// 	
// 	/* Might need to re-dimension, or it does mechanical RHS only */
// 	fTemp2.Dimension(NumElementNodes()*NumDOF());    
    
  }

  //
  // PROTECTED
  //

  //
  // construct a new material support and return a pointer
  //
  MaterialSupportT*
  FSDielectricElastomerT::NewMaterialSupport(MaterialSupportT* p) const
  {
//  	cout << "FSDIelectricElastomerT::NewMaterialSupport" << endl;
    //
    // allocate
    //
    if (!p) p = new FSDEMatSupportT(1, NumIP());
    
    //
    // inherited initializations
    //
    FiniteStrainT::NewMaterialSupport(p);

    //
    // set parent class fields
    //
    FSDEMatSupportT* ps = dynamic_cast<FSDEMatSupportT*> (p);

    if (ps != 0) {
//	  cout << "Setting Electric Field from FSDielectricElastomerT" << endl;
      ps->SetElectricField(&fE_List);
    }

    return p;
  }

  //
  // construct materials manager and read data
  //
  MaterialListT*
  FSDielectricElastomerT::NewMaterialList(const StringT& name, int size)
  {
//  	cout << "FSDIelectricElastomerT::NewMaterialList" << endl;
  	
  	/* resolve number of spatial dimensions */
  	int nsd = -1;
  	if (name == "large_strain_material_2D")
  		nsd = 2;
  	else if (name == "large_strain_material_3D")
  		nsd = 3;
  	
  	/* no match */
  	if (nsd == -1) return NULL;
  	
    MaterialListT* mlp = 0;

    if (size > 0) {

   	  	/* material support */
      	if (0 == fFSDEMatSupport) {
        	fFSDEMatSupport = dynamic_cast<FSDEMatSupportT*> (NewMaterialSupport());

        if (0 == fFSDEMatSupport) {
          ExceptionT::GeneralFail("FSDielectricElastomerT::NewMaterialList");
        }
      }
      if (nsd == 2)
      	mlp = new FSSolidMatList2DT(size, *fFSDEMatSupport);
	  else if (nsd == 3)
		mlp = new FSSolidMatList3DT(size, *fFSDEMatSupport); 
    }
	else {
		if (nsd == 2)
			mlp = new FSSolidMatList2DT;
	 	else if (nsd == 3)
	 		mlp = new FSSolidMatList3DT;
    } 
    return mlp;

  }

  //
  // form shape functions and derivatives; for ComputeOutput
  // IS THIS NEEDED FOR VALUES OF E?
  void FSDielectricElastomerT::SetGlobalShape()
  {
//	cout << "FSDielectricElastomerT::SetGlobalShape" << endl;
    //
    // inherited
    //
    FiniteStrainT::SetGlobalShape();

    //
    // what needs to be computed
    //
    SetLocalU(fLocScalarPotential);
//	cout << "fLocScalarPotential = " << fLocScalarPotential << endl;
    for (int i = 0; i < NumIP(); i++) {

      //
      // electric field
      //

        dArrayT& E = fE_List[i];
		dMatrixT E1(1, NumSD());

		/* Not sure if this is correct; would prefer Derivative_X */
		/* Sometimes z value of E1 is e-19....*/
		fShapes->GradU(fLocScalarPotential, E1, i);
		E1 *= -1.0;
		for (int i = 0; i < NumSD(); i++)
			E[i] = E1(0,i);
		
//		cout << "E-field derived from Psi = " << E << endl;
      }

  }

  //
  // write all current element information to the stream
  //
  void FSDielectricElastomerT::CurrElementInfo(ostream& out) const
  {
//	cout << "FSDielectricElastomerT::CurrElementInfo" << endl;
    //
    // inherited
    //
    FiniteStrainT::CurrElementInfo(out);

    //
    // write deformation gradients
    //
    out << std::endl;
    out << "electric field at IP:";
    out << std::endl;

    for (int i = 0; i < fE_List.Length(); ++i) {
      out << " ip: " << i + 1 << std::endl << fE_List[i] << std::endl;
    }

    out << std::endl;

  }

// 
//  increment current element - for ComputeOutput
// 
  bool FSDielectricElastomerT::NextElement()
  {
//	cout << "FSDielectricElastomerT::NextElement" << endl;
    bool isThereNext = FiniteStrainT::NextElement();

    if (isThereNext == true) {

      const int index = CurrentElement().MaterialNumber();

      ContinuumMaterialT* pMaterial = (*fMaterialList)[index];

      fCurrMaterial = dynamic_cast<FSDEMatT*> (pMaterial);
    }

    return isThereNext;
  }

  //
  // Initialize local arrays
  //
  void FSDielectricElastomerT::SetLocalArrays()
  {
//	cout << "FSDielectricElastomerT::SetLocalArrays" << endl;
    //
    // look for an electric scalar potential field
    //
    const FieldT* esp = 0;

    if (0 == fElectricScalarPotentialField) {
      esp = ElementSupport().Field("electric_scalar_potential");
      fElectricScalarPotentialField = esp;
    } else {
      esp = fElectricScalarPotentialField;
    }

    if (0 == esp) {

      std::cout << std::endl;
      std::cout << "FSDielectricElastomerT::SetLocalArrays: ";
      std::cout << "Electric scalar potential field not found.";
      std::cout << std::endl;

      throw ExceptionT::kGeneralFail;

    }

    //
    // Inherited
    //
    FiniteStrainT::SetLocalArrays();

    //
    // Allocate storage space:  1 DOF/node for scalar potential
    //
    const int nen = NumElementNodes();

    fLocScalarPotential.Dimension(nen, 1);

    //
    // Register fields
    //
    esp->RegisterLocal(fLocScalarPotential);

  }


  //
  //
  //
  void FSDielectricElastomerT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
      AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
  {
    for (int i = 0; i < fEqnos.Length(); ++i) {

      const int ndf = NumDOF();	// scalar potential
      const int nen = fConnectivities[i]->MinorDim();
      const int offset = ndf * nen;

      fElectricScalarPotentialField->SetLocalEqnos(*fConnectivities[i],
          fEqnos[i], offset);

    }
	
    ElementBaseT::Equations(eq_1, eq_2);
  }

/* accumulate the residual force on the specified node */
void FSDielectricElastomerT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
	/* not my field */
	if (&field != &(Field())) return;

	/* quick exit */
	bool hasnode = false;
	for (int i=0; i < fBlockData.Length() && !hasnode; i++)
		if (fConnectivities[i]->HasValue(node)) hasnode = true;
	if (!hasnode) return;

	/* set components and weights */
	double constMa = 0.0;
	double constKd = 0.0;

	/* components dicated by the algorithm */
	int formMa = fIntegrator->FormMa(constMa);
	int formKd = fIntegrator->FormKd(constKd);

	/* body forces */
	int formBody = 0;
	if (fMassType != kNoMass &&
	   (fBodySchedule && fBody.Magnitude() > kSmall))
	{
		cout << "\nWarning: Body forces not yet implemented in DielectricElastomerT";
		if (!formMa) constMa = 1.0; /* override */
	}

	/* override controller */
	if (fMassType == kNoMass) formMa = 0;

	/* temp for nodal force */
	dArrayT nodalforce;

	bool axisymmetric = Axisymmetric();
	Top();
	while (NextElement())
	{
		int nodeposition;
		const iArrayT& nodes_u = CurrentElement().NodesU();
		if (nodes_u.HasValue(node, nodeposition))
		{
			/* initialize */
			fRHS = 0.0;

			/* global shape function values */
			SetGlobalShape();

			/* internal force contribution */
			if (formKd) FormKd(constKd);

			/* inertia forces */
			if (formMa)
			{
				SetLocalU(fLocAcc);
				FormMa(fMassType, constMa*fCurrMaterial->Density(), axisymmetric, &fLocAcc, NULL, NULL);
			}
	
			double mr1, mr2, mr3, er1;
			dArrayT react(4);
			
			/* loop over nodes (double-noding OK) */
			int dex = 0;
			int dex2 = 0;
			for (int i = 0; i < nodes_u.Length(); i++)
			{
				if (nodes_u[i] == node)
				{
					mr1 = fRHS[dex];
					mr2 = fRHS[dex+1];
					mr3 = fRHS[dex+2];
					er1 = fRHS[dex2+3*NumElementNodes()];
					react[0] = mr1;
					react[1] = mr2;
					react[2] = mr3;
					react[3] = er1;
					
					/* components for node - mechanical + electrical DOFs */
//					nodalforce.Set(TotalNumDOF(), fRHS.Pointer(dex+2*NumElementNodes()));
					nodalforce.Set(TotalNumDOF(), react.Pointer(0));

					/* accumulate */
					force += nodalforce;
				}
				dex += NumDOF();
				dex2 += 1;
			}
		}
	}
}


  //
  // Strain displacement operator
  //
  void FSDielectricElastomerT::Set_B_C(const dArray2DT& DNaX, dMatrixT& B_C)
  {

    const int nsd = NumSD();
    const int nen = NumElementNodes();
    const int StrainDim = dSymMatrixT::NumValues(nsd);
    const int nnd = DNaX.MinorDim();

    assert(nen == nnd);

    B_C.Dimension(StrainDim, nsd * nen);
    const dMatrixT F = DeformationGradient();

    for (int ij = 0; ij < StrainDim; ++ij) {

      int i;
      int j;

      dSymMatrixT::ExpandIndex(nsd, ij, i, j);

      for (int a = 0; a < nen; ++a) {

        for (int k = 0; k < nsd; ++k) {

          const int ak = a * nsd + k;

          B_C(ij, ak) = DNaX(i, a) * F(k, j) + DNaX(j, a) * F(k, i);

          // Shear components doubled to conform with Voigt
          // convention
          if (ij >= nsd) {
            B_C(ij, ak) *= 2.0;
          }

        }

      }

    }

  }

  //
  // Geometric stiffness
  //
  void FSDielectricElastomerT::AccumulateGeometricStiffness(dMatrixT& Kg,
      const dArray2DT& DNaX, dSymMatrixT& S)
  {

    const int nsd = NumSD();
    const int nen = NumElementNodes();

    dMatrixT GradNa(nsd, nen);
    fShapes->GradNa(DNaX, GradNa);

    dMatrixT DNaS(nen);

    dMatrixT Sm(nsd);
    S.ToMatrix(Sm);

    DNaS.MultQTBQ(GradNa, Sm);

    dMatrixT SI(nsd);

    for (int i = 0; i < nen; ++i) {

      for (int j = 0; j < nen; ++j) {

        SI.Identity(DNaS(i, j));

        Kg.AddBlock(nsd * i, nsd * j, SI);

      }

    }

  }

/* calculate the LHS of residual, or element stiffness matrix */
  void FSDielectricElastomerT::FormStiffness(double constK)
  {
	/* Matrix format */
    dMatrixT::SymmetryFlagT format = (fLHS.Format()
        == ElementMatrixT::kNonSymmetric)
        ? dMatrixT::kWhole
        : dMatrixT::kUpperOnly;

	/* Element preliminaries */
    const int nsd = NumSD();
    const int nen = NumElementNodes();

    fAmm_mat = 0.0;
    fAmm_geo = 0.0;
    fAme = 0.0;
    fAem = 0.0;
    fAee = 0.0;

	/* definitions for new moduli */
	dMatrixT Cmech(2*NumSD(), 2*NumSD());
	dMatrixT Celecmech(2*NumSD(), NumSD());
	dArrayT  DE(NumSD());
	dMatrixT Celec(NumSD(), NumSD());
	dMatrixT GradShape(nsd, nen);

    fShapes->TopIP();
    while (fShapes->NextIP() != 0) 
    {
		/* integration weight */
		const double w = constK * fShapes->IPDet() * fShapes->IPWeight();

	  	/* LHS tangent stiffnesses */  
   	    fCurrMaterial->C_Mech_Elec(Cmech, Celecmech);
    	fCurrMaterial->S_C_Elec(DE, Celec);
      	dSymMatrixT S = fCurrMaterial->S_IJ();
//       	cout << "Cmech FormStiffness = " << Cmech << endl;
//       	cout << "Celecmech FormStiffness = " << Celecmech << endl;
//       	cout << "S FormStiffness = " << S << endl;
//       	cout << "Celec FormStiffness = " << Celec << endl;
        Cmech *= (0.25*w);
        Celecmech *= (0.5*w);
	  	Celec *= w;
        S *= w;

  	    /* prepare derivatives of shape functions in reference configuration */
    	const dArray2DT& DNaX = fShapes->Derivatives_X();
	  	fShapes->GradNa(DNaX, GradShape);	  

	  	/* B_C comes back 6 x 24 shape function gradient matrix */
	  	dMatrixT B_C;
	  	Set_B_C(DNaX, B_C);
 
	  	/* mechanical material stiffness (24 x 24 matrix for 8-node 3D element) */
 	  	fAmm_mat.MultQTBQ(B_C, Cmech, format, dMatrixT::kAccumulate);
	
	  	/* mechanical geometric stiffness (24 x 24 matrix for 8-node 3D element */ 
	  	AccumulateGeometricStiffness(fAmm_geo, DNaX, S);
  
       	/* mechanical-electrical stiffness (24 x 8 matrix for 8-node 3D element) */
       	/* What is the difference between format and dMatrixT::kWhole? */
       	fAme.MultATBC(B_C, Celecmech, GradShape, dMatrixT::kWhole, dMatrixT::kAccumulate); 
  
 	 	/* electrical-electrical stiffness (8 x 8 matrix for 8-node 3D element) */
  		fAee.MultQTBQ(GradShape, Celec, format, dMatrixT::kAccumulate);
    }

	/* Expand 24x24 geometric stiffness into material stiffness matrix */
	fAmm_mat.Expand(fAmm_geo, 1, dMatrixT::kAccumulate);

	/* Assemble into fLHS, or element stiffness matrix */
	fLHS.AddBlock(0, 0, fAmm_mat);
	fLHS.AddBlock(fAmm_mat.Rows(), fAmm_mat.Cols(), fAee);
	fLHS.AddBlock(0, fAmm_mat.Cols(), fAme);
  }


/* Compute RHS, or residual of element equations */
  void FSDielectricElastomerT::FormKd(double constK)
  {
	/* element preliminaries */
    const int nsd = NumSD();
    const int nen = NumElementNodes();
    
    /* Define mechanical and electrical residuals */
	dArrayT Rtotal((nsd+1)*nen);
	Rtotal = 0.0;
	dArrayT Rmech(nen*nsd);
	Rmech = 0.0;
	dArrayT Relec(nen);
	Relec = 0.0;
	dMatrixT GradShape(nsd, nen);
	dArrayT DE(NumSD());
	dMatrixT Celec(NumSD(), NumSD());

    fShapes->TopIP();
    while (fShapes->NextIP() != 0) 
    {

	  /* integration weight */
      const double w = constK * fShapes->IPDet() * fShapes->IPWeight();
      const dArray2DT& DNaX = fShapes->Derivatives_X();
      
      /* Now convert DNaX to a matrix instead of dArray2DT */
	  fShapes->GradNa(DNaX, GradShape);	
	  
	  /* Mechanical stress */
	  dSymMatrixT S = fCurrMaterial->S_IJ();
	  S *= (0.5*w);
	  dMatrixT B_C;
	  Set_B_C(DNaX, B_C);
	  B_C.MultTx(S, Rmech, 1.0, dMatrixT::kAccumulate);
	  
	  /* electrical stress */
	  fCurrMaterial->S_C_Elec(DE, Celec);
	  DE *= w;
	  
	  /* 3x1 vector of shape function gradient * D */
	  GradShape.MultTx(DE, Relec, 1.0, dMatrixT::kAccumulate);	
	  
	  /* NOTE:  mechanical inertia term, mechanical body force term, mechanical
	  surface traction term, electrical body force (charge), electrical surface 
	  traction not accounted for here */
	  
	}
  	Relec *= -1.0;	// need for right sign for residual
 	Rtotal.CopyIn(0, Rmech);
 	Rtotal.CopyIn(Rmech.Length(), Relec);
 	fRHS += Rtotal;
  }

/***********************************************************/
/* BELOW IS LIKE TOTALLAGRANGIANT.CPP */
/* calculate the LHS of residual, or element stiffness matrix */
//   void FSDielectricElastomerT::FormStiffness(double constK)
//   {
// 	/* Matrix format */
//     dMatrixT::SymmetryFlagT format = (fLHS.Format()
//         == ElementMatrixT::kNonSymmetric)
//         ? dMatrixT::kWhole
//         : dMatrixT::kUpperOnly;
// 
// 	/* integrate over the element */
//     const int nsd = NumSD();
//     const int nen = NumElementNodes();
// 
// 	/* Initialize portions of LHS */
//     fAmm_mat = 0.0;
//     fAmm_geo = 0.0;
//     fAme = 0.0;
//     fAem = 0.0;
//     fAee = 0.0;
//     fStressStiff = 0.0;
// 
// 	/* definitions for new moduli */
// 	dMatrixT Cmech(2*NumSD(), 2*NumSD());
// 	dMatrixT D(6);
// 	dMatrixT B;
// 	B.Dimension(6, nsd*nen);
// 
// 	/* integration */
// 	const double* Det    = fShapes->IPDets();
// 	const double* Weight = fShapes->IPWeights();
// 
//     fShapes->TopIP();
//     while (fShapes->NextIP() != 0) 
//     {
// 	/* STRESS STIFFNESS */
// 	
// 		/* Cauchy stress (and set deformation gradient) */
// 		(fCurrMaterial->s_ij()).ToMatrix(fStressMat);
// 
// 		/* chain rule shape function derivatives */
// 		fTempMat1 = DeformationGradient();
// 		double J = fTempMat1.Det();
// 		fTempMat1.Inverse();
// 		fShapes->TransformDerivatives(fTempMat1, fDNa_x);
// 
// 		/* get shape function gradients matrix */
// 		fShapes->GradNa(fDNa_x, fGradNa);	
// 				
// 		/* scale factor */
// 		double scale = constK*(*Det++)*(*Weight++)*J;
// 
// 		/* integration constants */		
// 		fStressMat *= scale;
// 	
// 		/* using the stress symmetry */
// 		fStressStiff.MultQTBQ(fGradNa, fStressMat, format,
// 			dMatrixT::kAccumulate);
// 			
// 	/* MATERIAL STIFFNESS */
// 		/* strain displacement matrix */
// 		Set_B(fDNa_x, B);
// 		
// 		/* get D matrix */
// 		D.SetToScaled(scale, fCurrMaterial->c_ijkl());
// 		
// 		/* accumulate */
// 		fAmm_mat.MultQTBQ(B, D, format, dMatrixT::kAccumulate);
//     }
// 	/* stress stiffness into material stiffness */
// 	fAmm_mat.Expand(fStressStiff, NumDOF(), dMatrixT::kAccumulate);
//  	fLHS.AddBlock(0, 0, fAmm_mat);
// 
//   }

/* Compute RHS, or residual of element equations */
//   void FSDielectricElastomerT::FormKd(double constK)
//   {
//     const int nsd = NumSD();
//     const int nen = NumElementNodes();
//     
//     /* Define mechanical and electrical residuals */
// 	dArrayT Rtotal((nsd+1)*nen);
// 	Rtotal = 0.0;
// 	dArrayT Rmech(nen*nsd), Rmech2(nen*nsd);
// 	Rmech = 0.0;
// 	Rmech2 = 0.0;
// 	dArrayT Relec(nen);
// 	Relec = 0.0;
// 	dArrayT DE(NumSD());
// 	dMatrixT Celec(NumSD(), NumSD());
// 
// 	/* matrix alias to ??? */
// 	dMatrixT fWP(NumSD(), fStressStiff.Rows(), Rmech2.Pointer());
// 	
// 	const double* Det    = fShapes->IPDets();
// 	const double* Weight = fShapes->IPWeights();
// 
//     fShapes->TopIP();
//     while (fShapes->NextIP() != 0) 
//     {
// 		/* get Cauchy stress */
// 		(fCurrMaterial->s_ij()).ToMatrix(fTempMat1);
// 
// 		/* F^(-1) */
// 		fTempMat2 = DeformationGradient();
// 	  
// 	  	double J = fTempMat2.Det();
// 	  	if (J <= 0.0)
// 	  		ExceptionT::BadJacobianDet("TotalLagrangianT::FormKd");
// 	  	else
// 	  		fTempMat2.Inverse();
// 	  
// 	  	/* compute PK1/J */
// 	  	fStressMat.MultABT(fTempMat1, fTempMat2);
// 	  
// 	  	/* get matrix of shape function gradients */
// 	  	fShapes->GradNa(fGradNa);
// 	  	
// 	  	/* Wi,J PiJ */
// 	  	fWP.MultAB(fStressMat, fGradNa);
// 	  	
// 	  	/* accumulate */
// 	  	Rmech.AddScaled(J*constK*(*Weight++)*(*Det++), Rmech2);		  
// 	}
// //	cout << "Rmech =  " << Rmech << endl;
// //  	Relec *= -1.0;	// need for right sign for residual
//  	Rtotal.CopyIn(0, Rmech);
// // 	Rtotal.CopyIn(Rmech.Length(), Relec);
//  	fRHS += Rtotal;
//   }

  //
  // extrapolate from integration points and compute output nodal/element
  // values
  //
  void FSDielectricElastomerT::ComputeOutput(const iArrayT& n_codes,
      dArray2DT& n_values, const iArrayT& e_codes, dArray2DT& e_values)
  {
    //
    // number of output values
    //
    int n_out = n_codes.Sum();
    int e_out = e_codes.Sum();

    // nothing to output
    if (n_out == 0 && e_out == 0) return;

    // dimensions
    int nsd = NumSD();
    int ndof = NumDOF();
    int nen = NumElementNodes();
    int nnd = ElementSupport().NumNodes();

    // reset averaging work space
    ElementSupport().ResetAverage(n_out);

    // allocate element results space
    e_values.Dimension(NumElements(), e_out);

    // nodal work arrays
    dArray2DT nodal_space(nen, n_out);
    dArray2DT nodal_all(nen, n_out);
    dArray2DT coords, disp;
    dArray2DT nodalstress, princstress, matdat;
    dArray2DT energy, speed;
    dArray2DT ndElectricScalarPotential;
//	dArrayT ndElectricScalarPotential;
    dArray2DT ndElectricDisplacement;
    dArray2DT ndElectricField;

    // ip values
    dArrayT ipmat(n_codes[iMaterialData]), ipenergy(1);
    dArrayT ipspeed(nsd), ipprincipal(nsd);
    dMatrixT ippvector(nsd);

    // set shallow copies
    double* pall = nodal_space.Pointer();
    coords.Alias(nen, n_codes[iNodalCoord], pall);
    pall += coords.Length();
    disp.Alias(nen, n_codes[iNodalDisp], pall);
    pall += disp.Length();

    nodalstress.Alias(nen, n_codes[iNodalStress], pall);
    pall += nodalstress.Length();
    princstress.Alias(nen, n_codes[iPrincipal], pall);
    pall += princstress.Length();
    energy.Alias(nen, n_codes[iEnergyDensity], pall);
    pall += energy.Length();
    speed.Alias(nen, n_codes[iWaveSpeeds], pall);
    pall += speed.Length();
    matdat.Alias(nen, n_codes[iMaterialData], pall);
    pall += matdat.Length();

    ndElectricDisplacement.Alias(nen, n_codes[ND_ELEC_DISP], pall);
    pall += ndElectricDisplacement.Length();

    ndElectricField.Alias(nen, n_codes[ND_ELEC_FLD], pall);
    pall += ndElectricField.Length();
    
    ndElectricScalarPotential.Alias(nen, n_codes[ND_ELEC_POT_SCALAR], pall);
    pall += ndElectricScalarPotential.Length();

    // element work arrays
    dArrayT element_values(e_values.MinorDim());
    pall = element_values.Pointer();
    dArrayT centroid, ip_centroid, ip_mass;
    dArrayT ip_coords(nsd);
    if (e_codes[iCentroid]) {
      centroid.Alias(nsd, pall);
      pall += nsd;
      ip_centroid.Dimension(nsd);
    }
    if (e_codes[iMass]) {
      ip_mass.Alias(NumIP(), pall);
      pall += NumIP();
    }
    double w_tmp, ke_tmp;
    double mass;
    double& strain_energy = (e_codes[iStrainEnergy])
        ? *pall++
        : w_tmp;
    double& kinetic_energy = (e_codes[iKineticEnergy])
        ? *pall++
        : ke_tmp;
    dArrayT linear_momentum, ip_velocity;

    if (e_codes[iLinearMomentum]) {
      linear_momentum.Alias(ndof, pall);
      pall += ndof;
      ip_velocity.Dimension(ndof);
    } else if (e_codes[iKineticEnergy]) ip_velocity.Dimension(ndof);

    dArray2DT ip_stress;
    if (e_codes[iIPStress]) {
      ip_stress.Alias(NumIP(), e_codes[iIPStress] / NumIP(), pall);
      pall += ip_stress.Length();
    }
    dArray2DT ip_material_data;
    if (e_codes[iIPMaterialData]) {
      ip_material_data.Alias(NumIP(), e_codes[iIPMaterialData] / NumIP(), pall);
      pall += ip_material_data.Length();
      ipmat.Dimension(ip_material_data.MinorDim());
    }

    dArray2DT ipElectricDisplacement;
    if (e_codes[IP_ELEC_DISP]) {
      ipElectricDisplacement.Alias(NumIP(), NumSD(), pall);
      pall += NumIP() * NumSD();
    }

    dArray2DT ipElectricField;
    if (e_codes[IP_ELEC_FLD]) {
      ipElectricField.Alias(NumIP(), NumSD(), pall);
      pall += NumIP() * NumSD();
    }

    // check that degrees are displacements
    int interpolant_DOF = InterpolantDOFs();

    Top();
    while (NextElement()) {

      if (CurrentElement().Flag() == ElementCardT::kOFF) continue;

      // initialize
      nodal_space = 0.0;

      // global shape function values
      SetGlobalShape();

      // collect nodal values
      if (e_codes[iKineticEnergy] || e_codes[iLinearMomentum]) {
        if (fLocVel.IsRegistered())
          SetLocalU(fLocVel);
        else
          fLocVel = 0.0;
      }

      // coordinates and displacements all at once
      if (n_codes[iNodalCoord]) fLocInitCoords.ReturnTranspose(coords);
      if (n_codes[iNodalDisp]) {
        if (interpolant_DOF)
          fLocDisp.ReturnTranspose(disp);
        else
          NodalDOFs(CurrentElement().NodesX(), disp);
      }

      if (n_codes[ND_ELEC_POT_SCALAR]) {
        if (interpolant_DOF) {
          fLocScalarPotential.ReturnTranspose(ndElectricScalarPotential);
        } else {
          NodalDOFs(CurrentElement().NodesX(), ndElectricScalarPotential);
        }
      }

      // initialize element values
      mass = strain_energy = kinetic_energy = 0;
      if (e_codes[iCentroid]) centroid = 0.0;
      if (e_codes[iLinearMomentum]) linear_momentum = 0.0;
      const double* j = fShapes->IPDets();
      const double* w = fShapes->IPWeights();

      // integrate
      dArray2DT Na_X_ip_w;
      fShapes->TopIP();
      while (fShapes->NextIP() != 0) {

        // density may change with integration point
        double density = fCurrMaterial->Density();

        // element integration weight
        double ip_w = (*j++) * (*w++);

        if (qNoExtrap) {
          Na_X_ip_w.Dimension(nen, 1);
          for (int k = 0; k < nen; k++) {
            Na_X_ip_w(k, 0) = 1.;
          }
        }

        // get Cauchy stress
        const dSymMatrixT& stress = fCurrMaterial->s_ij();
        dSymMatrixT strain;

        // stresses
        if (n_codes[iNodalStress]) {
          if (qNoExtrap) {
            for (int k = 0; k < nen; k++) {
              nodalstress.AddToRowScaled(k, Na_X_ip_w(k, 0), stress);
            }
          } else {
            fShapes->Extrapolate(stress, nodalstress);
          }
        }

        if (e_codes[iIPStress]) {
          double* row = ip_stress(fShapes->CurrIP());
          strain.Set(nsd, row);
          strain = stress;
          row += stress.Length();
          strain.Set(nsd, row);
          fCurrMaterial->Strain(strain);
        }

        // wave speeds
        if (n_codes[iWaveSpeeds]) {
          // acoustic wave speeds
          fCurrMaterial->WaveSpeeds(fNormal, ipspeed);
          if (qNoExtrap) {
            for (int k = 0; k < nen; k++) {
              speed.AddToRowScaled(k, Na_X_ip_w(k, 0), ipspeed);
            }
          } else {
            fShapes->Extrapolate(ipspeed, speed);
          }
        }

        // principal values - compute principal before smoothing
        if (n_codes[iPrincipal]) {
          // compute eigenvalues
          stress.PrincipalValues(ipprincipal);

          if (qNoExtrap) {
            for (int k = 0; k < nen; k++) {
              princstress.AddToRowScaled(k, Na_X_ip_w(k, 0), ipprincipal);
            }
          } else {
            fShapes->Extrapolate(ipprincipal, princstress);
          }
        }

        // strain energy density
        if (n_codes[iEnergyDensity] || e_codes[iStrainEnergy]) {
          double ip_strain_energy = fCurrMaterial->StrainEnergyDensity();

          // nodal average
          if (n_codes[iEnergyDensity]) {
            ipenergy[0] = ip_strain_energy;
            if (qNoExtrap) {
              for (int k = 0; k < nen; k++) {
                energy.AddToRowScaled(k, Na_X_ip_w(k, 0), ipenergy);
              }
            } else {
              fShapes->Extrapolate(ipenergy, energy);
            }
          }

          // integrate over element
          if (e_codes[iStrainEnergy]) {
            strain_energy += ip_w * ip_strain_energy;
          }

        }

        // material stuff
        if (n_codes[iMaterialData] || e_codes[iIPMaterialData]) {
          // compute material output
          fCurrMaterial->ComputeOutput(ipmat);

          // store nodal data
          if (n_codes[iMaterialData]) {
            if (qNoExtrap) {
              for (int k = 0; k < nen; k++) {
                matdat.AddToRowScaled(k, Na_X_ip_w(k, 0), ipmat);
              }
            } else {
              fShapes->Extrapolate(ipmat, matdat);
            }
          }

          // store element data
          if (e_codes[iIPMaterialData]) {
            ip_material_data.SetRow(fShapes->CurrIP(), ipmat);
          }
        }

        // mass averaged centroid
        if (e_codes[iCentroid] || e_codes[iMass]) {
          // mass
          mass += ip_w * density;

          // integration point mass
          if (e_codes[iMass]) ip_mass[fShapes->CurrIP()] = ip_w * density;

          // moment
          if (e_codes[iCentroid]) {
            fShapes->IPCoords(ip_centroid);
            centroid.AddScaled(ip_w * density, ip_centroid);
          }
        }

        // kinetic energy/linear momentum
        if (e_codes[iKineticEnergy] || e_codes[iLinearMomentum]) {
          // velocity at integration point
          fShapes->InterpolateU(fLocVel, ip_velocity);
          double ke_density = 0.5 * density * dArrayT::Dot(ip_velocity,
              ip_velocity);

          // kinetic energy
          if (e_codes[iKineticEnergy]) {
            kinetic_energy += ip_w * ke_density;
          }

          // linear momentum
          if (e_codes[iLinearMomentum]) {
            linear_momentum.AddScaled(ip_w * density, ip_velocity);
          }

        }

        // electric displacements
        const dArrayT& D = fCurrMaterial->D_I();
        if (n_codes[ND_ELEC_DISP]) {
          if (qNoExtrap) {
            for (int k = 0; k < nen; k++) {
              ndElectricDisplacement.AddToRowScaled(k, Na_X_ip_w(k, 0), D);
            }
          } else {
            fShapes->Extrapolate(D, ndElectricDisplacement);
          }
        }

        // electric field
        const dArrayT& E = fCurrMaterial->E_I();
        if (n_codes[ND_ELEC_FLD]) {
          if (qNoExtrap) {
            for (int k = 0; k < nen; k++) {
              ndElectricField.AddToRowScaled(k, Na_X_ip_w(k, 0), E);
            }
          } else {
            fShapes->Extrapolate(E, ndElectricField);
          }
        }

      }

      // copy in the cols
      int colcount = 0;
      nodal_all.BlockColumnCopyAt(disp, colcount);
      colcount += disp.MinorDim();

      nodal_all.BlockColumnCopyAt(coords, colcount);
      colcount += coords.MinorDim();

      if (qNoExtrap) {
        double nip(fShapes->NumIP());
        nodalstress /= nip;
        princstress /= nip;
        energy /= nip;
        speed /= nip;
        matdat /= nip;
        ndElectricDisplacement /= nip;
        ndElectricField /= nip;
        ndElectricScalarPotential /= nip;
      }
      nodal_all.BlockColumnCopyAt(nodalstress, colcount);
      colcount += nodalstress.MinorDim();

      nodal_all.BlockColumnCopyAt(princstress, colcount);
      colcount += princstress.MinorDim();

      nodal_all.BlockColumnCopyAt(energy, colcount);
      colcount += energy.MinorDim();

      nodal_all.BlockColumnCopyAt(speed, colcount);
      colcount += speed.MinorDim();

      nodal_all.BlockColumnCopyAt(matdat, colcount);
      colcount += matdat.MinorDim();

      nodal_all.BlockColumnCopyAt(ndElectricDisplacement, colcount);
      colcount += ndElectricDisplacement.MinorDim();

      nodal_all.BlockColumnCopyAt(ndElectricField, colcount);
      colcount += ndElectricField.MinorDim();

      nodal_all.BlockColumnCopyAt(ndElectricScalarPotential, colcount);
      colcount += ndElectricScalarPotential.MinorDim();

      // accumulate - extrapolation done from ip's to corners => X nodes
      ElementSupport().AssembleAverage(CurrentElement().NodesX(), nodal_all);

      // element values
      if (e_codes[iCentroid]) centroid /= mass;

      // store results
      e_values.SetRow(CurrElementNumber(), element_values);

    }

    // get nodally averaged values
    const OutputSetT& output_set = ElementSupport().OutputSet(fOutputID);
    const iArrayT& nodes_used = output_set.NodesUsed();
    dArray2DT extrap_values(nodes_used.Length(), n_out);
    extrap_values.RowCollect(nodes_used, ElementSupport().OutputAverage());

    int tmpDim = extrap_values.MajorDim();
    n_values.Dimension(tmpDim, n_out);
    n_values.BlockColumnCopyAt(extrap_values, 0);
   }

} // namespace Tahoe
