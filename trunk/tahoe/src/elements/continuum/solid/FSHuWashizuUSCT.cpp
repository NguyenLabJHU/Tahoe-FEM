//
// $Id: FSHuWashizuUSCT.cpp,v 1.1 2008-07-14 16:13:25 lxmota Exp $
//
// $Log: not supported by cvs2svn $
//

#if 0

#include "FSIsotropicMatT.h"
#include "FSHuWashizuUSCT.h"
#include "FSPZMatSupportT.h"
#include "ParameterContainerT.h"
#include "ShapeFunctionT.h"

//
// materials lists (3D only)
//
#include "FSSolidMatList3DT.h"


namespace Tahoe {

  //
  //
  //
  FSHuWashizuUSCT::~FSHuWashizuUSCT()
  {
  
    delete fFSMatSupport;
    delete fShapesHW;

  }

  //
  //
  //
  void
  FSHuWashizuUSCT::Initialize()
  {

    const int nen = NumElementNodes();
    const int nsd = NumSD();
    const int nel = nen * nsd;
    const int nme = nen * nsd;

    fMaterialTangent.Dimension(nme, nme);
    fGeometricTangent.Dimension(nme, nme);

  }

  //
  // information about subordinate parameter lists
  //
  void
  FSHuWashizuUSCT::DefineSubs(SubListT& sub_list) const
  {

    //
    // inherited
    //
    FiniteStrainT::DefineSubs(sub_list);	
    
    //
    // element block/material specification
    //
    sub_list.AddSub("HW-USC_element_block", ParameterListT::OnePlus);

  }

  //
  // return the description of the given inline subordinate parameter
  // list
  //
  void
  FSHuWashizuUSCT::DefineInlineSub(const StringT& name,
					 ParameterListT::ListOrderT& order, 
					 SubListT& sub_lists) const
  {

    if (name == "HW-USC_material_choice") {
      
      order = ParameterListT::Choice;
      
      //
      // list of choices
      //
      sub_lists.AddSub("HW-USC_material_3D");

    } else { 

      //
      // inherited
      //
      FiniteStrainT::DefineInlineSub(name, order, sub_lists);

    }

  }

  //
  // return the description of the given inline subordinate parameter
  // list
  //
  ParameterInterfaceT*
  FSHuWashizuUSCT::NewSub(const StringT& name) const
  {

    ParameterInterfaceT* pip = 0;
    
    if (name == "HW-USC_material_choice") {
      
      ParameterContainerT* block = new ParameterContainerT(name);
      
      //	
      // list of element block ID's (defined by ElementBaseT)
      //
      block->AddSub("block_ID_list", ParameterListT::Once);
	
      //
      // choice of materials lists (inline)
      //
      block->AddSub("HW-USC_material_choice",
		    ParameterListT::Once, true);
	
      //
      // set this as source of subs
      //
      block->SetSubSource(this);
		
      pip = block;

    } else {

      //
      // inherited
      //
      pip = FiniteStrainT::NewSub(name);

    }

    return pip;

  }

  //
  // accept parameter list
  //
  void
  FSHuWashizuUSCT::TakeParameterList(const ParameterListT& list)
  {

    //
    // inherited
    //
    FiniteStrainT::TakeParameterList(list);
    
    //
    // allocate right Cauchy-Green deformation and 2PK stress lists
    //
  
    const int nip = NumIP();
    const int nsd = NumSD();

    fC_all.Dimension(nip*nsd);
    fC_List.Dimension(nip);

    fS_all.Dimension(nip*nsd);
    fS_List.Dimension(nip);

    for (int i = 0; i < nip; ++i) {

      fC_List[i].Set(nsd, fC_all.Pointer(i*nsd));
      fS_List[i].Set(nsd, fS_all.Pointer(i*nsd));

    }

  }

  //
  // extract the list of material parameters
  //
  void
  FSHuWashizuUSCT::CollectMaterialInfo(const ParameterListT& all_params,
					     ParameterListT& mat_params) const
  {

    const char caller[] = "FiniteStrainT::CollectMaterialInfo";
	
    //
    // initialize
    //
    mat_params.Clear();
	
    //
    // collected material parameters
    //
    const int num_blocks = all_params.NumLists("HW-USC_element_block");

    for (int i = 0; i < num_blocks; ++i) {

      //
      // block information
      //
      const ParameterListT& block =
	all_params.GetList("HW-USC_element_block", i);

      if (i == 0) {

	//
	// resolve material list name
	//
	const ParameterListT& mat_list_params =
	  block.GetListChoice(*this, "HW-USC_material_choice");
	
	mat_params.SetName(mat_list_params.Name());
	
      }

      //	
      // collect material parameters
      //
      const ParameterListT& mat_list = block.GetList(mat_params.Name());
      const ArrayT<ParameterListT>& mat = mat_list.Lists();
      mat_params.AddList(mat[0]);

    }

  }


  //
  // Protected
  //

  //
  // construct a new material support and return a pointer
  //
  MaterialSupportT*
  FSHuWashizuUSCT::NewMaterialSupport(MaterialSupportT* p) const
  {

    //
    // allocate
    //
    if (!p) p = new FSPZMatSupportT(NumDOF(), NumIP());

    //
    // inherited initializations
    //
    FiniteStrainT::NewMaterialSupport(p);
	
    return p;

  }

  //
  // construct materials manager and read data
  //
  MaterialListT*
  FSHuWashizuUSCT::NewMaterialList(const StringT& name, int size)
  {

    //
    // determine number of spatial dimensions
    //
    const int nsd = NumSD();

    MaterialListT* mlp = 0;

    if (size > 0) {

      //
      // material support
      //
      if (0 == fFSMatSupport) {

	fFSMatSupport = dynamic_cast<FSMatSupportT*>(NewMaterialSupport());

	if (0 == fFSMatSupport) {

	  ExceptionT::GeneralFail("FSHuWashizuUSCT::NewMaterialList");

	}

      }

      mlp = new FSSolidMatList3DT(size, *fFSMatSupport);

    } else {

      mlp = new FSSolidMatList3DT;

    }
	
    return mlp;

  }

  //
  // form shape functions and derivatives
  //
  void
  FSHuWashizuUSCT::SetGlobalShape()
  {

    //
    // inherited
    //
    FiniteStrainT::SetGlobalShape();

    //
    // what needs to be computed
    //
    int material_number = CurrentElement().MaterialNumber();

    for (int i = 0; i < NumIP(); i++) {

      dArrayT& C = fC_List[i];

      fShapes->InterpolateU(fLocC, C);

      dArrayT& S = fS_List[i];

      fShapes->InterpolateU(fLocS, S);

    }
    
  }

  //
  // write all current element information to the stream
  //
  void
  FSHuWashizuUSCT::CurrElementInfo(ostream& out) const
  {

    //	
    // inherited
    //
    FiniteStrainT::CurrElementInfo(out);
	
    out << std::endl << "Mixed Cauchy-Green deformation at IP:" << std::endl;

    for (int i = 0; i < fC_List.Length(); ++i) {

      out << " ip: " << i+1 << std::endl << fC_List[i] << std::endl;

    }

    out << std::endl;

    out << std::endl << "Mixed 2PK stress at IP:" << std::endl;

    for (int i = 0; i < fS_List.Length(); ++i) {

      out << " ip: " << i+1 << std::endl << fS_List[i] << std::endl;

    }

    out << std::endl;

  }

  //
  // Initialize local arrays
  //
  void
  FSHuWashizuUSCT::SetLocalArrays()
  {

    const FieldT* Cbar =
      ElementSupport().Field("mixed Cauchy-Green deformation");

    if (0 == Cbar) {

      std::cout << std::endl;
      std::cout << "FSHuWashizuUSCT::SetLocalArrays: ";
      std::cout << "mixed Cauchy-Green deformation field not found.";
      std::cout << std::endl;

      throw ExceptionT::kGeneralFail;

    }

    const FieldT* Sbar =
      ElementSupport().Field("mixed 2PK stress");

    if (0 == Sbar) {

      std::cout << std::endl;
      std::cout << "FSHuWashizuUSCT::SetLocalArrays: ";
      std::cout << "mixed 2PK stress field not found.";
      std::cout << std::endl;

      throw ExceptionT::kGeneralFail;

    }

    //
    // Inherited
    //
    FiniteStrainT::SetLocalArrays();

    //
    // Allocate storage space
    //
    const int nen  = NumElementNodes();
    const int ndof = dSymMatrixT::NumValues( NumDOF() );

    fLocC.Dimension(nen, ndof);
    fLocS.Dimension(nen, ndof);

    //
    // Register fields
    //
    Cbar->RegisterLocal(fLocC);
    Sbar->RegisterLocal(fLocS);

  }

  //
  // Strain displacement operator
  //
  void
  FSHuWashizuUSCT::Set_B_C(const dArray2DT& DNaX, dMatrixT& B_C)
  {


    const int nsd = NumSD();
    const int nen = NumElementNodes();
    const int StrainDim = dSymMatrixT::NumValues(nsd);
    const int nnd = DNaX.MinorDim();

    assert(nen == nnd);

    B_C.Dimension(StrainDim, nsd*nnd);


    const int ij2i[] = {0,1,2,1,0,0};
    const int ij2j[] = {0,1,2,2,2,1};

    const dMatrixT F = DeformationGradient();


    for (int ij = 0; ij < StrainDim; ++ij) {
      
      const int i = ij2i[ij];
      const int j = ij2j[ij];
      
      for (int a = 0; a < nnd; ++a) {
	
	for (int k = 0; k < nsd; ++k) {
	  
	  int ka = a * nsd + k;
	  
	  B_C(ij, ka) = 0.5 * (DNaX(i, a) * F(k, j) + DNaX(j, a) * F(k, i));
	  
	  // Shear components doubled to conform with Voigt
	  // convention
	  if (ij >= nsd) {
	  
	    B_C(ij, ka) *= 2.0;
	    
	  }
	  
	}
	
      }
      
    }
    
  }

  //
  // Electric displacement - vector potential operator
  //
  void
  FSHuWashizuUSCT::Set_B_D(const dArray2DT& DNaX, dMatrixT& B_D)
  {

    //
    // Set the vector potential -> electric displacement operator as
    // a matrix instead of a vector to avoid using the cross
    // product.  Regular matrix - vector products are preferred.
    //
    // D = - Curl psi
    //

    const int nsd = NumSD();
    const int nen = NumElementNodes();
    const int nnd = DNaX.MinorDim();

    assert(nen == nnd);

    B_D.Dimension(nsd, nsd*nnd);

    for (int a = 0; a < nnd; ++a) {

      const int ka = nsd * a;
      
      B_D(0, ka + 0) = 0.0;
      B_D(0, ka + 1) = DNaX(2, a);
      B_D(0, ka + 2) = - DNaX(1, a);
      
      B_D(1, ka + 0) = - DNaX(2, a);
      B_D(1, ka + 1) = 0.0;
      B_D(1, ka + 2) = DNaX(0, a);
      
      B_D(2, ka + 0) = DNaX(1, a);
      B_D(2, ka + 1) = - DNaX(0, a);
      B_D(2, ka + 2) = 0.0;
      
    }
    
  }
  

  //
  // Geometric stiffness
  //
  void
  FSHuWashizuUSCT::AccumulateGeometricStiffness(dMatrixT& Kg,
						      const dArray2DT& DNaX,
						      dSymMatrixT& S)
  {

    const int nsd = NumSD();
    const int nen = NumElementNodes();
    const int nnd = DNaX.MinorDim();

    assert(nen == nnd);

    dMatrixT GradNa;
    GradNa.Dimension(nsd, nnd);
    fShapes->GradNa(DNaX, GradNa);

    dMatrixT DNaS;
    DNaS.Dimension(nnd, nnd);

    dMatrixT Sm;
    S.ToMatrix(Sm);

    DNaS.MultQTBQ(GradNa, Sm);

    dMatrixT SI;
    SI.Dimension(nsd);

    for (int i = 0; i < nen; ++i) {

      for (int j = 0; j < nen; ++j) {

	SI.Identity(DNaS(i,j));

	Kg.AddBlock(nsd*i, nsd*j, SI);

      }

    }

  }


  //
  // strain-displacement operator
  //
  void
  FSHuWashizuUSCT::Set_B(const dArray2DT& DNaX, dMatrixT& B)
  {

    const int nsd = NumSD();
    const int StrainDim = dSymMatrixT::NumValues(nsd);
    const int nen = NumElementNodes();

    B.Dimension(StrainDim + nsd, nsd * nen);

    dMatrixT B_C;

    Set_B_C(DNaX, B_C);

    B.SetBlock(0, 0, B_C);

    dMatrixT B_D;

    Set_B_D(DNaX, B_D);

    B.SetBlock(StrainDim, 0, B_C);

  }


  //
  // increment current element
  //
  bool
  FSHuWashizuUSCT::NextElement()
  {

    bool isThereNext = FiniteStrainT::NextElement();

    if (isThereNext == true) {

      const int index = CurrentElement().MaterialNumber();

      ContinuumMaterialT* pMaterial = (*fMaterialList)[index];

      fCurrMaterial = dynamic_cast<FSNeoHookePZLinT*>(pMaterial);

    }

    return isThereNext;

  }


  //
  // element stiffness matrix
  //
  void
  FSHuWashizuUSCT::FormStiffness(double constK)
  {

    //
    // matrix format
    //
    dMatrixT::SymmetryFlagT format =
      (fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
      dMatrixT::kWhole :
      dMatrixT::kUpperOnly;

    //
    // integration
    //
    const double* pDet    = fShapes->IPDets();
    const double* pWeight = fShapes->IPWeights();

    //
    // integrate over element
    //
    const int nsd = NumSD();
    const int nen = NumElementNodes();

    fGeometricTangent = 0.0;
    fMaterialTangent = 0.0;
    fMechanical2ElectricTangent = 0.0;
    fElectric2MechanicalTangent = 0.0;
    fElectricTangent = 0.0;

    fShapes->TopIP();

    while (fShapes->NextIP()) {

      //
      // scale/weighting factor for integration
      //
      const double sw = constK * (*pDet++) * (*pWeight);

      //
      // get material tangent moduli
      //
      dMatrixT C;
      C.SetToScaled(sw, fCurrMaterial->C_ijkl());

      dMatrixT H;
      H.SetToScaled(sw, fCurrMaterial->H_ijk());

      dMatrixT B;
      B.SetToScaled(sw, fCurrMaterial->B_ij());

      dSymMatrixT S;
      S.SetToScaled(sw, fCurrMaterial->S_ij());

      //
      // prepare derivatives of interpolation functions
      //
      const dArray2DT & DNa = fShapes->Derivatives_X();

      //
      // Primary variables -> gradient operators
      //

      dMatrixT B_C;

      Set_B_C(DNa, B_C);

      dMatrixT B_D;

      Set_B_D(DNa, B_D);

      //
      // Material stiffness
      //
      fMaterialTangent.MultQTBQ(B_C, C, format, dMatrixT::kAccumulate);

      //
      // Geometric stiffness
      //
      AccumulateGeometricStiffness(fGeometricTangent, DNa, S);

      //
      // HW-USC coupling tangents
      //
      fMechanical2ElectricTangent.MultATBC(B_D, H, B_C, dMatrixT::kWhole,
					   dMatrixT::kAccumulate);

      //
      // Purely electric tangent
      //
      fElectricTangent.MultQTBQ(B_D, B, format, dMatrixT::kAccumulate);

    }

    fElectric2MechanicalTangent.Transpose(fMechanical2ElectricTangent);

    //
    // Add geometric stiffness
    //
    fMaterialTangent.Expand(fGeometricTangent, 2*NumSD(),
			    dMatrixT::kAccumulate);

    //
    // Assemble into element stiffness matrix
    //
    fLHS.AddBlock(0, 0,
		  fMaterialTangent);

    fLHS.AddBlock(fMaterialTangent.Rows(), fMaterialTangent.Cols(),
		  fElectricTangent);

    fLHS.AddBlock(0, fMaterialTangent.Cols(),
		  fElectric2MechanicalTangent);

    //
    // non-symmetric
    //
    if (format != dMatrixT::kUpperOnly) {

      fLHS.AddBlock(fMaterialTangent.Rows(), 0,
		    fMechanical2ElectricTangent);

    }

  }

  //
  // internal force
  //
  void
  FSHuWashizuUSCT::FormKd(double constK)
  {

    //
    //
    //
    const int neq = NumElementNodes()*NumDOF();
    const int nsd = NumSD();
    dArrayT RHS(neq, fRHS.Pointer());

    //
    // integration scheme
    //
    const double* Det    = fShapes->IPDets();
    const double* Weight = fShapes->IPWeights();

    fShapes->TopIP();

    while ( fShapes->NextIP() ){

      //
      // integration weight
      //
      const double w = constK * (*Weight++) * (*Det++);

      dSymMatrixT S = fCurrMaterial->S_ij();

      dMatrixT s;
      s.Dimension(dSymMatrixT::NumValues(nsd), 1);
      s(0,0)=S(0,0);
      s(1,0)=S(1,1);
      s(2,0)=S(2,2);
      s(3,0)=S(2,1);
      s(4,0)=S(2,0);
      s(5,0)=S(1,0);

      const dArray2DT & DNa = fShapes->Derivatives_X();

      dMatrixT B_D;

      Set_B_D(DNa, B_D);

      B_D.Multx(s, RHS, dMatrixT::kAccumulate);

    }

  }

} // namespace Tahoe

#endif
