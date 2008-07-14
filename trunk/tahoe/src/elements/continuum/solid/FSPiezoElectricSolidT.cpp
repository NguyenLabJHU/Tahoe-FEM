//
// $Id: FSPiezoElectricSolidT.cpp,v 1.2 2008-07-14 17:37:23 lxmota Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.1  2008/06/16 18:15:51  lxmota
// Piezoelectric solid. Initial source.
//
//

#include "FSNeoHookePZLinT.h"
#include "FSPiezoElectricSolidT.h"
#include "FSPZMatSupportT.h"
#include "ParameterContainerT.h"
#include "ShapeFunctionT.h"
#include "XDOF_ManagerT.h"

//
// materials lists (3D only)
//
#include "FSSolidMatList3DT.h"


namespace Tahoe {

  //
  //
  //
  FSPiezoElectricSolidT::~FSPiezoElectricSolidT()
  {

    delete fFSPZMatSupport;

  }

  //
  //
  //
  void
  FSPiezoElectricSolidT::Initialize()
  {

    const int nen = NumElementNodes();
    const int nsd = NumSD();
    const int nel = nen * nsd;
    const int nme = nen * nsd;
    const int dof = ManifoldDim() + ElectricalDim();
    const int neq = nen * dof;

    fMaterialTangent.Dimension(nme, nme);
    fGeometricTangent.Dimension(nme, nme);
    fMechanical2ElectricTangent.Dimension(nel, nme);
    fElectric2MechanicalTangent.Dimension(nme, nel);
    fElectricTangent.Dimension(nel, nel);

    fLHS.Dimension(neq);
    fRHS.Dimension(neq);


  }

  //
  // specify parameters needed by the interface
  //
  void
  FSPiezoElectricSolidT::DefineParameters(ParameterListT& list) const
  {
    // inherited
    FiniteStrainT::DefineParameters(list);

    // additional fields
    list.AddParameter(ParameterT::Word, "electric_field_name");

  }

#if 0
  //
  // information about subordinate parameter lists
  //
  void
  FSPiezoElectricSolidT::DefineSubs(SubListT& sub_list) const
  {

    //
    // inherited
    //
    FiniteStrainT::DefineSubs(sub_list);

    //
    // element block/material specification
    //
    sub_list.AddSub("piezoelectric_element_block", ParameterListT::OnePlus);

  }

  //
  // return the description of the given inline subordinate parameter
  // list
  //
  void
  FSPiezoElectricSolidT::DefineInlineSub(const StringT& name,
      ParameterListT::ListOrderT& order,
      SubListT& sub_lists) const
      {

    if (name == "piezoelectric_material_choice") {

      order = ParameterListT::Choice;

      //
      // list of choices
      //
      sub_lists.AddSub("piezoelectric_material_3D");

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
  FSPiezoElectricSolidT::NewSub(const StringT& name) const
  {

    ParameterInterfaceT* pip = 0;

    if (name == "piezoelectric_material_choice") {

      ParameterContainerT* block = new ParameterContainerT(name);

      //
      // list of element block ID's (defined by ElementBaseT)
      //
      block->AddSub("block_ID_list", ParameterListT::Once);

      //
      // choice of materials lists (inline)
      //
      block->AddSub("piezoelectric_material_choice",
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
  // extract the list of material parameters
  //
  void
  FSPiezoElectricSolidT::CollectMaterialInfo(const ParameterListT& all_params,
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
    const int num_blocks = all_params.NumLists("piezoelectric_element_block");

    for (int i = 0; i < num_blocks; ++i) {

      //
      // block information
      //
      const ParameterListT& block =
        all_params.GetList("piezoelectric_element_block", i);

      if (i == 0) {

        //
        // resolve material list name
        //
        const ParameterListT& mat_list_params =
          block.GetListChoice(*this, "piezoelectric_material_choice");

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

#endif

  //
  // accept parameter list
  //
  void
  FSPiezoElectricSolidT::TakeParameterList(const ParameterListT& list)
  {

    //
    // inherited
    //
    FiniteStrainT::TakeParameterList(list);

    //
    // get electric vector potential field
    //
    // for now use same integration and interpolation schemes as primary field
    //
    const StringT& electric_field_name =
      list.GetParameter("electric_field_name");

    fElectricVectorPotentialField = ElementSupport().Field(electric_field_name);
    if (!fElectricVectorPotentialField) {
      ExceptionT::GeneralFail("FSPiezoElectricSolidT::TakeParameterList",
          "could not resolve \"%s\" field",
          electric_field_name.Pointer());
    }

#if 0
    //
    // offset to class needs flags
    //
    fNeedsOffset = fMaterialNeeds[0].Length();

    //
    // set material needs
    //
    for (int i = 0; i < fMaterialNeeds.Length(); ++i) {

      //
      // array of needs
      //
      ArrayT<bool>& needs = fMaterialNeeds[i];

      needs.Resize(needs.Length() + 2, true);

      //
      // casts are safe since class contructs materials list
      //
      ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
      FSNeoHookePZLinT* mat = dynamic_cast<FSNeoHookePZLinT*>(pcont_mat);

      //
      // collect needs
      //
      needs[fNeedsOffset + kD     ] = mat->Need_D();
      needs[fNeedsOffset + kD_last] = mat->Need_D_last();

      //
      // consistency
      //
      needs[kNeedVectorPotential] =
        needs[kNeedVectorPotential] || needs[fNeedsOffset + kD];

      needs[kNeedLastVectorPotential] =
        needs[kNeedLastVectorPotential] || needs[fNeedsOffset + kD_last];

    }

    //
    // what's needed
    //
    bool need_D      = false;
    bool need_D_last = false;

    for (int i = 0; i < fMaterialList->Length(); ++i) {

      need_D      = need_D || Needs_D(i);
      need_D_last = need_D_last || Needs_D_last(i);

    }


    //
    // allocate electric displacement list
    //
    if (need_D) {

      const int nip = NumIP();
      const int nsd = NumSD();

      fD_all.Dimension(nip*nsd);
      fD_List.Dimension(nip);

      for (int i = 0; i < nip; ++i) {

        fD_List[i].Set(nsd, fD_all.Pointer(i*nsd));

      }

    }

    //
    // allocate "last" electric displacement list
    //
    if (need_D_last) {

      const int nip = NumIP();
      const int nsd = NumSD();

      fD_last_all.Dimension(nip*nsd);
      fD_last_List.Dimension(nip);

      for (int i = 0; i < nip; ++i) {

        fD_last_List[i].Set(nsd, fD_last_all.Pointer(i*nsd));

      }

    }
#endif

    const int nip = NumIP();
    const int nsd = NumSD();

    fD_all.Dimension(nip*nsd);
    fD_List.Dimension(nip);

    for (int i = 0; i < nip; ++i) {

      fD_List[i].Set(nsd, fD_all.Pointer(i*nsd));

    }

    fD_last_all.Dimension(nip*nsd);
    fD_last_List.Dimension(nip);

    for (int i = 0; i < nip; ++i) {

      fD_last_List[i].Set(nsd, fD_last_all.Pointer(i*nsd));

    }

    // Allocate workspace
    Initialize();
  }

#if 0
  //
  //
  //
  void FSPiezoElectricSolidT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
    AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
  {

    // loop over connectivity blocks
    int element_count = 0;
    for (int i = 0; i < fConnectivities.Length(); i++)
    {
      // block connectivities and equations
      const iArray2DT& connects = *fConnectivities[i];
      iArray2DT& eqnos = fEqnos[i];

      // resize equations array
      const int neq = NumElementNodes() * ( ManifoldDim() + ElectricalDim() );
      eqnos.Dimension( connects.MajorDim(), neq );

      // set displacement equations
      Field().SetLocalEqnos(connects, eqnos);

      // xdof equations for the block
      const iArray2DT xdof_eqnos(connects.MajorDim(),
          xdof_eqnos_all.MinorDim(), xdof_eqnos_all(element_count));
      element_count += connects.MajorDim();

      // fill columns with xdof equations
      iArrayT tmp(xdof_eqnos.MajorDim());
      for (int i = neq; i < eqnos.MinorDim(); i++)
      {
        xdof_eqnos.ColumnCopy(i - neq, tmp);
        eqnos.SetColumn(i, tmp);
      }

      // add to list
      eq_1.AppendUnique(&eqnos);
    }

    // reset element cards
    SetElementCards(fBlockData, fConnectivities, fEqnos, fElementCards);
  }
#endif

  //
  // Protected
  //

  //
  // construct a new material support and return a pointer
  //
  MaterialSupportT*
  FSPiezoElectricSolidT::NewMaterialSupport(MaterialSupportT* p) const
  {

    //
    // allocate
    //
    if (!p) p = new FSPZMatSupportT(NumDOF(), NumIP());

    //
    // inherited initializations
    //
    FiniteStrainT::NewMaterialSupport(p);

    //
    // set parent class fields
    //
    FSPZMatSupportT* ps = dynamic_cast<FSPZMatSupportT*>(p);

    if (ps != 0) {

      ps->SetElectricDisplacement(&fD_List);
      ps->SetElectricDisplacement_last(&fD_last_List);

    }

    return p;

  }

  //
  // construct materials manager and read data
  //
  MaterialListT*
  FSPiezoElectricSolidT::NewMaterialList(const StringT& name, int size)
  {

    if (name != "large_strain_material_3D") {
      return 0;
    }

    MaterialListT* mlp = 0;

    if (size > 0) {

      //
      // material support
      //
      if (0 == fFSPZMatSupport) {

        fFSPZMatSupport = dynamic_cast<FSPZMatSupportT*>(NewMaterialSupport());

        if (0 == fFSPZMatSupport) {

          ExceptionT::GeneralFail("FSPiezoElectricSolidT::NewMaterialList");

        }

      }

      mlp = new FSSolidMatList3DT(size, *fFSPZMatSupport);

    } else {

      mlp = new FSSolidMatList3DT;

    }

    return mlp;

  }

  //
  // form shape functions and derivatives
  //
  void
  FSPiezoElectricSolidT::SetGlobalShape()
  {

    //
    // inherited
    //
    FiniteStrainT::SetGlobalShape();

    //
    // what needs to be computed
    //
#if 0
    int material_number = CurrentElement().MaterialNumber();
    bool needs_D        = Needs_D(material_number);
    bool needs_D_last   = Needs_D_last(material_number);
#endif

    SetLocalU(fLocVectorPotential);
    SetLocalU(fLocLastVectorPotential);

    for (int i = 0; i < NumIP(); i++) {

      //
      // electric displacement
      //
      { //if (needs_D) {

        dArrayT& D = fD_List[i];

        //
        // curl of vector potential
        //
        int m = fLocVectorPotential.NumberOfNodes();
        int n = fLocVectorPotential.MinorDim();

        ArrayT<dArrayT> vp(m);

        for (int j = 0; j < m; ++j) {

          vp[j].Dimension(n);

          for (int k = 0; k < n; ++k) {

            vp[j][k] = fLocVectorPotential(j,k);

          }

        }

        fShapes->CurlU(vp, D, i);

        //
        // reverse sign
        //
        D *= -1.0;

      }

      //
      // "last" electric displacement
      //

      { //if (needs_D_last)	{


        dArrayT& D = fD_last_List[i];

        //
        // curl of vector potential
        //
        int m = fLocLastVectorPotential.NumberOfNodes();
        int n = fLocLastVectorPotential.MinorDim();

        ArrayT<dArrayT> vp(m);

        for (int j = 0; j < m; ++j) {

          vp[j].Dimension(n);

          for (int k = 0; k < n; ++k) {

            vp[j][k] = fLocLastVectorPotential(j,k);

          }

        }

        fShapes->CurlU(vp, D, i);

        //
        // reverse sign
        //
        D *= -1.0;

      }

    }

  }

  //
  // write all current element information to the stream
  //
  void
  FSPiezoElectricSolidT::CurrElementInfo(ostream& out) const
  {

    //
    // inherited
    //
    FiniteStrainT::CurrElementInfo(out);

    //
    // write deformation gradients
    //
    out << "\n electric displacement at IP:\n";

    for (int i = 0; i < fD_List.Length(); ++i) {

      out << " ip: " << i+1 << '\n' << fD_List[i] << '\n';

    }

    out << '\n';

  }

  //
  // Initialize local arrays
  //
  void
  FSPiezoElectricSolidT::SetLocalArrays()
  {

    //
    // look for an electric vector potential field
    //
    const FieldT* evp = 0;

    if (0 == fElectricVectorPotentialField) {
      evp = ElementSupport().Field("electric_vector_potential");
      fElectricVectorPotentialField = evp;
    } else {
      evp = fElectricVectorPotentialField;
    }

    if (0 == evp) {

      std::cout << std::endl;
      std::cout << "FSPiezoElectricSolidT::SetLocalArrays: ";
      std::cout << "Electric vector potential field not found.";
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
    const int ndof = NumDOF();

    fLocLastVectorPotential.Dimension(nen, ndof);
    fLocVectorPotential.Dimension(nen, ndof);

    //
    // Register fields
    //
    evp->RegisterLocal(fLocLastVectorPotential);
    evp->RegisterLocal(fLocVectorPotential);

  }

  //
  // Strain displacement operator
  //
  void
  FSPiezoElectricSolidT::Set_B_C(const dArray2DT& DNaX, dMatrixT& B_C)
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
  FSPiezoElectricSolidT::Set_B_D(const dArray2DT& DNaX, dMatrixT& B_D)
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
  FSPiezoElectricSolidT::AccumulateGeometricStiffness(dMatrixT& Kg,
      const dArray2DT& DNaX,
      dSymMatrixT& S)
  {

    const int nsd = NumSD();
    const int nen = NumElementNodes();
    const int nnd = DNaX.MinorDim();

    assert(nen == nnd);

    dMatrixT GradNa(nsd, nnd);
    fShapes->GradNa(DNaX, GradNa);

    dMatrixT DNaS(nnd, nnd);

    dMatrixT Sm(ManifoldDim());
    S.ToMatrix(Sm);

    DNaS.MultQTBQ(GradNa, Sm);

    dMatrixT SI(nsd);

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
  FSPiezoElectricSolidT::Set_B(const dArray2DT& DNaX, dMatrixT& B)
  {

    const int nsd = NumSD();
    const int nen = NumElementNodes();

    B.Dimension(StrainDim() + ElectricalDim(), nsd * nen);

    dMatrixT B_C;

    Set_B_C(DNaX, B_C);

    B.SetBlock(0, 0, B_C);

    dMatrixT B_D;

    Set_B_D(DNaX, B_D);

    B.SetBlock(StrainDim(), 0, B_C);

  }


  //
  // increment current element
  //
  bool
  FSPiezoElectricSolidT::NextElement()
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
  FSPiezoElectricSolidT::FormStiffness(double constK)
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
      dMatrixT C(StrainDim());
      C.SetToScaled(sw, fCurrMaterial->C_ijkl());

      dMatrixT H(ElectricalDim(), StrainDim());
      H.SetToScaled(sw, fCurrMaterial->H_ijk());

      dMatrixT B(ElectricalDim());
      B.SetToScaled(sw, fCurrMaterial->B_ij());

      dSymMatrixT S(ManifoldDim());
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
      // Piezoelectric coupling tangents
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
    fMaterialTangent.Expand(fGeometricTangent, 1, dMatrixT::kAccumulate);

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
  FSPiezoElectricSolidT::FormKd(double constK)
  {

    //
    //
    //
    const int neq = NumElementNodes()*(ManifoldDim() + ElectricalDim());
    const int nsd = NumSD();

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
      dArrayT E     = fCurrMaterial->E_i();

      dArrayT s(dSymMatrixT::NumValues(nsd));
      s[0]=S(0,0);
      s[1]=S(1,1);
      s[2]=S(2,2);
      s[3]=S(2,1);
      s[4]=S(2,0);
      s[5]=S(1,0);

      const dArray2DT & DNa = fShapes->Derivatives_X();

      //
      // AQUI: Residuo necesita adicion de campo electrico.
      //

      dMatrixT B_C;

      Set_B_C(DNa, B_C);

      B_C.MultTx(s, fRHS, dMatrixT::kAccumulate);

    }

  }

} // namespace Tahoe
