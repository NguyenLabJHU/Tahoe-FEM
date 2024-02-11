//
// $Id: FSPiezoElectricSolidT.cpp,v 1.1 2009-05-05 23:34:55 beichuan Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.2  2009/04/02 00:55:37  amota
// Changed output to show divergence of electric vector potential.
//
// Revision 1.1  2008/12/12 18:56:28  amota
// Moved to new localtion.
//
// Revision 1.2  2008/07/14 17:37:23  lxmota
// Various corrections related to initialization.
//
// Revision 1.1  2008/06/16 18:15:51  lxmota
// Piezoelectric solid. Initial source.
//
//

#include "FSNeoHookePZLinT.h"
#include "FSPiezoElectricSolidT.h"
#include "FSPZMatSupportT.h"
#include "ParameterContainerT.h"
#include "OutputSetT.h"
#include "ShapeFunctionT.h"

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

    if (0 != fFSPZMatSupport) delete fFSPZMatSupport;

  }

  //
  // specify parameters needed by the interface
  //
  void FSPiezoElectricSolidT::DefineParameters(ParameterListT& list) const
  {
    // inherited
    FiniteStrainT::DefineParameters(list);

    // additional fields
    list.AddParameter(ParameterT::Word, "electric_field_name");

  }

  //
  // accept parameter list
  //
  void FSPiezoElectricSolidT::TakeParameterList(const ParameterListT& list)
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
    const StringT& electric_field_name = list.GetParameter(
        "electric_field_name");

    fElectricVectorPotentialField = ElementSupport().Field(electric_field_name);
    if (!fElectricVectorPotentialField) {
      ExceptionT::GeneralFail("FSPiezoElectricSolidT::TakeParameterList",
          "could not resolve \"%s\" field", electric_field_name.Pointer());
    }

    const int nip = NumIP();
    const int nsd = NumSD();

    fD_all.Dimension(nip * nsd);
    fD_List.Dimension(nip);

    for (int i = 0; i < nip; ++i) {
      fD_List[i].Alias(nsd, fD_all.Pointer(i * nsd));
    }

    fDivPhi_all.Dimension(nip);
    fDivPhi_List.Dimension(nip);
    fDivPhi_List.Alias(nip, fDivPhi_all.Pointer());

    // Allocate workspace
    Workspace();

  }

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
    FSPZMatSupportT* ps = dynamic_cast<FSPZMatSupportT*> (p);

    if (ps != 0) {

      ps->SetElectricDisplacement(&fD_List);
      ps->SetDivergenceVectorPotential(&fDivPhi_List);

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
        fFSPZMatSupport = dynamic_cast<FSPZMatSupportT*> (NewMaterialSupport());

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
  void FSPiezoElectricSolidT::SetGlobalShape()
  {

    //
    // inherited
    //
    FiniteStrainT::SetGlobalShape();

    //
    // what needs to be computed
    //
    SetLocalU(fLocVectorPotential);

    for (int i = 0; i < NumIP(); i++) {

      //
      // electric displacement and divergence of vector potential
      //
      {

        dArrayT& D = fD_List[i];
        double& DivPhi = fDivPhi_List[i];

        //
        // curl of vector potential
        //
        int m = fLocVectorPotential.NumberOfNodes();
        int n = fLocVectorPotential.MinorDim();

        ArrayT<dArrayT> vp(m);
        dMatrixT GradPhi(n);

        for (int j = 0; j < m; ++j) {
          vp[j].Dimension(n);

          for (int k = 0; k < n; ++k) {
            vp[j][k] = fLocVectorPotential(j, k);
          }
        }

        fShapes->CurlU(vp, D, i);
        fShapes->GradU(fLocVectorPotential, GradPhi, i);
        DivPhi = GradPhi.Trace();

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
  void FSPiezoElectricSolidT::CurrElementInfo(ostream& out) const
  {

    //
    // inherited
    //
    FiniteStrainT::CurrElementInfo(out);

    //
    // write deformation gradients
    //
    out << std::endl;
    out << "electric displacement";
    out << " and divergence of vector potential at IP:";
    out << std::endl;

    for (int i = 0; i < fD_List.Length(); ++i) {
      out << " ip: " << i + 1 << std::endl << fD_List[i] << std::endl;
      out << fDivPhi_List[i] << std::endl;
    }

    out << std::endl;

  }

  //
  // Initialize local arrays
  //
  void FSPiezoElectricSolidT::SetLocalArrays()
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
    const int nen = NumElementNodes();
    const int ndof = NumDOF();

    fLocVectorPotential.Dimension(nen, ndof);

    //
    // Register fields
    //
    evp->RegisterLocal(fLocVectorPotential);

  }

  //
  // Strain displacement operator
  //
  void FSPiezoElectricSolidT::Set_B_C(const dArray2DT& DNaX, dMatrixT& B_C)
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
  // Electric displacement - vector potential operator
  //
  void FSPiezoElectricSolidT::Set_B_D(const dArray2DT& DNaX, dMatrixT& B_D)
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

    B_D.Dimension(nsd, nsd * nen);

    for (int a = 0; a < nen; ++a) {

      const int ka = nsd * a;

      B_D(0, ka + 0) = 0.0;
      B_D(0, ka + 1) = DNaX(2, a);
      B_D(0, ka + 2) = -DNaX(1, a);

      B_D(1, ka + 0) = -DNaX(2, a);
      B_D(1, ka + 1) = 0.0;
      B_D(1, ka + 2) = DNaX(0, a);

      B_D(2, ka + 0) = DNaX(1, a);
      B_D(2, ka + 1) = -DNaX(0, a);
      B_D(2, ka + 2) = 0.0;

    }

  }

  //
  // Electric displacement - vector potential operator, penalty term
  //
  void FSPiezoElectricSolidT::Set_B_K(const dArray2DT& DNaX, dMatrixT& B_K)
  {

    //
    // Div Phi = Grad Na Phi_a
    //

    const int nsd = NumSD();
    const int nen = NumElementNodes();
    const int nnd = DNaX.MinorDim();

    assert(nen == nnd);

    B_K.Dimension(1, nsd * nen);

    for (int i = 0; i < nen; ++i) {
      for (int j = 0; j < nsd; ++j) {
        B_K(0, i * nsd + j) = DNaX(j, i);
      }
    }

  }

  //
  // Geometric stiffness
  //
  void FSPiezoElectricSolidT::AccumulateGeometricStiffness(dMatrixT& Kg,
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

  //
  // strain-displacement operator
  //
  void FSPiezoElectricSolidT::Set_B(const dArray2DT& DNaX, dMatrixT& B)
  {

    const int nen = NumElementNodes();

    B.Dimension(StrainDim() + ElectricalDim(), TotalNumDOF() * nen);

    dMatrixT B_C;

    Set_B_C(DNaX, B_C);

    B.SetBlock(0, 0, B_C);

    dMatrixT B_D;

    Set_B_D(DNaX, B_D);

    B.SetBlock(StrainDim(), NumDOF() * nen, B_D);

  }

  //
  // increment current element
  //
  bool FSPiezoElectricSolidT::NextElement()
  {

    bool isThereNext = FiniteStrainT::NextElement();

    if (isThereNext == true) {

      const int index = CurrentElement().MaterialNumber();

      ContinuumMaterialT* pMaterial = (*fMaterialList)[index];

      fCurrMaterial = dynamic_cast<FSNeoHookePZLinT*> (pMaterial);

    }

    return isThereNext;

  }

  //
  //
  //
  void FSPiezoElectricSolidT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
      AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
  {

    for (int i = 0; i < fEqnos.Length(); ++i) {

      const int ndf = NumDOF();
      const int nen = fConnectivities[i]->MinorDim();
      const int offset = ndf * nen;

      fElectricVectorPotentialField->SetLocalEqnos(*fConnectivities[i],
          fEqnos[i], offset);

    }

    ElementBaseT::Equations(eq_1, eq_2);

  }

  //
  // element stiffness matrix
  //
  void FSPiezoElectricSolidT::FormStiffness(double constK)
  {

    //
    // matrix format
    //
    dMatrixT::SymmetryFlagT format = (fLHS.Format()
        == ElementMatrixT::kNonSymmetric)
        ? dMatrixT::kWhole
        : dMatrixT::kUpperOnly;

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

    while (fShapes->NextIP() != 0) {

      //
      // scale/weighting factor for integration
      //
      const double w = constK * fShapes->IPDet() * fShapes->IPWeight();

      //
      // get material tangent moduli
      //
      dMatrixT C = fCurrMaterial->C_IJKL();
      dMatrixT H = fCurrMaterial->H_IJK();
      dMatrixT B = fCurrMaterial->B_IJ();
      dSymMatrixT S = fCurrMaterial->S_IJ();

      C *= (0.25 * w);
      H *= (0.5 * w);
      B *= w;
      S *= w;

      //
      // prepare derivatives of interpolation functions
      //
      const dArray2DT & DNaX = fShapes->Derivatives_X();

      //
      // Primary variables -> gradient operators
      //

      dMatrixT B_C;

      Set_B_C(DNaX, B_C);

      dMatrixT B_D;

      Set_B_D(DNaX, B_D);

      //
      // Material stiffness
      //
      fMaterialTangent.MultQTBQ(B_C, C, format, dMatrixT::kAccumulate);

      //
      // Geometric stiffness
      //
      AccumulateGeometricStiffness(fGeometricTangent, DNaX, S);

      //
      // Piezoelectric coupling tangents
      //
      fMechanical2ElectricTangent.MultATBC(B_D, H, B_C, dMatrixT::kWhole,
          dMatrixT::kAccumulate);

      //
      // Purely electric tangent
      //
      fElectricTangent.MultQTBQ(B_D, B, format, dMatrixT::kAccumulate);

      //
      // Penalty contribution
      //
      {
        dMatrixT TwoKay(1);
        TwoKay.Identity(w * 2.0 * fCurrMaterial->GetPenaltyCoefficient());

        dMatrixT B_K;
        Set_B_K(DNaX, B_K);

        fElectricTangent.MultQTBQ(B_K, TwoKay, format, dMatrixT::kAccumulate);
      }

    }

    fElectric2MechanicalTangent.Transpose(fMechanical2ElectricTangent);

    //
    // Add geometric stiffness
    //
    fMaterialTangent.Expand(fGeometricTangent, 1, dMatrixT::kAccumulate);

    //
    // Assemble into element stiffness matrix
    //
    fLHS.AddBlock(0, 0, fMaterialTangent);

    fLHS.AddBlock(fMaterialTangent.Rows(), fMaterialTangent.Cols(),
        fElectricTangent);

    fLHS.AddBlock(0, fMaterialTangent.Cols(), fElectric2MechanicalTangent);

    //
    // non-symmetric
    //
    if (format != dMatrixT::kUpperOnly) {

      fLHS.AddBlock(fMaterialTangent.Rows(), 0, fMechanical2ElectricTangent);

    }

  }

  //
  // internal force
  //
  void FSPiezoElectricSolidT::FormKd(double constK)
  {

    //
    //
    //
    const int nsd = NumSD();
    const int nen = NumElementNodes();

    //
    // integration scheme
    //
    dArrayT RS(nsd * nen);
    RS = 0.0;
    dArrayT RE(nsd * nen);
    RE = 0.0;
    dArrayT R(RS.Length() + RE.Length());
    R = 0.0;
    dArrayT R_K(2 * nsd * nen);
    R_K = 0.0;

    fShapes->TopIP();

    while (fShapes->NextIP() != 0) {

      //
      // integration weight
      //
      const double w = constK * fShapes->IPDet() * fShapes->IPWeight();

      //
      // Deformation power P_D := \int S : 0.5 \dot(C) dV
      //
      dSymMatrixT S = fCurrMaterial->S_IJ();
      dArrayT E = fCurrMaterial->E_I();

      S *= (0.5 * w);
      E *= w;

      const dArray2DT & DNaX = fShapes->Derivatives_X();

      dMatrixT B_C;
      Set_B_C(DNaX, B_C);
      B_C.MultTx(S, RS, 1.0, dMatrixT::kAccumulate);

      dMatrixT B_D;
      Set_B_D(DNaX, B_D);
      B_D.MultTx(E, RE, 1.0, dMatrixT::kAccumulate);

      //
      // Penalty contribution
      //
      {
        double DivPhi2K = fCurrMaterial->DivPhi();
        DivPhi2K *= (2.0 * w * fCurrMaterial->GetPenaltyCoefficient());

        dMatrixT B_K;
        Set_B_K(DNaX, B_K);

        for (int i = 0; i < nen; ++i) {
          for (int j = 0; j < nsd; j++) {
            R_K[(nen + i) * nsd + j] += DivPhi2K * B_K(0, i * nsd + j);
          }
        }

      }

    }

    R.CopyIn(0, RS);
    R.CopyIn(RS.Length(), RE);

    fRHS += R;
    fRHS += R_K;

  }

  //
  // extrapolate from integration points and compute output nodal/element
  // values
  //
  void FSPiezoElectricSolidT::ComputeOutput(const iArrayT& n_codes,
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
    int nstrs = StrainDim();

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
    dArray2DT ndElectricVectorPotential;
    dArray2DT ndDivergenceVectorPotential;
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

    ndElectricVectorPotential.Alias(nen, n_codes[ND_ELEC_POT], pall);
    pall += ndElectricVectorPotential.Length();

    ndDivergenceVectorPotential.Alias(nen, n_codes[ND_DIV_POT], pall);
    pall += ndDivergenceVectorPotential.Length();

    ndElectricDisplacement.Alias(nen, n_codes[ND_ELEC_DISP], pall);
    pall += ndElectricDisplacement.Length();

    ndElectricField.Alias(nen, n_codes[ND_ELEC_FLD], pall);
    pall += ndElectricField.Length();

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

      if (n_codes[ND_ELEC_POT]) {
        if (interpolant_DOF) {
          fLocVectorPotential.ReturnTranspose(ndElectricVectorPotential);
        } else {
          NodalDOFs(CurrentElement().NodesX(), ndElectricVectorPotential);
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

        // divergence of vector potential
        dArrayT DivPhi(1);
        DivPhi = fCurrMaterial->DivPhi();

        if (n_codes[ND_DIV_POT]) {
          if (qNoExtrap) {
            for (int k = 0; k < nen; k++) {
              ndDivergenceVectorPotential.AddToRowScaled(k, Na_X_ip_w(k, 0), DivPhi);
            }
          } else {
            fShapes->Extrapolate(DivPhi, ndDivergenceVectorPotential);
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
        ndElectricVectorPotential /= nip;
        ndDivergenceVectorPotential /= nip;
        ndElectricDisplacement /= nip;
        ndElectricField /= nip;
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

      nodal_all.BlockColumnCopyAt(ndElectricVectorPotential, colcount);
      colcount += ndElectricVectorPotential.MinorDim();

      nodal_all.BlockColumnCopyAt(ndDivergenceVectorPotential, colcount);
      colcount += ndDivergenceVectorPotential.MinorDim();

      nodal_all.BlockColumnCopyAt(ndElectricDisplacement, colcount);
      colcount += ndElectricDisplacement.MinorDim();

      nodal_all.BlockColumnCopyAt(ndElectricField, colcount);
      colcount += ndElectricField.MinorDim();

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
