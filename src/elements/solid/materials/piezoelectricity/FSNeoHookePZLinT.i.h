//
// $Id: FSNeoHookePZLinT.i.h,v 1.2 2008-12-12 18:58:15 amota Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.1  2008/09/03 18:40:50  beichuan
// Piezoelectricity
//
// Revision 1.2  2008/07/14 17:37:44  lxmota
// Various corrections related to initialization.
//
// Revision 1.1  2008/06/16 18:10:49  lxmota
// Piezoelectric material. Initial sources.
//
//

namespace Tahoe {

  inline FSNeoHookePZLinT::FSNeoHookePZLinT() :
    ParameterInterfaceT("Neohookean elastic linear piezoelectric"),
        fFSPZMatSupport(0)
  {
    SetName(FSNeoHookePZLinT::Name);
    Initialize();
  }

  //
  // Set electrical permittivity
  //
  inline void FSNeoHookePZLinT::SetElectricPermittivity(double epsilon)
  {
    fElectricPermittivity = epsilon;
  }

  //
  // Get electrical permittivity
  //
  inline double FSNeoHookePZLinT::GetElectricPermittivity() const
  {
    return fElectricPermittivity;
  }

  //
  // Set penalty coefficient
  //
  inline void FSNeoHookePZLinT::SetPenaltyCoeffifient(double k)
  {
    fPenaltyCoefficient = k;
  }

  //
  // Get penalty coefficient
  //
  inline double FSNeoHookePZLinT::GetPenaltyCoefficient() const
  {
    return fPenaltyCoefficient;
  }

  //
  //
  //
  inline void FSNeoHookePZLinT::SetFSPZMatSupport(
      const FSPZMatSupportT* support)
  {
    fFSPZMatSupport = support;
  }

  //
  //
  //
  inline const int FSNeoHookePZLinT::ManifoldDim() const
  {
    return FSPZMatSupportT::ManifoldDim();
  }

  //
  //
  //
  inline const int FSNeoHookePZLinT::StrainDim() const
  {
    return FSPZMatSupportT::StrainDim();
  }

  //
  //
  //
  inline const int FSNeoHookePZLinT::ElectricalDim() const
  {
    return FSPZMatSupportT::ElectricalDim();
  }

  //
  //
  //
  inline double FSNeoHookePZLinT::EnergyDensity(const dSymMatrixT& C,
      const dArrayT& D) const
  {

    const double Wm = EnergyDensityMechanical(C);
    const double Wr = EnergyDensityElectrical(C, D);
    const double Wz = EnergyDensityPiezoelectrical(C, D);

    return Wm + Wr + Wz;

  }

  //
  //
  //
  inline double FSNeoHookePZLinT::EnergyDensityMechanical(const dSymMatrixT& C) const
  {

    const double Wevol = EnergyDensityElasticVol(C);
    const double Wedev = EnergyDensityElasticDev(C);

    return Wevol + Wedev;

  }

  //
  //
  //
  inline double FSNeoHookePZLinT::EnergyDensityElasticVol(const dSymMatrixT& C) const
  {

    const double Jsqr = C.Det();
    const double twotheta = log(Jsqr);
    const double Wevol = 0.25 * fKappa * (Jsqr - 1.0 - twotheta);

    return Wevol;

  }

  //
  //
  //
  inline double FSNeoHookePZLinT::EnergyDensityElasticDev(const dSymMatrixT& C) const
  {

    const double Jsqr = C.Det();
    const double Jm23 = pow(Jsqr, -1.0 / 3.0);

    dSymMatrixT Cbar(ManifoldDim());
    Cbar.SetToScaled(Jm23, C);

    const double Wedev = 0.5 * fMu * (Cbar.Trace() - 3.0);

    return Wedev;

  }

  //
  //
  //
  inline double FSNeoHookePZLinT::EnergyDensityElectrical(const dSymMatrixT& C,
      const dArrayT& D) const
  {

    double Wr = 0.0;

    if (dependsOnC == true) {

      const double J = sqrt(C.Det());
      Wr = 0.5 * C.MultmBn(D, D) / J / fElectricPermittivity;

    } else {

      Wr = 0.5 * dArrayT::Dot(D,D) / fElectricPermittivity;

    }

    return Wr;

  }

  //
  //
  //
  inline double FSNeoHookePZLinT::EnergyDensityPiezoelectrical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    const dArrayT Ez = ElectricFieldPiezoelectrical(C, D);
    const double Wz = 0.5 * dArrayT::Dot(D, Ez);

    return Wz;

  }

  //
  //
  //
  inline const dSymMatrixT FSNeoHookePZLinT::Stress(const dSymMatrixT& C,
      const dArrayT& D) const
  {

    const dSymMatrixT Sm = StressMechanical(C);
    const dSymMatrixT Sr = StressElectrical(C, D);
    const dSymMatrixT Sz = StressPiezoelectrical(C, D);

    dSymMatrixT S = Sm;
    S += Sr;
    S += Sz;

    return S;

  }

  //
  //
  //
  inline const dSymMatrixT FSNeoHookePZLinT::StressMechanical(
      const dSymMatrixT& C) const
  {

    const dSymMatrixT Sevol = StressElasticVol(C);
    const dSymMatrixT Sedev = StressElasticDev(C);

    dSymMatrixT Sm = Sevol;
    Sm += Sedev;

    return Sm;

  }

  //
  //
  //
  inline const dSymMatrixT FSNeoHookePZLinT::StressElasticVol(
      const dSymMatrixT& C) const
  {

    dSymMatrixT Sevol = C;
    Sevol.Inverse();

    const double Jsqr = C.Det();
    Sevol *= (0.5 * fKappa * (Jsqr - 1.0));

    return Sevol;

  }

  //
  //
  //
  inline const dSymMatrixT FSNeoHookePZLinT::StressElasticDev(
      const dSymMatrixT& C) const
  {

    dSymMatrixT Sedev = C;
    Sedev.Inverse();

    Sedev *= (-C.Trace() / 3.0);
    Sedev.PlusIdentity(1.0);

    const double Jsqr = C.Det();
    const double Jm23 = pow(Jsqr, -1.0 / 3.0);
    Sedev *= (fMu * Jm23);

    return Sedev;

  }

  //
  //
  //
  inline const dSymMatrixT FSNeoHookePZLinT::StressElectrical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    dSymMatrixT Sr(ManifoldDim());

    if (dependsOnC == true) {

      dMatrixT D_x_D(ElectricalDim());
      D_x_D.Outer(D, D);

      Sr.FromMatrix(D_x_D);

      const double J = sqrt(C.Det());

      Sr *= (1.0 / J / fElectricPermittivity);

      const double Wr = EnergyDensityElectrical(C, D);
      dSymMatrixT WrCinv = C;
      WrCinv.Inverse();
      WrCinv *= Wr;

      Sr -= WrCinv;

    } else {

      Sr = 0.0;

    }

    return Sr;

  }

  //
  //
  //
  inline const dSymMatrixT FSNeoHookePZLinT::StressPiezoelectrical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    dArrayT DG(StrainDim());

    for (int i = 0; i < StrainDim(); ++i) {

      DG[i] = fPiezoelectricTensor(0, i) * D[0] + fPiezoelectricTensor(1, i)
          * D[1] + fPiezoelectricTensor(2, i) * D[2];

    }

    dSymMatrixT Sz(ManifoldDim());

    Sz(0, 0) = DG[0];
    Sz(1, 1) = DG[1];
    Sz(2, 2) = DG[2];

    Sz(0, 1) = Sz(1, 0) = DG[5];
    Sz(0, 2) = Sz(2, 0) = DG[4];
    Sz(1, 2) = Sz(2, 1) = DG[3];

    return Sz;

  }

  //
  //
  //
  inline const dArrayT FSNeoHookePZLinT::ElectricField(const dSymMatrixT& C,
      const dArrayT& D) const
  {

    const dArrayT Er = ElectricFieldElectrical(C, D);
    const dArrayT Ez = ElectricFieldPiezoelectrical(C, D);

    dArrayT E = Er;
    E += Ez;

    return E;

  }

  //
  //
  //
  inline const dArrayT FSNeoHookePZLinT::ElectricFieldElectrical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    dArrayT Er(ElectricalDim());

    if (dependsOnC == true) {

      const double J = sqrt(C.Det());
      C.Multx(D, Er);
      Er /= (J * fElectricPermittivity);

    } else {

      Er = D;
      Er /= fElectricPermittivity;

    }

    return Er;

  }

  //
  //
  //
  inline const dArrayT FSNeoHookePZLinT::ElectricFieldPiezoelectrical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    dSymMatrixT E = C;
    E.PlusIdentity(-1.0);
    E *= 0.5;

    dArrayT Ez(ElectricalDim());

    for (int i = 0; i < ElectricalDim(); ++i) {

      Ez[i] = fPiezoelectricTensor(i, 0) * E(0, 0) + fPiezoelectricTensor(i, 1)
          * E(1, 1) + fPiezoelectricTensor(i, 2) * E(2, 2)
          + fPiezoelectricTensor(i, 3) * E(1, 2) + fPiezoelectricTensor(i, 4)
          * E(0, 2) + fPiezoelectricTensor(i, 5) * E(1, 1);

    }

    return Ez;

  }

  //
  //
  //
  inline const dMatrixT FSNeoHookePZLinT::TangentMechanical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    const dMatrixT Cevol = TangentMechanicalElasticVol(C);
    const dMatrixT Cedev = TangentMechanicalElasticDev(C);
    const dMatrixT Cr = TangentMechanicalElectrical(C, D);

    dMatrixT Cm = Cevol;
    Cm += Cedev;
    Cm += Cr;

    return Cm;

  }

  //
  //
  //
  inline const dMatrixT FSNeoHookePZLinT::TangentMechanicalElasticVol(
      const dSymMatrixT& C) const
  {

    dSymMatrixT Cinv = C;
    Cinv.Inverse();

    dMatrixT Cinv_x_Cinv(StrainDim());
    Cinv_x_Cinv.DyadAB(Cinv, Cinv);

    dMatrixT Cinv_o_Cinv(StrainDim());
    Cinv_o_Cinv.ReducedI_C(Cinv);

    const double Jsqr = C.Det();
    Cinv_x_Cinv *= Jsqr;
    Cinv_o_Cinv *= (Jsqr - 1.0);
    dMatrixT Cevol = Cinv_x_Cinv;
    Cevol -= Cinv_o_Cinv;
    Cevol *= fKappa;

    return Cevol;

  }

  //
  //
  //
  inline const dMatrixT FSNeoHookePZLinT::TangentMechanicalElasticDev(
      const dSymMatrixT& C) const
  {

    dSymMatrixT Cinv = C;
    Cinv.Inverse();

    dMatrixT Cinv_x_Cinv(StrainDim());
    Cinv_x_Cinv.DyadAB(Cinv, Cinv);

    dMatrixT Cinv_o_Cinv(StrainDim());
    Cinv_o_Cinv.ReducedI_C(Cinv);

    dSymMatrixT I = C;
    I.Identity(1.0);
    dMatrixT Cinv_x_I(StrainDim());
    Cinv_x_I.DyadAB(Cinv, I);

    dMatrixT I_x_Cinv(StrainDim());
    I_x_Cinv.DyadAB(I, Cinv);

    Cinv_x_Cinv /= 3.0;
    dMatrixT Cedev = Cinv_x_Cinv;
    Cedev += Cinv_o_Cinv;
    Cedev *= C.Trace();
    Cedev -= Cinv_x_I;
    Cedev -= I_x_Cinv;
    const double Jsqr = C.Det();
    const double Jm23 = pow(Jsqr, -1.0 / 3.0);
    Cedev *= (2.0 * fMu * Jm23 / 3.0);

    return Cedev;

  }

  //
  //
  //
  inline const dMatrixT FSNeoHookePZLinT::TangentMechanicalElectrical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    dMatrixT Cr(StrainDim());

    if (dependsOnC == true) {

      dSymMatrixT Cinv = C;
      Cinv.Inverse();

      dMatrixT Cinv_x_Cinv(StrainDim());
      Cinv_x_Cinv.DyadAB(Cinv, Cinv);

      dMatrixT Cinv_o_Cinv(StrainDim());
      Cinv_o_Cinv.ReducedI_C(Cinv);

      dMatrixT DxD(ElectricalDim());
      DxD.Outer(D, D);
      dSymMatrixT D_x_D(ManifoldDim());
      D_x_D.FromMatrix(DxD);

      dMatrixT Cinv_x_D_x_D(StrainDim());
      Cinv_x_D_x_D.DyadAB(Cinv, D_x_D);

      dMatrixT D_x_D_x_Cinv(StrainDim());
      D_x_D_x_Cinv.DyadAB(D_x_D, Cinv);

      const double Wr = EnergyDensityElectrical(C, D);
      Cinv_x_Cinv *= Wr;
      Cr = Cinv_x_Cinv;
      Cinv_o_Cinv *= (2.0 * Wr);
      Cr += Cinv_o_Cinv;
      const double J = sqrt(C.Det());
      Cinv_x_D_x_D /= (J * fElectricPermittivity);
      D_x_D_x_Cinv /= (J * fElectricPermittivity);
      Cr -= Cinv_x_D_x_D;
      Cr -= D_x_D_x_Cinv;

    } else {

      Cr = 0.0;

    }

    return Cr;

  }

  //
  //
  //
  inline const dMatrixT FSNeoHookePZLinT::TangentElectrical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    dMatrixT beta(ElectricalDim());

    if (dependsOnC == true) {

      C.ToMatrix(beta);
      const double J = sqrt(C.Det());
      beta /= (J * fElectricPermittivity);

    } else {

      beta.Identity(1.0 / fElectricPermittivity);

    }

    return beta;

  }

  //
  //
  //
  inline const dMatrixT FSNeoHookePZLinT::TangentPiezoelectrical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    dMatrixT tangent(ElectricalDim(), StrainDim());

    if (dependsOnC == true) {

      dSymMatrixT Cinv = C;
      Cinv.Inverse();

      dMatrixT C_x_Cinv(StrainDim());
      C_x_Cinv.DyadAB(C, Cinv);

      dMatrixT S2mCoCinv(StrainDim());
      S2mCoCinv.ReducedIndexI();
      S2mCoCinv *= 2.0;
      S2mCoCinv -= C_x_Cinv;

      const double J = sqrt(C.Det());

      S2mCoCinv *= (1.0 / J / fElectricPermittivity);

      for (int j = 0; j < StrainDim(); ++j) {

        tangent(0, j) = D[0] * S2mCoCinv(0, j) + D[1] * S2mCoCinv(5, j) + D[2]
            * S2mCoCinv(4, j);

        tangent(1, j) = D[0] * S2mCoCinv(5, j) + D[1] * S2mCoCinv(1, j) + D[2]
            * S2mCoCinv(3, j);

        tangent(2, j) = D[0] * S2mCoCinv(4, j) + D[1] * S2mCoCinv(4, j) + D[2]
            * S2mCoCinv(2, j);

      }

      tangent += fPiezoelectricTensor;

    } else {

      tangent = fPiezoelectricTensor;

    }

    return tangent;

  }

  //
  //
  //
  inline const dArrayT FSNeoHookePZLinT::ElectricDisplacement()
  {
    fElectricDisplacement = fFSPZMatSupport->ElectricDisplacement();
    return fElectricDisplacement;
  }

  //
  //
  //
  inline const dArrayT FSNeoHookePZLinT::ElectricDisplacement(int ip)
  {
    fElectricDisplacement = fFSPZMatSupport->ElectricDisplacement(ip);
    return fElectricDisplacement;
  }

  //
  //
  //
  inline double FSNeoHookePZLinT::DivergenceVectorPotential()
  {
    fDivergenceVectorPotential = fFSPZMatSupport->DivergenceVectorPotential();
    return fDivergenceVectorPotential;
  }

  //
  //
  //
  inline double FSNeoHookePZLinT::DivergenceVectorPotential(int ip)
  {
    fDivergenceVectorPotential = fFSPZMatSupport->DivergenceVectorPotential(ip);
    return fDivergenceVectorPotential;
  }

  //
  //
  //
  inline const dSymMatrixT FSNeoHookePZLinT::RightCauchyGreenDeformation()
  {

    const dMatrixT F = F_mechanical();
    dMatrixT FTF(ManifoldDim());
    FTF.MultATB(F, F);
    dSymMatrixT C(ManifoldDim());
    C.Symmetrize(FTF);

    return C;

  }

  //
  // material energy density
  //
  inline double FSNeoHookePZLinT::StrainEnergyDensity()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const dArrayT D = ElectricDisplacement();
    fEnergyDensity = EnergyDensity(C, D);

    return fEnergyDensity;

  }

  //
  // material mechanical tangent modulus
  //
  inline const dMatrixT&
  FSNeoHookePZLinT::C_IJKL()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const dArrayT D = ElectricDisplacement();
    fTangentMechanical = TangentMechanical(C, D);

    return fTangentMechanical;

  }

  //
  //
  //
  inline const dMatrixT&
  FSNeoHookePZLinT::C_IJKL_Num()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const dArrayT D = ElectricDisplacement();
    fTangentMechanical = TangentMechanicalNum(C, D);

    return fTangentMechanical;

  }

  //
  // material piezoelectric tangent modulus
  //
  inline const dMatrixT&
  FSNeoHookePZLinT::H_IJK()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const dArrayT D = ElectricDisplacement();
    fTangentPiezoelectrical = TangentPiezoelectrical(C, D);

    return fTangentPiezoelectrical;

  }

  //
  // material electric tangent modulus
  //
  inline const dMatrixT&
  FSNeoHookePZLinT::B_IJ()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const dArrayT D = ElectricDisplacement();
    fTangentElectrical = TangentElectrical(C, D);

    return fTangentElectrical;

  }

  //
  // Second Piola-Kirchhoff stress
  //
  inline const dSymMatrixT&
  FSNeoHookePZLinT::S_IJ()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const dArrayT D = ElectricDisplacement();
    fStress = Stress(C, D);

    return fStress;

  }

  //
  //
  //
  inline const dSymMatrixT&
  FSNeoHookePZLinT::S_IJ_Num()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const dArrayT D = ElectricDisplacement();
    fStress = StressNum(C, D);

    return fStress;

  }


  //
  // Second Piola-Kirchhoff stress
  //
  inline const dSymMatrixT&
  FSNeoHookePZLinT::S_IJ(const dSymMatrixT& C)
  {

    const dArrayT D = ElectricDisplacement();
    fStress = Stress(C, D);

    return fStress;

  }

  //
  // Electric displacement
  //
  inline const dArrayT&
  FSNeoHookePZLinT::D_I()
  {
    return fElectricDisplacement;
  }

  //
  // Electric field
  //
  inline const dArrayT&
  FSNeoHookePZLinT::E_I()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const dArrayT D = ElectricDisplacement();
    fElectricField = ElectricField(C, D);

    return fElectricField;

  }

  //
  // Divergence of vector potential
  //
  inline double FSNeoHookePZLinT::DivPhi()
  {
    return DivergenceVectorPotential();
  }

  //
  // spatial tangent modulus
  //
  inline const dMatrixT&
  FSNeoHookePZLinT::c_ijkl()
  {

    const dMatrixT F = F_mechanical();
    const double J = F.Det();

    // prevent aliasing
    const dMatrixT CIJKL = C_IJKL();
    fTangentMechanical.SetToScaled(1.0 / J, PushForward(F, CIJKL));

    return fTangentMechanical;

  }

  //
  // Cauchy stress
  //
  inline const dSymMatrixT&
  FSNeoHookePZLinT::s_ij()
  {

    const dMatrixT F = F_mechanical();
    const double J = F.Det();

    // prevent aliasing
    const dSymMatrixT S = S_IJ();
    fStress.SetToScaled(1.0 / J, PushForward(F, S));

    return fStress;

  }

  //
  // pressure associated with the last computed stress
  //
  inline double FSNeoHookePZLinT::Pressure() const
  {

    return fStress.Trace() / 3.0;

  }

} //namespace Tahoe
