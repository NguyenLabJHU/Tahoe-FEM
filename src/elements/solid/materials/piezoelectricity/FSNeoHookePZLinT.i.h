//
// $Id: FSNeoHookePZLinT.i.h,v 1.1 2008-09-03 18:40:50 beichuan Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.2  2008/07/14 17:37:44  lxmota
// Various corrections related to initialization.
//
// Revision 1.1  2008/06/16 18:10:49  lxmota
// Piezoelectric material. Initial sources.
//
//

namespace Tahoe{

  inline
  FSNeoHookePZLinT::FSNeoHookePZLinT() :
    ParameterInterfaceT("Neohookean elastic linear piezoelectric"),
    fFSPZMatSupport(0)
  {
    SetName(FSNeoHookePZLinT::Name);
    initialize();
  }

  //
  //
  //
  inline bool
  FSNeoHookePZLinT::Need_D() const
  {

    return true;

  }

  //
  //
  //
  inline bool
  FSNeoHookePZLinT::Need_D_last() const
  {

    return false;

  }

  //
  // Set shear modulus
  //
  inline void
  FSNeoHookePZLinT::setShearModulus(double mu)
  {
    fShearModulus = mu;
  }

  //
  // Get shear modulus
  //
  inline double
  FSNeoHookePZLinT::getShearModulus() const
  {
    return fShearModulus;
  }

  //
  // Set bulk modulus
  //
  inline void
  FSNeoHookePZLinT::setBulkModulus(double kappa)
  {
    fBulkModulus = kappa;
  }

  //
  // Get bulk modulus
  //
  inline double
  FSNeoHookePZLinT::getBulkModulus() const
  {
    return fBulkModulus;
  }

  //
  // Set electrical permittivity
  //
  inline void
  FSNeoHookePZLinT::setElectricPermittivity(double epsilon)
  {
    fElectricPermittivity = epsilon;
  }

  //
  // Get electrical permittivity
  //
  inline double
  FSNeoHookePZLinT::getElectricPermittivity() const
  {
    return fElectricPermittivity;
  }

  //
  //
  //
  inline void
  FSNeoHookePZLinT::setFSPZMatSupport(const FSPZMatSupportT* support)
  {
    fFSPZMatSupport = support;
  }

  //
  //
  //
  inline const int
  FSNeoHookePZLinT::ManifoldDim() const
  {
    return FSPZMatSupportT::ManifoldDim();
  }

  //
  //
  //
  inline const int
  FSNeoHookePZLinT::StrainDim() const
  {
    return FSPZMatSupportT::StrainDim();
  }

  //
  //
  //
  inline const int
  FSNeoHookePZLinT::ElectricalDim() const
  {
    return FSPZMatSupportT::ElectricalDim();
  }

  //
  //
  //
  inline double
  FSNeoHookePZLinT::energyDensity(const dSymMatrixT& C, const dArrayT& D) const
  {

    const double Wm = energyDensityMechanical(C);
    const double Wr = energyDensityElectrical(C, D);
    const double Wz = energyDensityPiezoelectrical(C, D);

    return Wm + Wr + Wz;

  }

  //
  //
  //
  inline double
  FSNeoHookePZLinT::energyDensityMechanical(const dSymMatrixT& C) const
  {

    const double Wevol = energyDensityElasticVol(C);
    const double Wedev = energyDensityElasticDev(C);

    return Wevol + Wedev;

  }

  //
  //
  //
  inline double
  FSNeoHookePZLinT::energyDensityElasticVol(const dSymMatrixT& C) const
  {

    const double Jsqr     = C.Det();
    const double twotheta = log(Jsqr);
    const double Wevol    = 0.25 * fBulkModulus * (Jsqr - 1.0 - twotheta);

    return Wevol;

  }

  //
  //
  //
  inline double
  FSNeoHookePZLinT::energyDensityElasticDev(const dSymMatrixT& C) const
  {

    const double Jsqr = C.Det();
    const double Jm23 = pow(Jsqr, -1.0 / 3.0);

    dSymMatrixT Cbar(ManifoldDim());
    Cbar.SetToScaled(Jm23, C);

    const double Wedev = 0.5 * fShearModulus * (Cbar.Trace() - 3.0);

    return Wedev;

  }

  //
  //
  //
  inline double
  FSNeoHookePZLinT::energyDensityElectrical(const dSymMatrixT& C,
					    const dArrayT& D) const
  {

    const double J  = sqrt(C.Det());
    const double Wr = 0.5 * C.MultmBn(D, D) / J / fElectricPermittivity;

    return Wr;

  }

  //
  //
  //
  inline double
  FSNeoHookePZLinT::energyDensityPiezoelectrical(const dSymMatrixT& C,
						 const dArrayT& D) const
  {

    const dArrayT Ez = electricFieldPiezoelectrical(C, D);
    const double Wz  = 0.5 * dArrayT::Dot(D, Ez);

    return Wz;

  }

  //
  //
  //
  inline const dSymMatrixT
  FSNeoHookePZLinT::stress2PK(const dSymMatrixT& C, const dArrayT& D) const
  {

    const dSymMatrixT Sm = stressMechanical(C);
    const dSymMatrixT Sr = stressElectrical(C, D);
    const dSymMatrixT Sz = stressPiezoelectrical(C, D);

    dSymMatrixT S = Sm;
    S += Sr;
    S += Sz;

    return S;

  }

  //
  //
  //
  inline const dSymMatrixT
  FSNeoHookePZLinT::stressMechanical(const dSymMatrixT& C) const
  {

    const dSymMatrixT Sevol = stressElasticVol(C);
    const dSymMatrixT Sedev = stressElasticDev(C);

    dSymMatrixT Sm = Sevol;
    Sm += Sedev;

    return Sm;

  }

  //
  //
  //
  inline const dSymMatrixT
  FSNeoHookePZLinT::stressElasticVol(const dSymMatrixT& C) const
  {

    dSymMatrixT Sevol = C;
    Sevol.Inverse();

    const double Jsqr = C.Det();
    Sevol *= ( 0.5 * fBulkModulus * (Jsqr - 1.0) );

    return Sevol;

  }

  //
  //
  //
  inline const dSymMatrixT
  FSNeoHookePZLinT::stressElasticDev(const dSymMatrixT& C) const
  {

    dSymMatrixT Sedev = C;
    Sedev.Inverse();

    Sedev *= ( - C.Trace() / 3.0 );
    Sedev.PlusIdentity(1.0);

    const double Jsqr = C.Det();
    const double Jm23 = pow(Jsqr, -1.0 / 3.0);
    Sedev *= ( fShearModulus * Jm23 );

    return Sedev;

  }

  //
  //
  //
  inline const dSymMatrixT
  FSNeoHookePZLinT::stressElectrical(const dSymMatrixT& C,
				     const dArrayT& D) const
  {

    dMatrixT D_x_D( ElectricalDim() );
    D_x_D.Outer(D, D);

    dSymMatrixT Sr( ManifoldDim() );
    Sr.FromMatrix(D_x_D);

    const double J  = sqrt(C.Det());

    Sr *= ( 1.0 / J / fElectricPermittivity );

    const double Wr = energyDensityElectrical(C, D);
    dSymMatrixT WrCinv = C;
    WrCinv.Inverse();
    WrCinv *= Wr;

    Sr -= WrCinv;

    return Sr;

  }

  //
  //
  //
  inline const dSymMatrixT
  FSNeoHookePZLinT::stressPiezoelectrical(const dSymMatrixT& C,
					  const dArrayT& D) const
  {

    dArrayT DG( StrainDim() );

    for (int i = 0; i < StrainDim(); ++i) {

      DG[i] =
        fPiezoelectricTensor(0,i) * D[0] +
        fPiezoelectricTensor(1,i) * D[1] +
        fPiezoelectricTensor(2,i) * D[2];


    }

    dSymMatrixT Sz( ManifoldDim() );

    Sz(0,0) = DG[0];
    Sz(1,1) = DG[1];
    Sz(2,2) = DG[2];

    Sz(0,1) = Sz(1,0) = DG[5];
    Sz(0,2) = Sz(2,0) = DG[4];
    Sz(1,2) = Sz(2,1) = DG[3];

    return Sz;

  }

  //
  //
  //
  inline const dArrayT
  FSNeoHookePZLinT::electricField(const dSymMatrixT& C,
				  const dArrayT& D) const
  {

    const dArrayT Er = electricFieldElectrical(C, D);
    const dArrayT Ez = electricFieldPiezoelectrical(C, D);

    dArrayT E = Er;
    E+= Ez;

    return E;

  }

  //
  //
  //
  inline const dArrayT
  FSNeoHookePZLinT::electricFieldElectrical(const dSymMatrixT& C,
					    const dArrayT& D) const
  {

    dArrayT Er( ElectricalDim() );
    C.Multx(D, Er);

    const double J  = sqrt(C.Det());

    Er /= (J * fElectricPermittivity);

    return Er;

  }

  //
  //
  //
  inline const dArrayT
  FSNeoHookePZLinT::electricFieldPiezoelectrical(const dSymMatrixT& C,
      const dArrayT& D) const
  {

    dSymMatrixT E = C;
    E.PlusIdentity(-1.0);
    E *= 0.5;

    dArrayT Ez( ElectricalDim() );

    for (int i = 0; i < ElectricalDim(); ++i) {

      Ez[i] =
        fPiezoelectricTensor(i,0) * E(0,0) +
        fPiezoelectricTensor(i,1) * E(1,1) +
        fPiezoelectricTensor(i,2) * E(2,2) +
        fPiezoelectricTensor(i,3) * E(1,2) +
        fPiezoelectricTensor(i,4) * E(0,2) +
        fPiezoelectricTensor(i,5) * E(1,1);

    }

    return Ez;

  }

  //
  //
  //
  inline const dMatrixT
  FSNeoHookePZLinT::tangentMechanical(const dSymMatrixT& C,
				      const dArrayT& D) const
  {

    const dMatrixT Cevol = tangentMechanicalElasticVol(C);
    const dMatrixT Cedev = tangentMechanicalElasticDev(C);
    const dMatrixT Cr    = tangentMechanicalElectrical(C, D);

    dMatrixT Cm = Cevol;
    Cm += Cedev;
    Cm += Cr;

    return Cm;

  }

  //
  //
  //
  inline const dMatrixT
  FSNeoHookePZLinT::tangentMechanicalElasticVol(const dSymMatrixT& C) const
  {

    dSymMatrixT Cinv = C;
    Cinv.Inverse();

    dMatrixT Cinv_x_Cinv( StrainDim() );
    Cinv_x_Cinv.DyadAB(Cinv, Cinv);

    dMatrixT Cinv_o_Cinv( StrainDim() );
    Cinv_o_Cinv.ReducedI_C(Cinv);

    const double Jsqr = C.Det();
    Cinv_x_Cinv *= Jsqr;
    Cinv_o_Cinv *= (Jsqr - 1.0);
    dMatrixT Cevol = Cinv_x_Cinv;
    Cevol -= Cinv_o_Cinv;
    Cevol *= fBulkModulus;

    return Cevol;

  }

  //
  //
  //
  inline const dMatrixT
  FSNeoHookePZLinT::tangentMechanicalElasticDev(const dSymMatrixT& C) const
  {

    dSymMatrixT Cinv = C;
    Cinv.Inverse();

    dMatrixT Cinv_x_Cinv( StrainDim() );
    Cinv_x_Cinv.DyadAB(Cinv, Cinv);

    dMatrixT Cinv_o_Cinv( StrainDim() );
    Cinv_o_Cinv.ReducedI_C(Cinv);

    dSymMatrixT I = C;
    I.Identity(1.0);
    dMatrixT Cinv_x_I( StrainDim() );
    Cinv_x_I.DyadAB(Cinv, I);

    dMatrixT I_x_Cinv( StrainDim() );
    I_x_Cinv.DyadAB(I, Cinv);

    Cinv_x_Cinv /= 3.0;
    dMatrixT Cedev = Cinv_x_Cinv;
    Cedev += Cinv_o_Cinv;
    Cedev *= C.Trace();
    Cedev -= Cinv_x_I;
    Cedev -= I_x_Cinv;
    const double Jsqr = C.Det();
    const double Jm23 = pow(Jsqr, -1.0 / 3.0);
    Cedev *= (2.0 * fShearModulus * Jm23 / 3.0);

    return Cedev;

  }

  //
  //
  //
  inline const dMatrixT
  FSNeoHookePZLinT::tangentMechanicalElectrical(const dSymMatrixT& C,
						const dArrayT& D) const
  {

    dSymMatrixT Cinv = C;
    Cinv.Inverse();

    dMatrixT Cinv_x_Cinv( StrainDim() );
    Cinv_x_Cinv.DyadAB(Cinv, Cinv);


    dMatrixT Cinv_o_Cinv( StrainDim() );
    Cinv_o_Cinv.ReducedI_C(Cinv);

    dMatrixT DxD( ElectricalDim() );
    DxD.Outer(D, D);
    dSymMatrixT D_x_D( ManifoldDim() );
    D_x_D.FromMatrix(DxD);

    dMatrixT Cinv_x_D_x_D( StrainDim() );
    Cinv_x_D_x_D.DyadAB(Cinv, D_x_D);

    dMatrixT D_x_D_x_Cinv( StrainDim() );
    D_x_D_x_Cinv.DyadAB(D_x_D, Cinv);

    const double Wr = energyDensityElectrical(C, D);
    Cinv_x_Cinv *= Wr;
    dMatrixT Cr = Cinv_x_Cinv;
    Cinv_o_Cinv *= (2.0 * Wr);
    Cr += Cinv_o_Cinv;
    const double J  = sqrt(C.Det());
    Cinv_x_D_x_D /= (J * fElectricPermittivity);
    D_x_D_x_Cinv /= (J * fElectricPermittivity);
    Cr -= Cinv_x_D_x_D;
    Cr -= D_x_D_x_Cinv;

    return Cr;

  }

  //
  //
  //
  inline const dMatrixT
  FSNeoHookePZLinT::tangentElectrical(const dSymMatrixT& C,
				      const dArrayT& D) const
  {

    dMatrixT beta( ElectricalDim() );
    C.ToMatrix(beta);
    const double J  = sqrt(C.Det());
    beta /= (J * fElectricPermittivity);

    return beta;

  }

  //
  //
  //
  inline const dMatrixT
  FSNeoHookePZLinT::tangentPiezoelectrical(const dSymMatrixT& C,
					   const dArrayT& D) const
  {
    dSymMatrixT Cinv = C;
    Cinv.Inverse();

    dMatrixT C_x_Cinv( StrainDim() );
    C_x_Cinv.DyadAB(C, Cinv);

    dMatrixT S2mCoCinv( StrainDim() );
    S2mCoCinv.ReducedIndexI();
    S2mCoCinv *= 2.0;
    S2mCoCinv -= C_x_Cinv;

    const double J  = sqrt(C.Det());

    S2mCoCinv *= (1.0 / J / fElectricPermittivity);

    dMatrixT tangent( ElectricalDim(), StrainDim() );

    for (int j = 0; j < StrainDim(); ++j) {

      tangent(0,j) =
        D[0] * S2mCoCinv(0,j) + D[1] * S2mCoCinv(5,j) + D[2] * S2mCoCinv(4,j);

      tangent(1,j) =
        D[0] * S2mCoCinv(5,j) + D[1] * S2mCoCinv(1,j) + D[2] * S2mCoCinv(3,j);

      tangent(2,j) =
        D[0] * S2mCoCinv(4,j) + D[1] * S2mCoCinv(4,j) + D[2] * S2mCoCinv(2,j);

    }

    tangent += fPiezoelectricTensor;

    return fPiezoelectricTensor;

  }

  //
  //
  //
  inline const dArrayT
  FSNeoHookePZLinT::electricDisplacement()
  {

    fElectricDisplacement = fFSPZMatSupport->ElectricDisplacement();
    return fElectricDisplacement;

  }

  //
  //
  //
  inline const dArrayT
  FSNeoHookePZLinT::electricDisplacement(int ip)
  {

    fElectricDisplacement = fFSPZMatSupport->ElectricDisplacement(ip);
    return fElectricDisplacement;
  }

  //
  //
  //
  inline const dArrayT
  FSNeoHookePZLinT::electricDisplacement_last()
  {

    fElectricDisplacement = fFSPZMatSupport->ElectricDisplacement_last();
    return fElectricDisplacement;

  }

  //
  //
  //
  inline const dArrayT
  FSNeoHookePZLinT::electricDisplacement_last(int ip)
  {

    fElectricDisplacement = fFSPZMatSupport->ElectricDisplacement_last(ip);
    return fElectricDisplacement;

  }

  //
  //
  //
  inline const dSymMatrixT
  FSNeoHookePZLinT::rightCauchyGreenDeformation()
  {

    const dMatrixT F = F_mechanical();
    dMatrixT FTF( ManifoldDim() );
    FTF.MultATB(F, F);
    dSymMatrixT C( ManifoldDim() );
    C.Symmetrize(FTF);

    return C;

  }

  //
  // material energy density
  //
  inline double
  FSNeoHookePZLinT::StrainEnergyDensity()
  {

    const dSymMatrixT C = rightCauchyGreenDeformation();
    const dArrayT D     = electricDisplacement();
    fEnergyDensity      = energyDensity(C, D);

    return fEnergyDensity;

  }

  //
  // material mechanical tangent modulus
  //
  inline const dMatrixT&
  FSNeoHookePZLinT::C_ijkl()
  {

    const dSymMatrixT C  = rightCauchyGreenDeformation();
    const dArrayT D      = electricDisplacement();
    fTangentMechanical   = tangentMechanical(C, D);

    return fTangentMechanical;

  }

  //
  // material piezoelectric tangent modulus
  //
  inline const dMatrixT&
  FSNeoHookePZLinT::H_ijk()
  {

    const dSymMatrixT C     = rightCauchyGreenDeformation();
    const dArrayT D         = electricDisplacement();
    fTangentPiezoelectrical = tangentPiezoelectrical(C, D);

    return fTangentPiezoelectrical;

  }

  //
  // material electric tangent modulus
  //
  inline const dMatrixT&
  FSNeoHookePZLinT::B_ij()
  {

    const dSymMatrixT C  = rightCauchyGreenDeformation();
    const dArrayT D      = electricDisplacement();
    fTangentElectrical   = tangentElectrical(C, D);

    return fTangentElectrical;

  }

  //
  // Second Piola-Kirchhoff stress
  //
  inline const dSymMatrixT&
  FSNeoHookePZLinT::S_ij()
  {

    const dSymMatrixT C = rightCauchyGreenDeformation();
    const dArrayT D     = electricDisplacement();
    fStress             = stress2PK(C, D);

    return fStress;

  }

  //
  // Electric field
  //
  inline const dArrayT&
  FSNeoHookePZLinT::E_i()
  {

    const dSymMatrixT C = rightCauchyGreenDeformation();
    const dArrayT D     = electricDisplacement();
    fElectricField      = electricField(C, D);

    return fElectricField;

  }

  //
  // spatial tangent modulus
  //
  inline const dMatrixT&
  FSNeoHookePZLinT::c_ijkl()
  {

    const dMatrixT F     = F_mechanical();
    const double J       = F.Det();

    // prevent aliasing
    const dMatrixT Cijkl = C_ijkl();
    fTangentMechanical.SetToScaled( 1.0 / J, PushForward(F, Cijkl) );

    return fTangentMechanical;

  }

  //
  // Cauchy stress
  //
  inline const dSymMatrixT&
  FSNeoHookePZLinT::s_ij()
  {

    const dMatrixT F     = F_mechanical();
    const double J       = F.Det();

    // prevent aliasing
    const dSymMatrixT S  = S_ij();
    fStress.SetToScaled( 1.0 / J, PushForward(F, S) );

    return fStress;

  }

  //
  // pressure associated with the last computed stress
  //
  inline double
  FSNeoHookePZLinT::Pressure() const
  {

    return fStress.Trace() / 3.0;

  }

} //namespace Tahoe
