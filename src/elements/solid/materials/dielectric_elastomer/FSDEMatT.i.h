
namespace Tahoe {

  inline FSDEMatT::FSDEMatT() :
    ParameterInterfaceT("Finite Strain Dielectric Elastomer"),
        fFSDEMatSupport(0)
  {
    SetName(FSDEMatT::Name);
    Initialize();
  }

  //
  // Set electrical permittivity
  //
  inline void FSDEMatT::SetElectricPermittivity(double epsilon)
  {
    fElectricPermittivity = epsilon;
  }

  //
  // Get electrical permittivity
  //
  inline double FSDEMatT::GetElectricPermittivity() const
  {
    return fElectricPermittivity;
  }

  //
  //
  //
  inline void FSDEMatT::SetFSDEMatSupport(
      const FSDEMatSupportT* support)
  {
    fFSDEMatSupport = support;
  }

  //
  //
  //
  inline const dSymMatrixT FSDEMatT::StressMechanical(
      const dSymMatrixT& C, const dArrayT& D) const
  {


//    dSymMatrixT Sm = Sevol;
//    Sm += Sedev;

//    return Sm;

  }

  //
  //
  //
  inline const dSymMatrixT FSDEMatT::StressElectrical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    dSymMatrixT Sr;


    return Sr;

  }

  //
  //
  //
  inline const dSymMatrixT FSDEMatT::StressElectromechanical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    dSymMatrixT Sz;



    return Sz;

  }

  //
  //
  //
  inline const dArrayT FSDEMatT::ElectricField(const dSymMatrixT& C,
      const dArrayT& D) const
  {

    const dArrayT Er = ElectricFieldElectrical(C, D);

    dArrayT E = Er;

    return E;

  }

  //
  //
  //
  inline const dArrayT FSDEMatT::ElectricFieldElectrical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    dArrayT Er;

	/* HSP - why is J necessary? */
    const double J = sqrt(C.Det());
    C.Multx(D, Er);
    Er /= (J * fElectricPermittivity);

    return Er;

  }

  //
  //
  //
  inline const dMatrixT FSDEMatT::TangentMechanical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

//    dMatrixT Cm = Cevol;

//    return Cm;

  }

  //
  //
  //
  inline const dMatrixT FSDEMatT::TangentElectrical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    dMatrixT beta;

    return beta;

  }

  //
  //
  //
  inline const dMatrixT FSDEMatT::TangentElectromechanical(
      const dSymMatrixT& C, const dArrayT& D) const
  {

    dMatrixT tangent;

    return tangent;

  }

  //
  //
  //
  inline const dArrayT FSDEMatT::ElectricDisplacement()
  {
    fElectricDisplacement = fFSDEMatSupport->ElectricDisplacement();
    return fElectricDisplacement;
  }

  //
  //
  //
  inline const dArrayT FSDEMatT::ElectricDisplacement(int ip)
  {
    fElectricDisplacement = fFSDEMatSupport->ElectricDisplacement(ip);
    return fElectricDisplacement;
  }

  //
  //
  //
  inline const dSymMatrixT FSDEMatT::RightCauchyGreenDeformation()
  {

    const dMatrixT F = F_mechanical();
    dMatrixT FTF(3);
    FTF.MultATB(F, F);
    dSymMatrixT C(3);
    C.Symmetrize(FTF);

    return C;

  }

  //
  // material energy density
  //
  inline double FSDEMatT::StrainEnergyDensity()
  {

//    fEnergyDensity = EnergyDensity(C, D);
	  return 0.0;
//    return fEnergyDensity;

  }

  //
  // material mechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMatT::C_IJKL()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const dArrayT D = ElectricDisplacement();
    fTangentMechanical = TangentMechanical(C, D);

    return fTangentMechanical;

  }

  //
  // material electromechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMatT::E_IJK()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const dArrayT D = ElectricDisplacement();
    fTangentElectromechanical = TangentElectromechanical(C, D);

    return fTangentElectromechanical;

  }

  //
  // material electric tangent modulus
  //
  inline const dMatrixT&
  FSDEMatT::B_IJ()
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
  FSDEMatT::S_IJ()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const dArrayT D = ElectricDisplacement();
    fStress = StressMechanical(C, D);

    return fStress;

  }

  //
  // Electric displacement
  //
  inline const dArrayT&
  FSDEMatT::D_I()
  {
    return fElectricDisplacement;
  }

  //
  // Electric field
  //
  inline const dArrayT&
  FSDEMatT::E_I()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const dArrayT D = ElectricDisplacement();
    fElectricField = ElectricField(C, D);

    return fElectricField;

  }

  //
  // spatial tangent modulus
  //
  inline const dMatrixT&
  FSDEMatT::c_ijkl()
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
  FSDEMatT::s_ij()
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
  inline double FSDEMatT::Pressure() const
  {

    return 0.0;

  }

  //
  // compute symmetric Cij reduced index matrix */
  //
  inline void FSDEMatT::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli) 
  {


  }

  //
  // compute symmetric 2nd PK2 reduced index vector */
  //
  inline void FSDEMatT::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2) 
  {


  }

  //
  // compute strain energy density for the specified strain */
  //
  inline double FSDEMatT::ComputeEnergyDensity(const dSymMatrixT& E) 
  {
	 return 0.0;

  }

} //namespace Tahoe
