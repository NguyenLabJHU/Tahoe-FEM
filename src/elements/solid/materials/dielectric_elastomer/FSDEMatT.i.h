#include "FSDE_inc.h"

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
  inline const dArrayT FSDEMatT::ElectricField()
  {
    fElectricField = fFSDEMatSupport->ElectricField();
    return fElectricField;
  }

  //
  //
  //
  inline const dArrayT FSDEMatT::ElectricField(int ip)
  {
    fElectricField = fFSDEMatSupport->ElectricField(ip);
    return fElectricField;
  }

  //
  //
  //
  inline const dMatrixT FSDEMatT::RightCauchyGreenDeformation()
  {
    const dMatrixT& F = F_mechanical();
    dMatrixT FTF(3);
    FTF.MultATB(F, F);

    return FTF;

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


  /* Mechanical and electromechanical tangent modulus */
  inline void FSDEMatT::C_Mech_Elec(dMatrixT& mech, dMatrixT& elec)
  {
// //  	cout << "FSDEMatT::C_Mech_Elec" << endl;
//   	const dMatrixT& C = RightCauchyGreenDeformation();
//   	const dArrayT& E = ElectricField();
//  	cout << "C_Mech_Elec C = " << C << endl;
  
  	/* call C function for both tangent moduli */
//   	get_ddCmech_elec(fParams.Pointer(), E.Pointer(), C.Pointer(),
//   		mech.Pointer(), elec.Pointer());
//   
//   	mech *= 4.0;
//   	elec *= -2.0;
// 	cout << "mech = " << mech << endl;
// 	cout << "elec = " << elec << endl;
  }

  /* Electrical displacement and tangent modulus */
  inline void FSDEMatT::S_C_Elec(dArrayT& D, dMatrixT& CE)
  {
//  	cout << "FSDEMatT::S_C_Elec" << endl;
//   	const dMatrixT& C = RightCauchyGreenDeformation();
//   	const dArrayT& E = ElectricField();
//  	cout << "S_C_Elec C = " << C << endl;
  
  	/* call C function for both tangent moduli */
//   	get_ddC_sc_elec(fParams.Pointer(), E.Pointer(), C.Pointer(),
//   		D.Pointer(), CE.Pointer());
//   
//   	D *= -1.0;
//   	CE *= -1.0;
  }

  //
  // material mechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMatT::C_IJKL()
  {
     const dMatrixT& C = RightCauchyGreenDeformation();
     const dArrayT& E = ElectricField();
     
	/* call C function for mechanical tangent modulus */
	get_ddCmech(fParams.Pointer(), E.Pointer(),  
		C.Pointer(), fTangentMechanical.Pointer()); 

	fTangentMechanical*=4.0;
	
    return fTangentMechanical;
  }

  //
  // material electromechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMatT::E_IJK()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();

	/* call C function for electromechanical tangent modulus */
 	get_ddCE(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), fTangentElectromechanical.Pointer()); 
 
 	fTangentElectromechanical*=-2.0;
    return fTangentElectromechanical;

  }

  //
  // material electric tangent modulus
  //
  inline const dMatrixT&
  FSDEMatT::B_IJ()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();

	/* call C function for electrical tangent modulus */
 	get_ddE(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), fTangentElectrical.Pointer()); 
 
 	fTangentElectrical *= -1.0;
    return fTangentElectrical;

  }

  //
  // Second Piola-Kirchhoff stress (mechanical)
  //
  inline const dSymMatrixT&
  FSDEMatT::S_IJ()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
   	const dArrayT& E = ElectricField();
    dMatrixT stress_temp(3);
 
	/* call C function for mechanical stress */
	get_dUdCmech(fParams.Pointer(), E.Pointer(),  
		C.Pointer(), stress_temp.Pointer()); 

    fStress.FromMatrix(stress_temp);
	fStress*=2.0;

    return fStress;

  }

  //
  // Electric displacement - is it necessary to pass Efield?
  //
  inline const dArrayT&
  FSDEMatT::D_I()
  {
  	const dMatrixT& C = RightCauchyGreenDeformation();
  	const dArrayT& E = ElectricField();
  	dMatrixT CE(3);
  
  	/* call C function for both tangent moduli */
  	get_dUdE(fParams.Pointer(), E.Pointer(), C.Pointer(),
  		fElectricDisplacement.Pointer());
  
  	fElectricDisplacement *= -1.0;
  	return fElectricDisplacement;
  }

  //
  // Electric field
  //
  inline const dArrayT&
  FSDEMatT::E_I()
  {
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
