#include "FSDE_inc_2D.h"

namespace Tahoe {

  inline FSDEMat2DT::FSDEMat2DT() :
    ParameterInterfaceT("Finite Strain Dielectric Elastomer 2D"),
        fFSDEMatSupport2D(0)
  {
    SetName(FSDEMat2DT::Name);
    Initialize();
  }

  //
  // Set electrical permittivity
  //
  inline void FSDEMat2DT::SetElectricPermittivity(double epsilon)
  {	
    fElectricPermittivity = epsilon;
  }

  //
  // Get electrical permittivity
  //
  inline double FSDEMat2DT::GetElectricPermittivity() const
  {
    return fElectricPermittivity;
  }

  //
  //
  //
  inline void FSDEMat2DT::SetFSDEMatSupport(
      const FSDEMatSupportT* support)
  {
    fFSDEMatSupport = support;
  }


  //
  //
  //
  inline const dArrayT FSDEMat2DT::ElectricField()
  {
    fElectricField = fFSDEMatSupport->ElectricField();
    return fElectricField;
  }

  //
  //
  //
  inline const dArrayT FSDEMat2DT::ElectricField(int ip)
  {
    fElectricField = fFSDEMatSupport->ElectricField(ip);
    return fElectricField;
  }

  //
  //
  //
  inline const dMatrixT FSDEMat2DT::RightCauchyGreenDeformation()
  {
    const dMatrixT F = F_mechanical();
    dMatrixT FTF(2);
    FTF.MultATB(F, F);

    return FTF;

  }

  //
  // material energy density
  //
  inline double FSDEMat2DT::StrainEnergyDensity()
  {

//    fEnergyDensity = EnergyDensity(C, D);
	  return 0.0;
//    return fEnergyDensity;

  }


  //
  // material mechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMat2DT::C_IJKL()
  {
	cout << "FSDEMat2DT::C_IJKL" << endl;
    const dMatrixT& C = RightCauchyGreenDeformation();
    const dArrayT& E = ElectricField();
 	
	/* call C function for mechanical tangent modulus */
 	get_ddC_2D(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), fTangentElectromechanical.Pointer()); 
	fTangentMechanical*=4.0;
	cout << "C_IJKL = " << fTangentMechanical << endl;
    return fTangentMechanical;
  }

  //
  // material electromechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMat2DT::E_IJK()
  {
	cout << "FSDEMat2DT::E_IJK" << endl;
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();

	/* call C function for electromechanical tangent modulus */
 	get_ddCE_2D(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), fTangentElectromechanical.Pointer()); 
 
 	fTangentElectromechanical*=-2.0;
 	cout << "E_IJK = " << fTangentElectromechanical << endl;
    return fTangentElectromechanical;

  }

  //
  // material electric tangent modulus
  //
  inline const dMatrixT&
  FSDEMat2DT::B_IJ()
  {
	cout << "FSDEMat2DT::B_IJ" << endl;
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();

	/* call C function for electrical tangent modulus */
 	get_ddE_2D(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), fTangentElectrical.Pointer()); 
 
 	fTangentElectrical *= -1.0;
 	cout << "B_IJ = " << fTangentElectrical << endl;
    return fTangentElectrical;

  }

  //
  // Second Piola-Kirchhoff stress (mechanical)
  //
  inline const dSymMatrixT&
  FSDEMat2DT::S_IJ()
  {
	cout << "FSDEMat2DT::S_IJ" << endl;
    const dMatrixT& C = RightCauchyGreenDeformation();
   	const dArrayT& E = ElectricField();
    
	/* call C function for mechanical stress */
	get_dUdCmech_2D(fParams.Pointer(), E.Pointer(),  
		C.Pointer(), stress_temp.Pointer()); 

    fStress.FromMatrix(stress_temp);
	fStress*=2.0;
	cout << "S_IJ = " << fStress << endl;
    return fStress;

  }

  //
  // Electric displacement - is it necessary to pass Efield?
  //
  inline const dArrayT&
  FSDEMat2DT::D_I()
  {
  	cout << "FSDEMat2DT::D_I" << endl;
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();

	/* call C function for electric displacement */
 	get_dUdE_2D(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), fElectricDisplacement.Pointer());     
    
    fElectricDisplacement *= -1.0;
    cout << "fElectricDisplacement = " << fElectricDisplacement << endl;
    return fElectricDisplacement;
  }

  //
  // Electric field
  //
  inline const dArrayT&
  FSDEMat2DT::E_I()
  {
    return fElectricField;
  }

  //
  // spatial tangent modulus
  //
  inline const dMatrixT&
  FSDEMat2DT::c_ijkl()
  {
	cout << "FSDEMat2DT::c_ijkl" << endl;
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
  FSDEMat2DT::s_ij()
  {
	cout << "FSDEMat2DT::s_ij" << endl;
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
  inline double FSDEMat2DT::Pressure() const
  {

    return 0.0;

  }

  //
  // compute symmetric Cij reduced index matrix */
  //
  inline void FSDEMat2DT::ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli) 
  {


  }

  //
  // compute symmetric 2nd PK2 reduced index vector */
  //
  inline void FSDEMat2DT::ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2) 
  {


  }

  //
  // compute strain energy density for the specified strain */
  //
  inline double FSDEMat2DT::ComputeEnergyDensity(const dSymMatrixT& E) 
  {
	 return 0.0;

  }

} //namespace Tahoe
