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

  //
  // material mechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMatT::C_IJKL()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
    const dArrayT& E = ElectricField();

	double I1 = C(0,0)+C(1,1)+C(2,2);
	double J = C.Det();
	J = sqrt(J);
	
	/* call C function for mechanical part of tangent modulus */
 	mech_tanmod_ab(fParams.Pointer(), E.Pointer(), C.Pointer(), J, I1, fTangentMechanical.Pointer()); 
 	me_tanmod_ab(fParams.Pointer(), E.Pointer(), C.Pointer(), J, fTangentMechanicalElec.Pointer());
 	fTangentMechanical+=fTangentMechanicalElec;
    return fTangentMechanical;
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
	dMatrixT stress_temp2(3);
	double I1 = C(0,0)+C(1,1)+C(2,2);
	double J = C.Det();
	J = sqrt(J);
	
	/* call C function for mechanical part of PK2 stress */
 	mech_pk2_ab(fParams.Pointer(), E.Pointer(), C.Pointer(), J, I1, stress_temp.Pointer()); 
	me_pk2_ab(fParams.Pointer(), E.Pointer(), C.Pointer(), J, stress_temp2.Pointer());
	stress_temp+=stress_temp2;
	
	fStress.FromMatrix(stress_temp);
    return fStress;
  }

  //
  // material electromechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMatT::E_IJK()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
	double J = C.Det();
	J = sqrt(J);

	/* call C function for electromechanical tangent modulus */
 	me_mixedmodulus_ab(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), J, fTangentElectromechanical.Pointer()); 
 
    return fTangentElectromechanical;

  }

  //
  // material electric tangent modulus
  //
  inline const dMatrixT&
  FSDEMatT::B_IJ()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	double J = C.Det();
	J = sqrt(J);

	dMatrixT Cinv(3);
	Cinv.Inverse(C);
	fTangentElectrical = Cinv;
	fTangentElectrical *= fElectricPermittivity;
	fTangentElectrical *= J;
    return fTangentElectrical;

  }

  //
  // Electric displacement 
  //
  inline const dArrayT&
  FSDEMatT::D_I()
  {
  	const dMatrixT& C = RightCauchyGreenDeformation();
  	const dArrayT& E = ElectricField();
  	
  	double J = C.Det();
  	J = sqrt(J);
  
	/* call C function for electric stress (i.e. electric displacement D_{I}) */
 	elec_pk2(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), J, fElectricDisplacement.Pointer()); 
 		
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

}	// namespace Tahoe