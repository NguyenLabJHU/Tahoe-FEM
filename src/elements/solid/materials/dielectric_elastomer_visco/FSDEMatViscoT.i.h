#include "FSDEVisco_inc.h"

namespace Tahoe {

  inline FSDEMatViscoT::FSDEMatViscoT() :
    ParameterInterfaceT("Finite Strain Dielectric Elastomer Viscoelastic"),
        fFSDEMatSupportVisco(0)
  {
    SetName(FSDEMatViscoT::Name);
    Initialize();
  }

  //
  // Set electrical permittivity
  //
  inline void FSDEMatViscoT::SetElectricPermittivity(double epsilon)
  {	
    fElectricPermittivity = epsilon;
  }

  //
  // Get electrical permittivity
  //
  inline double FSDEMatViscoT::GetElectricPermittivity() const
  {
    return fElectricPermittivity;
  }

  //
  //
  //
  inline void FSDEMatViscoT::SetFSDEMatSupportVisco(
      const FSDEMatSupportViscoT* support)
  {
    fFSDEMatSupportVisco = support;
  }


  //
  //
  //
  inline const dArrayT FSDEMatViscoT::ElectricField()
  {
    fElectricField = fFSDEMatSupportVisco->ElectricField();
    return fElectricField;
  }

  //
  //
  //
  inline const dArrayT FSDEMatViscoT::ElectricField(int ip)
  {
    fElectricField = fFSDEMatSupportVisco->ElectricField(ip);
    return fElectricField;
  }

  //
  //
  //
  inline const dMatrixT FSDEMatViscoT::RightCauchyGreenDeformation()
  {
    const dMatrixT& F = F_mechanical();
    dMatrixT FTF(3);
    FTF.MultATB(F, F);

    return FTF;

  }

  //
  // material energy density
  //
  inline double FSDEMatViscoT::StrainEnergyDensity()
  {

//    fEnergyDensity = EnergyDensity(C, D);
	  return 0.0;
//    return fEnergyDensity;

  }

  //
  // material mechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMatViscoT::C_IJKL()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
    const dArrayT& E = ElectricField();

	double I1 = C(0,0)+C(1,1)+C(2,2);
	double J = C.Det();
	J = sqrt(J);
	
	/* call C function for mechanical part of tangent modulus */
 	mech_tanmod_ab_visco(fParams.Pointer(), E.Pointer(), C.Pointer(), J, I1, fTangentMechanical.Pointer()); 
 	me_tanmod_ab_visco(fParams.Pointer(), E.Pointer(), C.Pointer(), J, fTangentMechanicalElec.Pointer());
 	fTangentMechanical+=fTangentMechanicalElec;
    return fTangentMechanical;
  }

  //
  // Second Piola-Kirchhoff stress (mechanical)
  //
  inline const dSymMatrixT&
  FSDEMatViscoT::S_IJ()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
	
	dMatrixT stress_temp(3);
	dMatrixT stress_temp2(3);
	double I1 = C(0,0)+C(1,1)+C(2,2);
	double J = C.Det();
	J = sqrt(J);
	
	/* call C function for mechanical part of PK2 stress */
 	mech_pk2_ab_visco(fParams.Pointer(), E.Pointer(), C.Pointer(), J, I1, stress_temp.Pointer()); 
	me_pk2_ab_visco(fParams.Pointer(), E.Pointer(), C.Pointer(), J, stress_temp2.Pointer());
	stress_temp+=stress_temp2;
	
	fStress.FromMatrix(stress_temp);
    return fStress;
  }

  //
  // material electromechanical tangent modulus
  //
  inline const dMatrixT&
  FSDEMatViscoT::E_IJK()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
	double J = C.Det();
	J = sqrt(J);

	/* call C function for electromechanical tangent modulus */
 	me_mixedmodulus_ab_visco(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), J, fTangentElectromechanical.Pointer()); 
 
    return fTangentElectromechanical;

  }

  //
  // material electric tangent modulus
  //
  inline const dMatrixT&
  FSDEMatViscoT::B_IJ()
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
  FSDEMatViscoT::D_I()
  {
  	const dMatrixT& C = RightCauchyGreenDeformation();
  	const dArrayT& E = ElectricField();
  	
  	double J = C.Det();
  	J = sqrt(J);
  
	/* call C function for electric stress (i.e. electric displacement D_{I}) */
 	elec_pk2_visco(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), J, fElectricDisplacement.Pointer()); 
 		
  	return fElectricDisplacement;
  }

  //
  // Electric field
  //
  inline const dArrayT&
  FSDEMatViscoT::E_I()
  {
    return fElectricField;
  }

  //
  // spatial tangent modulus
  //
  inline const dMatrixT&
  FSDEMatViscoT::c_ijkl()
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
  FSDEMatViscoT::s_ij()
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
  inline double FSDEMatViscoT::Pressure() const
  {

    return 0.0;

  }

}	// namespace Tahoe