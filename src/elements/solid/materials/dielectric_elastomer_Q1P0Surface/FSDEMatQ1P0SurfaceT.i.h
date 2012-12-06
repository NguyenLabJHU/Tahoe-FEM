#include "FSDE_incQ1P0Surface.h"

namespace Tahoe {

  inline FSDEMatQ1P0SurfaceT::FSDEMatQ1P0SurfaceT() :
    ParameterInterfaceT("Finite Strain Dielectric Elastomer Q1P0Surface"),
        fFSDEMatSupportQ1P0Surface(0)
  {
    SetName(FSDEMatQ1P0SurfaceT::Name);
    Initialize();
  }

  // Set electrical permittivity
  inline void FSDEMatQ1P0SurfaceT::SetElectricPermittivity(double epsilon)
  {	
    fElectricPermittivity = epsilon;
  }

  // Get electrical permittivity
  inline double FSDEMatQ1P0SurfaceT::GetElectricPermittivity() const
  {
    return fElectricPermittivity;
  }

  //
  inline void FSDEMatQ1P0SurfaceT::SetFSDEMatSupportQ1P0Surface(
      const FSDEMatSupportQ1P0SurfaceT* support)
  {
    fFSDEMatSupportQ1P0Surface = support;
  }

  //
  inline const dArrayT FSDEMatQ1P0SurfaceT::ElectricField()
  {
    fElectricField = fFSDEMatSupportQ1P0Surface->ElectricField();
    return fElectricField;
  }

  //
  inline const dArrayT FSDEMatQ1P0SurfaceT::ElectricField(int ip)
  {
    fElectricField = fFSDEMatSupportQ1P0Surface->ElectricField(ip);
    return fElectricField;
  }

  //
  inline const dMatrixT FSDEMatQ1P0SurfaceT::RightCauchyGreenDeformation()
  {
    const dMatrixT& F = F_mechanical();
    dMatrixT FTF(3);
    FTF.MultATB(F, F);

    return FTF;
  }

  // material energy density
  inline double FSDEMatQ1P0SurfaceT::StrainEnergyDensity()
  {
	  return 0.0;
  }

  // material mechanical tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P0SurfaceT::C_IJKL()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
    const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();

	double I1 = C(0,0)+C(1,1)+C(2,2);
	double J = C.Det();
	J = sqrt(J);
	
	/* call C function for mechanical part of tangent modulus */
// 	mech_tanmod_Q1P0Surface(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, I1, fTangentMechanical.Pointer()); 
// 	me_tanmod_Q1P0Surface(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, fTangentMechanicalElec.Pointer());
 	fTangentMechanical+=fTangentMechanicalElec;
    return fTangentMechanical;
  }

  // Second Piola-Kirchhoff stress (mechanical)
  inline const dSymMatrixT&
  FSDEMatQ1P0SurfaceT::S_IJ()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();
	
	dMatrixT stress_temp(3);
	dMatrixT stress_temp2(3);
	double I1 = C(0,0)+C(1,1)+C(2,2);
	double J = C.Det();
	J = sqrt(J);
	
	/* call C function for mechanical part of PK2 stress */
// 	mech_pk2_Q1P0Surface(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, I1, stress_temp.Pointer()); 
//	me_pk2_Q1P0Surface(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, stress_temp2.Pointer());
	stress_temp+=stress_temp2;
	
	fStress.FromMatrix(stress_temp);
    return fStress;
  }

  // material electromechanical tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P0SurfaceT::E_IJK()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();	
	double J = C.Det();
	J = sqrt(J);

	/* call C function for electromechanical tangent modulus */
// 	me_mixedmodulus_Q1P0Surface(fParams.Pointer(), E.Pointer(),  
// 		C.Pointer(), F.Pointer(), J, fTangentElectromechanical.Pointer()); 
 
    return fTangentElectromechanical;
  }

  // spatial electromechanical tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P0SurfaceT::e_ijk()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
	double J = C.Det();
	J = sqrt(J);  
    const dMatrixT& F = F_mechanical();
	
	/* call C function for (spatial) electromechanical tangent modulus */
// 	me_mixedmodulus_Q1P0spatialSurface(fParams.Pointer(), E.Pointer(),  
// 		C.Pointer(), F.Pointer(), J, fTangentElectromechanicalSpatial.Pointer()); 
 	
 	fTangentElectromechanicalSpatial /= J;
    return fTangentElectromechanicalSpatial;
  }

  // material electric tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P0SurfaceT::B_IJ()
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

  // spatial electric tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P0SurfaceT::b_ij()
  {
    const dMatrixT& F = F_mechanical();
    const dMatrixT& C = RightCauchyGreenDeformation();    
    const double J = F.Det();
	
	// repeat B_IJ first, then push forward
	dMatrixT Cinv(3), bij(3);
	Cinv.Inverse(C);
	bij = Cinv;
	bij *= fElectricPermittivity;
	bij *= J;
	
    // prevent aliasing
//    const dMatrixT b = B_IJ();
    fTangentElectrical.MultABCT(F, bij, F);
    fTangentElectrical /= J;
    return fTangentElectrical;
  }

  // Electric displacement 
  inline const dArrayT&
  FSDEMatQ1P0SurfaceT::D_I()
  {
  	const dMatrixT& C = RightCauchyGreenDeformation();
  	const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();  	
  	
  	double J = C.Det();
  	J = sqrt(J);
  
	/* call C function for electric stress (i.e. electric displacement D_{I}) */
// 	elec_pk2_Q1P0Surface(fParams.Pointer(), E.Pointer(),  
// 		C.Pointer(), F.Pointer(), J, fElectricDisplacement.Pointer()); 
 		
  	return fElectricDisplacement;
  }

  // spatial electric tangent modulus
  inline const dArrayT&
  FSDEMatQ1P0SurfaceT::d_i()
  {
    const dMatrixT& F = F_mechanical();
    const double J = F.Det();
	
    // prevent aliasing
    const dArrayT D = D_I();
	F.Multx(D, fElectricDisplacement);
	fElectricDisplacement /= J;
    return fElectricDisplacement;		// need to divide by J
  }

  // Electric field
  inline const dArrayT&
  FSDEMatQ1P0SurfaceT::E_I()
  {
    return fElectricField;
  }

  // spatial tangent modulus
  inline const dMatrixT&
  FSDEMatQ1P0SurfaceT::c_ijkl()
  {
//     const dMatrixT& F = F_mechanical();
//     const double J = F.Det();
// 
//     // prevent aliasing
//     const dMatrixT CIJKL = C_IJKL();
//     fTangentMechanical.SetToScaled(1.0 / J, PushForward(F, CIJKL));

	fTangentMechanical = FSSolidMatT::c_ijkl();

    return fTangentMechanical;

  }

  // Cauchy stress
  inline const dSymMatrixT&
  FSDEMatQ1P0SurfaceT::s_ij()
  {
    const dMatrixT& F = F_mechanical();
    const double J = F.Det();
	
    // prevent aliasing
    const dSymMatrixT S = S_IJ();

    fStress.SetToScaled(1.0 / J, PushForward(F, S));
    return fStress;
  }

  // pressure associated with the last computed stress
  inline double FSDEMatQ1P0SurfaceT::Pressure() const
  {

    return 0.0;

  }

}	// namespace Tahoe