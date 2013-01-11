#include "ExceptionT.h"
#include "FSDEMatQ1P0SurfaceT.h"
#include "FSDE_incQ1P0Surface.h"

namespace Tahoe {

  FSDEMatQ1P0SurfaceT::FSDEMatQ1P0SurfaceT() :
    ParameterInterfaceT("Dielectric_Elastomer_Q1P0Surface"),
        fFSDEMatSupportQ1P0Surface(0)
  {
    SetName(FSDEMatQ1P0SurfaceT::Name);
  }

  //
  static const char DE[] = "Dielectric_Elastomer_Q1P0Surface";
  const char* FSDEMatQ1P0SurfaceT::Name = DE;
  const int kNSD       = 3;
  const int kNumDOF    = 3;
  const int kStressDim =dSymMatrixT::NumValues(kNSD);

  void FSDEMatQ1P0SurfaceT::Initialize()
  {
    fMu = 0.0;
    fElectricPermittivity = 0.0;
    fNrig = 0.0;
    fLambda = 0.0;  
    
    /* Isotropic surface stuff */
    fE = 0.0;
    fNu = 0.0;
  }

  //
  void FSDEMatQ1P0SurfaceT::DefineParameters(ParameterListT& list) const
  {
	FSSolidMatT::DefineParameters(list);
	
	list.AddParameter(fMu, "mu");
	list.AddParameter(fElectricPermittivity, "epsilon");
 	list.AddParameter(fNrig, "Nrig");
 	list.AddParameter(fLambda, "lambda");
 	
 	/* Isotropic surface stuff */
 	list.AddParameter(fNu, "Poisson");
	list.AddParameter(fE, "Young_Modulus");
   }

  //
  void FSDEMatQ1P0SurfaceT::TakeParameterList(const ParameterListT& list)
  {
	FSSolidMatT::TakeParameterList(list);
	
	fMu = list.GetParameter("mu");
	fElectricPermittivity = list.GetParameter("epsilon");
 	fNrig = list.GetParameter("Nrig");
 	fLambda = list.GetParameter("lambda");

 	/* Isotropic surface stuff */
	fNu = list.GetParameter("Poisson");
	fE = list.GetParameter("Young_Modulus");

	/* write into vector to pass to C code for stress/modulus calculations */
	fParams.Dimension(4);
	fParams[0] = fMu;
	fParams[1] = fLambda;
 	fParams[2] = fElectricPermittivity;
 	fParams[3] = fNrig;
	
	/* dimension work space */
	fTangentMechanical.Dimension(kStressDim);
	fTangentMechanicalElec.Dimension(kStressDim);
	fStress.Dimension(kNumDOF);
	fTangentElectrical.Dimension(kNumDOF);
	fTangentElectromechanical.Dimension(kStressDim, kNumDOF);
	fTangentElectromechanicalSpatial.Dimension(kStressDim, kNumDOF);	
	fElectricDisplacement.Dimension(kNumDOF);
	fElectricField.Dimension(kNumDOF);
  }

  // information about subordinate parameter lists
  void FSDEMatQ1P0SurfaceT::DefineSubs(SubListT& sub_list) const
  {
    /* Inherited */
	FSSolidMatT::DefineSubs(sub_list);
	
    return;
  }

  // Set electrical permittivity
   void FSDEMatQ1P0SurfaceT::SetElectricPermittivity(double epsilon)
  {	
    fElectricPermittivity = epsilon;
  }

  // Get electrical permittivity
   double FSDEMatQ1P0SurfaceT::GetElectricPermittivity() const
  {
    return fElectricPermittivity;
  }

  //
   void FSDEMatQ1P0SurfaceT::SetFSDEMatSupportQ1P0Surface(
      const FSDEMatSupportQ1P0SurfaceT* support)
  {
    fFSDEMatSupportQ1P0Surface = support;
  }

  //
   const dArrayT FSDEMatQ1P0SurfaceT::ElectricField()
  {
    fElectricField = fFSDEMatSupportQ1P0Surface->ElectricField();
    return fElectricField;
  }

  //
   const dArrayT FSDEMatQ1P0SurfaceT::ElectricField(int ip)
  {
    fElectricField = fFSDEMatSupportQ1P0Surface->ElectricField(ip);
    return fElectricField;
  }

  //
   const dMatrixT FSDEMatQ1P0SurfaceT::RightCauchyGreenDeformation()
  {
    const dMatrixT& F = F_mechanical();
    dMatrixT FTF(3);
    FTF.MultATB(F, F);

    return FTF;
  }

  // material energy density
   double FSDEMatQ1P0SurfaceT::StrainEnergyDensity()
  {
	  return 0.0;
  }

  // material mechanical tangent modulus
   const dMatrixT&
  FSDEMatQ1P0SurfaceT::C_IJKL()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
    const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();

	double I1 = C(0,0)+C(1,1)+C(2,2);
	double J = C.Det();
	J = sqrt(J);
	
	/* call C function for mechanical part of tangent modulus */
 	mech_tanmod_q1p0Surface(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, I1, fTangentMechanical.Pointer()); 
 	me_tanmod_q1p0Surface(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, fTangentMechanicalElec.Pointer());
 	fTangentMechanical+=fTangentMechanicalElec;
    return fTangentMechanical;
  }

  // Second Piola-Kirchhoff stress (mechanical)
   const dSymMatrixT&
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
 	mech_pk2_q1p0Surface(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, I1, stress_temp.Pointer()); 
	me_pk2_q1p0Surface(fParams.Pointer(), E.Pointer(), C.Pointer(), F.Pointer(), J, stress_temp2.Pointer());
	stress_temp+=stress_temp2;
	
	fStress.FromMatrix(stress_temp);
    return fStress;
  }

  // material electromechanical tangent modulus
   const dMatrixT&
  FSDEMatQ1P0SurfaceT::E_IJK()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();	
	double J = C.Det();
	J = sqrt(J);

	/* call C function for electromechanical tangent modulus */
 	me_mixedmodulus_q1p0Surface(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), F.Pointer(), J, fTangentElectromechanical.Pointer()); 
 
    return fTangentElectromechanical;
  }

  // spatial electromechanical tangent modulus
   const dMatrixT&
  FSDEMatQ1P0SurfaceT::e_ijk()
  {
    const dMatrixT& C = RightCauchyGreenDeformation();
	const dArrayT& E = ElectricField();
	double J = C.Det();
	J = sqrt(J);  
    const dMatrixT& F = F_mechanical();
	
	/* call C function for (spatial) electromechanical tangent modulus */
 	me_mixedmodulus_q1p0spatialSurface(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), F.Pointer(), J, fTangentElectromechanicalSpatial.Pointer()); 
 	
 	fTangentElectromechanicalSpatial /= J;
    return fTangentElectromechanicalSpatial;
  }

  // material electric tangent modulus
   const dMatrixT&
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
   const dMatrixT&
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
   const dArrayT&
  FSDEMatQ1P0SurfaceT::D_I()
  {
  	const dMatrixT& C = RightCauchyGreenDeformation();
  	const dArrayT& E = ElectricField();
    const dMatrixT& F = F_mechanical();  	
  	
  	double J = C.Det();
  	J = sqrt(J);
  
	/* call C function for electric stress (i.e. electric displacement D_{I}) */
 	elec_pk2_q1p0Surface(fParams.Pointer(), E.Pointer(),  
 		C.Pointer(), F.Pointer(), J, fElectricDisplacement.Pointer()); 
 		
  	return fElectricDisplacement;
  }

  // spatial electric tangent modulus
   const dArrayT&
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
   const dArrayT&
  FSDEMatQ1P0SurfaceT::E_I()
  {
    return fElectricField;
  }

  // spatial tangent modulus
   const dMatrixT&
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
   const dSymMatrixT&
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
   double FSDEMatQ1P0SurfaceT::Pressure() const
  {

    return 0.0;

  }


} //namespace Tahoe
