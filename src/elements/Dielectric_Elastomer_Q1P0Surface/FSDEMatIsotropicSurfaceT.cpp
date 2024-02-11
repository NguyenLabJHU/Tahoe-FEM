#include "ExceptionT.h"
#include "FSDEMatIsotropicSurfaceT.h"

namespace Tahoe {

  FSDEMatIsotropicSurfaceT::FSDEMatIsotropicSurfaceT() :
    ParameterInterfaceT("Dielectric_Elastomer_IsotropicSurface"),
        fLastCall(kNone)
  {
    SetName(FSDEMatIsotropicSurfaceT::Name);
    Initialize();
  }

  //
  static const char DE[] = "Dielectric_Elastomer_IsotropicSurface";
  const char* FSDEMatIsotropicSurfaceT::Name = DE;
  const int kNSD       = 3;
  const int kNumDOF    = 3;
  const int kStressDim =dSymMatrixT::NumValues(kNSD);

  void FSDEMatIsotropicSurfaceT::Initialize()
  {
    fE = 0.0;
    fNu = 0.0;
    
    /* FSDEMatQ1P0SurfaceT stuff */
    fMu = 0.0;
    fElectricPermittivity = 0.0;
    fNrig = 0.0;
    fLambda = 0.0;      
  }

  //
  void FSDEMatIsotropicSurfaceT::DefineParameters(ParameterListT& list) const
  {
	FSSolidMatT::DefineParameters(list);

	/* FSDEMatQ1P0SurfaceT stuff */
	list.AddParameter(fMu, "mu");
	list.AddParameter(fElectricPermittivity, "epsilon");
 	list.AddParameter(fNrig, "Nrig");
 	list.AddParameter(fLambda, "lambda");	

	/* For isotropic "surface" material */
	list.AddParameter(fNu, "Poisson");
	list.AddParameter(fE, "Young_Modulus");
   }

  //
  void FSDEMatIsotropicSurfaceT::TakeParameterList(const ParameterListT& list)
  {
	FSSolidMatT::TakeParameterList(list);

	/* FSDEMatQ1P0SurfaceT stuff */
	fMu = list.GetParameter("mu");
	fElectricPermittivity = list.GetParameter("epsilon");
 	fNrig = list.GetParameter("Nrig");
 	fLambda = list.GetParameter("lambda");
 	
 	/* Isotropic surface stuff */
	fNu = list.GetParameter("Poisson");
	fE = list.GetParameter("Young_Modulus");
	
	/* dimension work space */
	fModulusKStV.Dimension(kStressDim);
	fIsotropicStress.Dimension(kNumDOF);
	fStress1.Dimension(kNumDOF);
	fTM.Dimension(kStressDim);

	/* Calculate (constant) Isotropic modulus */
	SetModulus();
  }

  // information about subordinate parameter lists
  void FSDEMatIsotropicSurfaceT::DefineSubs(SubListT& sub_list) const
  {
    /* Inherited */
	FSSolidMatT::DefineSubs(sub_list);
	
    return;
  }

/* set (material) tangent modulus - from FSKStV.cpp */
void FSDEMatIsotropicSurfaceT::SetModulus()
{	
	double mu     = 0.5*fE/(1.0 + fNu);
	double lambda = 2.0*mu*fNu/(1.0 - 2.0*fNu);
	fModulusKStV = 0.0;
	fModulusKStV(2,2) = fModulusKStV(1,1) = fModulusKStV(0,0) = lambda + 2.0*mu;
	fModulusKStV(1,2) = fModulusKStV(0,1) = fModulusKStV(0,2) = lambda;
	fModulusKStV(5,5) = fModulusKStV(4,4) = fModulusKStV(3,3) = mu;

	/* symmetric */
	fModulusKStV.CopySymmetric();
}

  // material energy density
   double FSDEMatIsotropicSurfaceT::StrainEnergyDensity()
  {
	  return 0.0;
  }

double FSDEMatIsotropicSurfaceT::Pressure() const
{
	return 0.0;
}

  // material mechanical tangent modulus
   const dMatrixT&
  FSDEMatIsotropicSurfaceT::C_IJKL()
  {
  	return fModulusKStV;
  }

  // Second Piola-Kirchhoff stress (mechanical)
   const dSymMatrixT&
  FSDEMatIsotropicSurfaceT::S_IJ()
  {
	dSymMatrixT strain(NumSD());
	
	/* get mechanical part of the deformation gradient */
	const dMatrixT& F_mech = F_mechanical();

	/* strain */
	Compute_E(F_mech, strain);

	/* compute stress */
	fIsotropicStress.A_ijkl_B_kl(fModulusKStV, strain);
	return fIsotropicStress;
  }

  // spatial tangent modulus
   const dMatrixT&
  FSDEMatIsotropicSurfaceT::c_ijkl()
  {
  	const dMatrixT& F = F_mechanical();
  	const double J = F.Det();
  	
  	// prevent aliasing
  	const dMatrixT CIJKL = C_IJKL();
  	
	/* Calculate push forward of modulus */
	fTM.SetToScaled(1.0 / J, PushForward(F, CIJKL));
	return fTM;
  }

  // Cauchy stress
   const dSymMatrixT&
  FSDEMatIsotropicSurfaceT::s_ij()
  {
  	const dMatrixT& F = F_mechanical();
  	const double J = F.Det();
  	
  	// prevent aliasing
  	const dSymMatrixT S = S_IJ();
  	
  	fIsotropicStress.SetToScaled(1.0 / J, PushForward(F, S));
//    fLastCall = kSpatial;
    return fIsotropicStress;
  }

} //namespace Tahoe
