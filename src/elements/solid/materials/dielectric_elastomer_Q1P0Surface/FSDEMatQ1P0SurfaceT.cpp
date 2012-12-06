#include "ExceptionT.h"
#include "FSDEMatQ1P0SurfaceT.h"
#include "FSDE_incQ1P0Surface.h"

namespace Tahoe {

  //
  //
  //
  static const char DE[] = "Dielectric_Elastomer_Q1P0Surface";
  const char* FSDEMatQ1P0SurfaceT::Name = DE;
  const int kNSD       = 3;
  const int kNumDOF    = 3;
  const int kStressDim =dSymMatrixT::NumValues(kNSD);
  //
  //
  //
  void FSDEMatQ1P0SurfaceT::Initialize()
  {
    fMu = 0.0;
    fElectricPermittivity = 0.0;
    fNrig = 0.0;
    fLambda = 0.0;
  }

  //
  //
  //
  void FSDEMatQ1P0SurfaceT::DefineParameters(ParameterListT& list) const
  {
	FSSolidMatT::DefineParameters(list);
	
	list.AddParameter(fMu, "mu");
	list.AddParameter(fElectricPermittivity, "epsilon");
 	list.AddParameter(fNrig, "Nrig");
 	list.AddParameter(fLambda, "lambda");
  }

  //
  //
  //
  void FSDEMatQ1P0SurfaceT::TakeParameterList(const ParameterListT& list)
  {
	FSSolidMatT::TakeParameterList(list);
	fMu = list.GetParameter("mu");
	fElectricPermittivity = list.GetParameter("epsilon");
 	fNrig = list.GetParameter("Nrig");
 	fLambda = list.GetParameter("lambda");

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

  //
  // information about subordinate parameter lists
  //
  void FSDEMatQ1P0SurfaceT::DefineSubs(SubListT& sub_list) const
  {
	FSSolidMatT::DefineSubs(sub_list);
    return;
  }



} //namespace Tahoe
