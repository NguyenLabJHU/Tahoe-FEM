#include "ExceptionT.h"
#include "FSDEMatQ1P0T.h"
#include "FSDE_incQ1P0.h"

namespace Tahoe {

  //
  //
  //
  static const char DE[] = "Dielectric_Elastomer_Q1P0";
  const char* FSDEMatQ1P0T::Name = DE;
  const int kNSD       = 3;
  const int kNumDOF    = 3;
  const int kStressDim =dSymMatrixT::NumValues(kNSD);
  //
  //
  //
  void FSDEMatQ1P0T::Initialize()
  {
    fMu = 0.0;
    fElectricPermittivity = 0.0;
    fNrig = 0.0;
    fLambda = 0.0;
    fYoung = 0.0;
    fPoisson = 0.0;
  }

  //
  //
  //
  void FSDEMatQ1P0T::DefineParameters(ParameterListT& list) const
  {
	FSSolidMatT::DefineParameters(list);
	
	list.AddParameter(fMu, "mu");
	list.AddParameter(fElectricPermittivity, "epsilon");
 	list.AddParameter(fNrig, "Nrig");
 	list.AddParameter(fLambda, "lambda");
 	list.AddParameter(fYoung, "Young_Modulus");
 	list.AddParameter(fPoisson, "Poisson");
  }

  //
  //
  //
  void FSDEMatQ1P0T::TakeParameterList(const ParameterListT& list)
  {
	FSSolidMatT::TakeParameterList(list);
	fMu = list.GetParameter("mu");
	fElectricPermittivity = list.GetParameter("epsilon");
 	fNrig = list.GetParameter("Nrig");
 	fLambda = list.GetParameter("lambda");
 	fYoung = list.GetParameter("Young_Modulus");
 	fPoisson = list.GetParameter("Poisson");

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
  void FSDEMatQ1P0T::DefineSubs(SubListT& sub_list) const
  {
	FSSolidMatT::DefineSubs(sub_list);
    return;
  }



} //namespace Tahoe
