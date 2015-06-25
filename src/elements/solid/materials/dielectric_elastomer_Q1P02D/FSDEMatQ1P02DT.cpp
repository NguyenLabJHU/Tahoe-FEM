#include "ExceptionT.h"
#include "FSDEMatQ1P02DT.h"
#include "FSDE_incQ1P02D.h"

namespace Tahoe {

  //
  //
  //
  static const char DE[] = "Dielectric_Elastomer_Q1P02D";
  const char* FSDEMatQ1P02DT::Name = DE;
  const int kNSD       = 2;
  const int kNumDOF    = 2;
  const int kStressDim =dSymMatrixT::NumValues(kNSD);
  //
  //
  //
  void FSDEMatQ1P02DT::Initialize()
  {
    fMu = 0.0;
    fElectricPermittivity = 0.0;
    fNrig = 0.0;
    fLambda = 0.0;
    fGamma = 0.0;
    fT_0 = 0.0;
  }

  //
  //
  //
  void FSDEMatQ1P02DT::DefineParameters(ParameterListT& list) const
  {
	FSSolidMatT::DefineParameters(list);
	
	list.AddParameter(fMu, "mu");
	list.AddParameter(fElectricPermittivity, "epsilon");
 	list.AddParameter(fNrig, "Nrig");
 	list.AddParameter(fLambda, "lambda");
 	list.AddParameter(fGamma, "gamma");
 	list.AddParameter(fT_0, "t_0");
  }

  //
  //
  //
  void FSDEMatQ1P02DT::TakeParameterList(const ParameterListT& list)
  {
	FSSolidMatT::TakeParameterList(list);
	fMu = list.GetParameter("mu");
	fElectricPermittivity = list.GetParameter("epsilon");
 	fNrig = list.GetParameter("Nrig");
 	fLambda = list.GetParameter("lambda");
 	fGamma = list.GetParameter("gamma");
 	fT_0 = list.GetParameter("t_0");

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
  void FSDEMatQ1P02DT::DefineSubs(SubListT& sub_list) const
  {
	FSSolidMatT::DefineSubs(sub_list);
    return;
  }



} //namespace Tahoe
