
#include "ExceptionT.h"
#include "FSDEMatT.h"
#include "FSDE_inc.h"

namespace Tahoe {

  //
  //
  //
  static const char DE[] = "Dielectric_Elastomer";
  const char* FSDEMatT::Name = DE;
  const int kNSD       = 3;
  const int kNumDOF    = 3;
  const int kStressDim =dSymMatrixT::NumValues(kNSD);
  //
  //
  //
  void FSDEMatT::Initialize()
  {
    fMu = 0.0;
//    fKappa = 0.0;
    fElectricPermittivity = 0.0;
    fNrig = 0.0;
    fLambda = 0.0;
//     fJm = 0.0;
	fTs = 0.0;
  }

  //
  //
  //
  void FSDEMatT::DefineParameters(ParameterListT& list) const
  {
    NL_E_MatT::DefineParameters(list);

	list.AddParameter(fMu, "mu");
//	list.AddParameter(fKappa, "kappa");
	list.AddParameter(fElectricPermittivity, "epsilon");
 	list.AddParameter(fNrig, "Nrig");
 	list.AddParameter(fLambda, "lambda");
 	list.AddParameter(fTs, "TimeStep");
// 	list.AddParameter(fJm, "Jm");
  }

  //
  //
  //
  void FSDEMatT::TakeParameterList(const ParameterListT& list)
  {
//  	cout << "FSDEMatT::TakeParameterList" << endl;
    NL_E_MatT::TakeParameterList(list);

	fMu = list.GetParameter("mu");
//	fKappa = list.GetParameter("kappa");
	fElectricPermittivity = list.GetParameter("epsilon");
// 	fJm = list.GetParameter("Jm");
 	fNrig = list.GetParameter("Nrig");
 	fLambda = list.GetParameter("lambda");
	fTs = list.GetParameter("TimeStep");

	/* write into vector to pass to C code for stress/modulus calculations */
	fParams.Dimension(3);
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
	fElectricDisplacement.Dimension(kNumDOF);
	fElectricField.Dimension(kNumDOF);
  }

  //
  // information about subordinate parameter lists
  //
  void FSDEMatT::DefineSubs(SubListT& sub_list) const
  {
    NL_E_MatT::DefineSubs(sub_list);
    return;
  }



} //namespace Tahoe
