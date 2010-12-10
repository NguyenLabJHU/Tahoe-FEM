
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
    fElectricPermittivity = 0.0;
    fMu = 0.0;
    fNrig = 0.0;
    fLambda = 0.0;
    fKappa = 0.0;
    fJm = 0.0;
  }

  //
  //
  //
  void FSDEMatT::DefineParameters(ParameterListT& list) const
  {
    NL_E_MatT::DefineParameters(list);

	list.AddParameter(fElectricPermittivity, "epsilon");
	list.AddParameter(fMu, "mu");
	list.AddParameter(fNrig, "Nrig");
	list.AddParameter(fLambda, "lambda");
	list.AddParameter(fKappa, "kappa");
	list.AddParameter(fJm, "Jm");
  }

  //
  //
  //
  void FSDEMatT::TakeParameterList(const ParameterListT& list)
  {
//  	cout << "FSDEMatT::TakeParameterList" << endl;
    NL_E_MatT::TakeParameterList(list);

    fElectricPermittivity = list.GetParameter("epsilon");
	fMu = list.GetParameter("mu");
	fNrig = list.GetParameter("Nrig");
	fLambda = list.GetParameter("lambda");
	fKappa = list.GetParameter("kappa");
	fJm = list.GetParameter("Jm");

	/* write into vector to pass to C code for stress/modulus calculations */
	fParams.Dimension(6);
	fParams[0] = fElectricPermittivity;
	fParams[1] = fMu;
	fParams[2] = fNrig;
	fParams[3] = fLambda;
	fParams[4] = fKappa;
	fParams[5] = fJm;
	
	/* dimension work space */
	fTangentMechanical.Dimension(kStressDim);
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
