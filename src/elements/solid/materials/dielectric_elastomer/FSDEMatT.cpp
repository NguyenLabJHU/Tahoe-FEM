
#include "ExceptionT.h"
#include "FSDEMatT.h"

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
  }

  //
  //
  //
  void FSDEMatT::TakeParameterList(const ParameterListT& list)
  {
    NL_E_MatT::TakeParameterList(list);

    fElectricPermittivity = list.GetParameter("epsilon");
	fMu = list.GetParameter("mu");
	fNrig = list.GetParameter("Nrig");

	/* write into vector to pass to C code for stress/modulus calculations */
	fParams.Dimension(3);
	fParams[0] = fElectricPermittivity;
	fParams[1] = fMu;
	fParams[2] = fNrig;
	
	/* dimension work space */
	fTangentMechanical.Dimension(kStressDim);
	fStress.Dimension(kNumDOF);
	stress_temp.Dimension(kNumDOF);
	fTangentElectrical.Dimension(kNumDOF);
	fTangentElectromechanical.Dimension(kStressDim, kNumDOF);
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
