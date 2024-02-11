#include "ExceptionT.h"
#include "FSDEMat2DT.h"
#include "FSDE_inc_2D.h"

namespace Tahoe {

  //
  //
  //
  static const char DE[] = "Dielectric_Elastomer_2D";
  const char* FSDEMat2DT::Name = DE;
  const int kNSD       = 2;
  const int kNumDOF    = 2;
  const int kStressDim =dSymMatrixT::NumValues(kNSD);
  //
  //
  //
  void FSDEMat2DT::Initialize()
  {
    fElectricPermittivity = 0.0;
    fMu = 0.0;
    fNrig = 0.0;
    fLambda = 0.0;
  }

  //
  //
  //
  void FSDEMat2DT::DefineParameters(ParameterListT& list) const
  {
    NL_E_MatT::DefineParameters(list);

	list.AddParameter(fElectricPermittivity, "epsilon");
	list.AddParameter(fMu, "mu");
	list.AddParameter(fNrig, "Nrig");
	list.AddParameter(fLambda, "lambda");
  }

  //
  //
  //
  void FSDEMat2DT::TakeParameterList(const ParameterListT& list)
  {
//  	cout << "FSDEMat2DT::TakeParameterList" << endl;
    NL_E_MatT::TakeParameterList(list);

    fElectricPermittivity = list.GetParameter("epsilon");
	fMu = list.GetParameter("mu");
	fNrig = list.GetParameter("Nrig");
	fLambda = list.GetParameter("lambda");

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
  void FSDEMat2DT::DefineSubs(SubListT& sub_list) const
  {
    NL_E_MatT::DefineSubs(sub_list);
    return;
  }



} //namespace Tahoe
