#include "ExceptionT.h"
#include "FSDEMatQ1P0ElastocapillaryT.h"

namespace Tahoe {

  FSDEMatQ1P0ElastocapillaryT::FSDEMatQ1P0ElastocapillaryT():
    ParameterInterfaceT("Finite Strain Dielectric Elastomer Elastocapillary Q1P0")
  {
    SetName(FSDEMatQ1P0ElastocapillaryT::Name);
    Initialize();	
  }
  //
  static const char DE[] = "Dielectric_Elastomer_Q1P0_Elastocapillary";
  const char* FSDEMatQ1P0ElastocapillaryT::Name = DE;
  const int kNSD       = 3;
  const int kNumDOF    = 3;
  const int kStressDim =dSymMatrixT::NumValues(kNSD);

  //
  void FSDEMatQ1P0ElastocapillaryT::Initialize()
  {
		fSurfTension = 0.0;
  }

  //
  void FSDEMatQ1P0ElastocapillaryT::DefineParameters(ParameterListT& list) const
  {
	FSSolidMatT::DefineParameters(list);

  }

  //
  void FSDEMatQ1P0ElastocapillaryT::TakeParameterList(const ParameterListT& list)
  {
	FSSolidMatT::TakeParameterList(list);
	fSurfTension = list.GetParameter("gamma");
	cout << "fSurfTension = " << fSurfTension << endl;
	
	/* dimension work space for stress and stiffness */
	fCapillaryStress.Dimension(kNumDOF);
	fCapillaryTangent.Dimension(kStressDim);
  }

  // information about subordinate parameter lists
  void FSDEMatQ1P0ElastocapillaryT::DefineSubs(SubListT& sub_list) const
  {
	FSSolidMatT::DefineSubs(sub_list);
    return;
  }

const dMatrixT& FSDEMatQ1P0ElastocapillaryT::c_ijkl()
{


	return fCapillaryTangent;
}

const dSymMatrixT& FSDEMatQ1P0ElastocapillaryT::s_ij()
{


	return fCapillaryStress;
}

} //namespace Tahoe
