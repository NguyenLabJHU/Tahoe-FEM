
#include "ExceptionT.h"
#include "FSDEMatT.h"

namespace Tahoe {

  //
  //
  //
  static const char DE[] = "Dielectric-elastomer";
  const char* FSDEMatT::Name = DE;

  //
  //
  //
  void FSDEMatT::Initialize()
  {

    fElectricPermittivity = 0.0;

  }

  //
  //
  //
  void FSDEMatT::DefineParameters(ParameterListT& list) const
  {

    FSIsotropicMatT::DefineParameters(list);

    list.AddParameter(fElectricPermittivity, "epsilon");

    //
    // set the description
    //
//    list.SetDescription("Psi(C)=0.5*mu*(I1bar-3)+0.25*kappa*(J^2-1-2*log(J))");

  }

  //
  //
  //
  void FSDEMatT::TakeParameterList(const ParameterListT& list)
  {

    FSIsotropicMatT::TakeParameterList(list);

    fElectricPermittivity = list.GetParameter("epsilon");

    //
    // check
    //
    if (fElectricPermittivity < -kSmall) {
      ExceptionT::BadInputValue("FSDEMatT::TakeParameterList",
          "expecting a non-negative epsilon: %e", fElectricPermittivity);
    }


  }

  //
  // information about subordinate parameter lists
  //
  void FSDEMatT::DefineSubs(SubListT& sub_list) const
  {
    FSIsotropicMatT::DefineSubs(sub_list);
    return;
  }



} //namespace Tahoe
