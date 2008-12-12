//
// $Id: FSNeoHookeDamageT.cpp,v 1.1 2008-12-12 18:59:06 amota Exp $
//
// $Log: not supported by cvs2svn $
//

#include "ExceptionT.h"
#include "FSNeoHookeDamageT.h"

namespace Tahoe {

  //
  //
  //
  static const char NPD[] = "Neohookean-damage";
  const char* FSNeoHookeDamageT::Name = NPD;
  const int FSNeoHookeDamageT::kNumInternal = 2;
  const int FSNeoHookeDamageT::kNumLabels = 4;
  const char* FSNeoHookeDamageT::kLabels[] = { "damage", "energy density",
      "effective energy density", "discontinuous damage" };

  //
  //
  //
  void FSNeoHookeDamageT::Initialize()
  {

    fMaximumDamage = 0.0;
    fDamageSaturation = 0.0;

  }

  //
  //
  //
  void FSNeoHookeDamageT::DefineParameters(ParameterListT& list) const
  {

    FSIsotropicMatT::DefineParameters(list);

    list.AddParameter(fMaximumDamage, "zeta_inf");
    list.AddParameter(fDamageSaturation, "iota");

    //
    // set the description
    //
    list.SetDescription("Psi(C)=0.5*mu*(I1bar-3)+0.25*kappa*(J^2-1-2*log(J))");

  }

  //
  //
  //
  void FSNeoHookeDamageT::TakeParameterList(const ParameterListT& list)
  {

    FSIsotropicMatT::TakeParameterList(list);

    fMaximumDamage = list.GetParameter("zeta_inf");
    fDamageSaturation = list.GetParameter("iota");

    //
    // check
    //
    bool maxDamageInRange = (0.0 <= fMaximumDamage) && (fMaximumDamage <= 1.0);
    if (maxDamageInRange == false) {
      ExceptionT::BadInputValue("FSNeoHookeDamageT::TakeParameterList",
          "expecting 0 <= zeta_inf <= 1: %e", fMaximumDamage);
    }

    bool posDefDamageSat = 0.0 < fDamageSaturation;
    if (posDefDamageSat == false) {
      ExceptionT::BadInputValue("FSNeoHookeDamageT::TakeParameterList",
          "expecting positive definite iota: %e", fDamageSaturation);
    }

  }

} //namespace Tahoe
