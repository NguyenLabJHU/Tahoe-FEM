//
// $Id: FSNeoHookeDamageT.i.h,v 1.2 2009-04-02 00:51:15 amota Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.1  2008/12/12 18:59:06  amota
// Initial sources.
//
//

#include <algorithm>

namespace Tahoe {

  inline FSNeoHookeDamageT::FSNeoHookeDamageT() :
    ParameterInterfaceT("Neohookean elastic with simple damage")
  {
    SetName(FSNeoHookeDamageT::Name);
    Initialize();
  }

  //
  //
  //
  inline iArrayT FSNeoHookeDamageT::InitialIntegerData() const
  {

    //empty array
    iArrayT ia;

    return ia;

  }

  //
  //
  //
  inline dArrayT FSNeoHookeDamageT::InitialDoubleData() const
  {
    // internal variables are last and current discontinuous damage,
    // with zero as initial value
    dArrayT da(2);
    da[0] = 0.0;
    da[1] = 0.0;

    return da;
  }

  //
  //
  //
  inline bool FSNeoHookeDamageT::HasHistory() const
  {
    return true;
  }

  //
  //
  //
  inline void FSNeoHookeDamageT::ResetHistory(iArrayT& iv, dArrayT& dv,
      int nip, int ip)
  {

    // Nothing to do. Internal variables not modified by default.

  }

  //
  // apply pre-conditions at the current time step
  //
  inline void FSNeoHookeDamageT::BeginStep(iArrayT& iv, dArrayT& dv, int nip,
      int ip)
  {

    // Nothing to do. Internal variables not modified by default.

  }

  //
  //
  //
  inline void FSNeoHookeDamageT::UpdateHistory(iArrayT& iv, dArrayT& dv,
      int nip, int ip)
  {

    // compute current effective energy density
    dSymMatrixT C = RightCauchyGreenDeformation();
    double W0 = EnergyDensityEffective(C);
    double& discontinuousDamage = dv[kAlpha * nip + ip];

    discontinuousDamage = std::max(discontinuousDamage, W0);

  }

  //
  // finalize the current time step
  //
  inline void FSNeoHookeDamageT::EndStep(iArrayT& iv, dArrayT& dv, int nip,
      int ip)
  {

    double& lastDiscontinuousDamage = dv[kPrevAlpha * nip + ip];
    double& discontinuousDamage = dv[kAlpha * nip + ip];

    lastDiscontinuousDamage = discontinuousDamage;

  }

  //
  //
  //
  inline double FSNeoHookeDamageT::DiscontinuousDamage()
  {

    ElementCardT& element = CurrentElement();

    if (element.IsAllocated() == 0) {
      return 0.0;
    }

    dArrayT& dv = element.DoubleData();
    const int nip = NumIP();
    const int ip = CurrIP();

    dSymMatrixT C = RightCauchyGreenDeformation();
    double W0 = EnergyDensityEffective(C);
    double& discontinuousDamage = dv[kAlpha * nip + ip];

    return std::max(discontinuousDamage, W0);

  }

  //
  //
  //
  inline double FSNeoHookeDamageT::Damage()
  {

    return Damage(DiscontinuousDamage());

  }

  //
  //
  //
  inline bool FSNeoHookeDamageT::DamageIncreasing() const
  {

    ElementCardT& element = CurrentElement();

    if (element.IsAllocated() == 0) {
      return false;
    }

    dArrayT& dv = element.DoubleData();
    const int nip = NumIP();
    const int ip = CurrIP();
    double& discontinuousDamage = dv[kAlpha * nip + ip];
    double& lastDiscontinuousDamage = dv[kPrevAlpha * nip + ip];

    return discontinuousDamage > lastDiscontinuousDamage;

  }

  //
  //
  //
  inline void FSNeoHookeDamageT::ComputeOutput(dArrayT& output)
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();

    output[0] = Damage();
    output[1] = StrainEnergyDensity();
    output[2] = EnergyDensityEffective(C);
    output[3] = DiscontinuousDamage();

  }

  //
  //
  //
  inline int FSNeoHookeDamageT::NumOutputVariables() const
  {
    return kNumLabels;
  }

  //
  //
  //
  inline void FSNeoHookeDamageT::OutputLabels(ArrayT<StringT>& labels) const
  {
    labels.Dimension(kNumLabels);

    for (int i = 0; i < kNumLabels; ++i) {
      labels[i] = kLabels[i];
    }

  }

  //
  //
  //
  inline double FSNeoHookeDamageT::EnergyDensity(const dSymMatrixT& C,
      double damage) const
  {

    const double W0 = EnergyDensityEffective(C);

    return (1.0 - damage) * W0;

  }

  inline double FSNeoHookeDamageT::EnergyDensityEffective(const dSymMatrixT& C) const
  {

    const double Wm = EnergyDensityMechanical(C);

    return Wm;

  }

  //
  //
  //
  inline double FSNeoHookeDamageT::EnergyDensityMechanical(const dSymMatrixT& C) const
  {

    const double Wevol = EnergyDensityElasticVol(C);
    const double Wedev = EnergyDensityElasticDev(C);

    return Wevol + Wedev;

  }

  //
  //
  //
  inline double FSNeoHookeDamageT::EnergyDensityElasticVol(const dSymMatrixT& C) const
  {

    const double Jsqr = C.Det();
    const double twotheta = log(Jsqr);
    const double Wevol = 0.25 * fKappa * (Jsqr - 1.0 - twotheta);

    return Wevol;

  }

  //
  //
  //
  inline double FSNeoHookeDamageT::EnergyDensityElasticDev(const dSymMatrixT& C) const
  {

    const double Jsqr = C.Det();
    const double Jm23 = pow(Jsqr, -1.0 / 3.0);

    const int nsd = NumSD();
    dSymMatrixT Cbar(nsd);
    Cbar.SetToScaled(Jm23, C);

    const double Wedev = 0.5 * fMu * (Cbar.Trace() - 3.0);

    return Wedev;

  }

  //
  //
  //
  inline const dSymMatrixT FSNeoHookeDamageT::Stress(const dSymMatrixT& C,
      double damage) const
  {

    dSymMatrixT S = StressEffective(C);
    S *= (1.0 - damage);

    return S;

  }

  //
  //
  //
  inline const dSymMatrixT FSNeoHookeDamageT::StressEffective(
      const dSymMatrixT& C) const
  {

    const dSymMatrixT Sm = StressMechanical(C);

    return Sm;

  }

  //
  //
  //
  inline const dSymMatrixT FSNeoHookeDamageT::StressMechanical(
      const dSymMatrixT& C) const
  {

    const dSymMatrixT Sevol = StressElasticVol(C);
    const dSymMatrixT Sedev = StressElasticDev(C);

    dSymMatrixT Sm = Sevol;
    Sm += Sedev;

    return Sm;

  }

  //
  //
  //
  inline const dSymMatrixT FSNeoHookeDamageT::StressElasticVol(
      const dSymMatrixT& C) const
  {

    dSymMatrixT Sevol = C;
    Sevol.Inverse();

    const double Jsqr = C.Det();
    Sevol *= (0.5 * fKappa * (Jsqr - 1.0));

    return Sevol;

  }

  //
  //
  //
  inline const dSymMatrixT FSNeoHookeDamageT::StressElasticDev(
      const dSymMatrixT& C) const
  {

    dSymMatrixT Sedev = C;
    Sedev.Inverse();

    Sedev *= (-C.Trace() / 3.0);
    Sedev.PlusIdentity(1.0);

    const double Jsqr = C.Det();
    const double Jm23 = pow(Jsqr, -1.0 / 3.0);
    Sedev *= (fMu * Jm23);

    return Sedev;

  }

  //
  //
  //
  inline const dMatrixT FSNeoHookeDamageT::Tangent(const dSymMatrixT& C,
      double discontinuousDamage, bool damageIncreasing) const
  {

    dMatrixT CC = TangentEffective(C);
    CC *= (1.0 - Damage(discontinuousDamage));

    if (damageIncreasing == true) {

      dSymMatrixT S0 = StressEffective(C);

      const int nsd = NumSD();
      const int ned = dSymMatrixT::NumValues(nsd);
      dMatrixT S0_x_S0(ned);
      S0_x_S0.DyadAB(S0, S0);
      S0_x_S0 *= -DamageDerivative(discontinuousDamage);

      CC += S0_x_S0;

    }

    return CC;

  }

  //
  //
  //
  inline const dMatrixT FSNeoHookeDamageT::TangentEffective(
      const dSymMatrixT& C) const
  {

    const dMatrixT Cm = TangentMechanical(C);

    return Cm;

  }

  //
  //
  //
  inline const dMatrixT FSNeoHookeDamageT::TangentMechanical(
      const dSymMatrixT& C) const
  {

    const dMatrixT Cevol = TangentMechanicalElasticVol(C);
    const dMatrixT Cedev = TangentMechanicalElasticDev(C);

    dMatrixT Cm = Cevol;
    Cm += Cedev;

    return Cm;

  }

  //
  //
  //
  inline const dMatrixT FSNeoHookeDamageT::TangentMechanicalElasticVol(
      const dSymMatrixT& C) const
  {

    dSymMatrixT Cinv = C;
    Cinv.Inverse();

    const int nsd = NumSD();
    const int ned = dSymMatrixT::NumValues(nsd);
    dMatrixT Cinv_x_Cinv(ned);
    Cinv_x_Cinv.DyadAB(Cinv, Cinv);

    dMatrixT Cinv_o_Cinv(ned);
    Cinv_o_Cinv.ReducedI_C(Cinv);

    const double Jsqr = C.Det();
    Cinv_x_Cinv *= Jsqr;
    Cinv_o_Cinv *= (Jsqr - 1.0);
    dMatrixT Cevol = Cinv_x_Cinv;
    Cevol -= Cinv_o_Cinv;
    Cevol *= fKappa;

    return Cevol;

  }

  //
  //
  //
  inline const dMatrixT FSNeoHookeDamageT::TangentMechanicalElasticDev(
      const dSymMatrixT& C) const
  {

    dSymMatrixT Cinv = C;
    Cinv.Inverse();

    const int nsd = NumSD();
    const int ned = dSymMatrixT::NumValues(nsd);

    dMatrixT Cinv_x_Cinv(ned);
    Cinv_x_Cinv.DyadAB(Cinv, Cinv);

    dMatrixT Cinv_o_Cinv(ned);
    Cinv_o_Cinv.ReducedI_C(Cinv);

    dSymMatrixT I = C;
    I.Identity(1.0);
    dMatrixT Cinv_x_I(ned);
    Cinv_x_I.DyadAB(Cinv, I);

    dMatrixT I_x_Cinv(ned);
    I_x_Cinv.DyadAB(I, Cinv);

    Cinv_x_Cinv /= 3.0;
    dMatrixT Cedev = Cinv_x_Cinv;
    Cedev += Cinv_o_Cinv;
    Cedev *= C.Trace();
    Cedev -= Cinv_x_I;
    Cedev -= I_x_Cinv;
    const double Jsqr = C.Det();
    const double Jm23 = pow(Jsqr, -1.0 / 3.0);
    Cedev *= (2.0 * fMu * Jm23 / 3.0);

    return Cedev;

  }

  //
  //
  //
  inline const dSymMatrixT FSNeoHookeDamageT::RightCauchyGreenDeformation()
  {

    const dMatrixT F = F_mechanical();
    dMatrixT FTF(NumSD());
    FTF.MultATB(F, F);
    dSymMatrixT C(NumSD());
    C.Symmetrize(FTF);

    return C;

  }

  //
  // material energy density
  //
  inline double FSNeoHookeDamageT::StrainEnergyDensity()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const double damage = Damage();
    double W = EnergyDensity(C, damage);

    return W;

  }

  //
  // material mechanical tangent modulus
  //
  inline const dMatrixT&
  FSNeoHookeDamageT::C_IJKL()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const double damage = Damage();
    const bool damageIncreasing = DamageIncreasing();
    fTangent = Tangent(C, damage, damageIncreasing);

    return fTangent;

  }

  //
  // Second Piola-Kirchhoff stress
  //
  inline const dSymMatrixT&
  FSNeoHookeDamageT::S_IJ()
  {

    const dSymMatrixT C = RightCauchyGreenDeformation();
    const double damage = Damage();
    fStress = Stress(C, damage);

    return fStress;

  }

  //
  // Second Piola-Kirchhoff stress
  //
  inline const dSymMatrixT&
  FSNeoHookeDamageT::S_IJ(const dSymMatrixT& C, const dArrayT& iv)
  {

    const double alpha = iv[kAlpha];
    const double damage = Damage(alpha);

    fStress = Stress(C, damage);

    return fStress;

  }

  //
  // spatial tangent modulus
  //
  inline const dMatrixT&
  FSNeoHookeDamageT::c_ijkl()
  {

    const dMatrixT F = F_mechanical();
    const double J = F.Det();

    // prevent aliasing
    const dMatrixT CIJKL = C_IJKL();

    fTangent.SetToScaled(1.0 / J, PushForward(F, CIJKL));

    return fTangent;

  }

  //
  // Cauchy stress
  //
  inline const dSymMatrixT&
  FSNeoHookeDamageT::s_ij()
  {

    const dMatrixT F = F_mechanical();
    const double J = F.Det();

    // prevent aliasing
    const dSymMatrixT S = S_IJ();

    fStress.SetToScaled(1.0 / J, PushForward(F, S));

    return fStress;

  }

  //
  // pressure associated with the last computed stress
  //
  inline double FSNeoHookeDamageT::Pressure() const
  {

    return fStress.Trace() / 3.0;

  }

  //
  // Auxiliary methods that compute the damage variable in terms of the
  // discontinuous damage.
  // \zeta(\alpha) = \zeta_\inf [1 - \exp(\alpha / \iota)]
  // \zeta      : damage variable
  // \alpha     : discontinuous damage
  // \zeta_\inf : maximum damage
  // \iota      : damage saturation parameter
  //
  inline double FSNeoHookeDamageT::Damage(double alpha) const
  {
    return fMaximumDamage * (1.0 - exp(-alpha / fDamageSaturation));
  }

  //
  //
  //
  inline double FSNeoHookeDamageT::DamageDerivative(double alpha) const
  {
    return fMaximumDamage / fDamageSaturation * exp(-alpha / fDamageSaturation);
  }

} //namespace Tahoe
