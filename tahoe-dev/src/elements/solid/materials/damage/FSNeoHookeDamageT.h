//
// $Id: FSNeoHookeDamageT.h,v 1.2 2009-04-02 00:51:15 amota Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.1  2008/12/12 18:59:06  amota
// Initial sources.
//
//

#if !defined(_FSNeoHookeDamageT_)
#define _FSNeoHookeDamageT_

#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "FSIsotropicMatT.h"
#include "FSMatSupportT.h"

namespace Tahoe {

  //
  // Neo-Hookean elastic with simple damage
  //
  class FSNeoHookeDamageT: public FSIsotropicMatT {

    //
    // methods
    //

  public:

    //
    // constructors
    //
    FSNeoHookeDamageT();

    //
    // \name Interface to Tahoe. Member functions inherited from Tahoe
    // to serve as interface with other parts of the code.
    //

    //
    // @{
    //

    //
    // set parameters
    //
    void DefineParameters(ParameterListT& list) const;
    void TakeParameterList(const ParameterListT& list);

    //
    // Interface required by Tahoe
    //
    double StrainEnergyDensity();

    //
    // material mechanical tangent modulus
    //
    virtual const dMatrixT& C_IJKL();

    //
    // Second Piola-Kirchhoff stress
    //
    virtual const dSymMatrixT& S_IJ();
    virtual const dSymMatrixT& S_IJ(const dSymMatrixT& C, const dArrayT& iv);

    //
    // spatial mechanical tangent modulus
    //
    virtual const dMatrixT& c_ijkl();

    //
    // Cauchy stress
    //
    virtual const dSymMatrixT& s_ij();

    //
    // pressure associated with the last computed stress
    //
    double Pressure() const;

    //
    // Auxiliary methods that compute the damage variable in terms of the
    // discontinuous damage.
    // \zeta(\alpha) = \zeta_\inf [1 - \exp(\alpha / \iota)]
    // \zeta      : damage variable
    // \alpha     : discontinuous damage
    // \zeta_\inf : maximum damage
    // \iota      : damage saturation parameter
    //
    double Damage(double alpha) const;

    double DamageDerivative(double alpha) const;

    //
    // @}
    //

    //
    //
    //
    bool HasHistory() const;
    void UpdateHistory(iArrayT& iv, dArrayT& dv, int nip, int ip);
    void ResetHistory(iArrayT& iv, dArrayT& dv, int nip, int ip);

    // apply pre-conditions at the current time step
    void BeginStep(iArrayT& iv, dArrayT& dv, int nip, int ip);

    // finalize the current time step
    void EndStep(iArrayT& iv, dArrayT& dv, int nip, int ip);

    virtual iArrayT InitialIntegerData() const;
    virtual dArrayT InitialDoubleData() const;

    virtual void ComputeOutput(dArrayT& output);

    virtual int NumOutputVariables() const;
    virtual void OutputLabels(ArrayT<StringT>& labels) const;

  protected:

  private:

    void Initialize();

    //
    // Helmholtz free energy density
    //
    double EnergyDensity(const dSymMatrixT& C, double damage) const;

    double EnergyDensityEffective(const dSymMatrixT& C) const;

    double EnergyDensityMechanical(const dSymMatrixT& C) const;

    double EnergyDensityElasticVol(const dSymMatrixT& C) const;
    double EnergyDensityElasticDev(const dSymMatrixT& C) const;

    //
    // 2nd Piola-Kirchhoff stress measures
    //
    const dSymMatrixT Stress(const dSymMatrixT& C, double damage) const;

    const dSymMatrixT StressEffective(const dSymMatrixT& C) const;
    const dSymMatrixT StressMechanical(const dSymMatrixT& C) const;

    const dSymMatrixT StressElasticVol(const dSymMatrixT& C) const;
    const dSymMatrixT StressElasticDev(const dSymMatrixT& C) const;

    //
    // Tangent moduli
    //
    const dMatrixT Tangent(const dSymMatrixT& C, double discontinuousDamage,
        bool damageIncreasing) const;
    const dMatrixT TangentEffective(const dSymMatrixT& C) const;

    const dMatrixT TangentMechanical(const dSymMatrixT& C) const;
    const dMatrixT TangentMechanicalElasticVol(const dSymMatrixT& C) const;
    const dMatrixT TangentMechanicalElasticDev(const dSymMatrixT& C) const;

    const dSymMatrixT RightCauchyGreenDeformation();

    double Damage();
    double DiscontinuousDamage();
    bool DamageIncreasing() const;

    //
    // data
    //

  public:

    static const char* Name;

  protected:

    // number of internal variables
    static const int kNumInternal;
    static const int kNumLabels;
    static const char* kLabels[];

    // indices to internal variables
    enum VariablesT {
      kPrevAlpha = 0, kAlpha = 1
    };

  private:

    dSymMatrixT fStress;
    dMatrixT fTangent;

    // parameters of damage function
    double fMaximumDamage;
    double fDamageSaturation;
    bool fDamageIncreasing;

  };

} // namespace Tahoe

#include "FSNeoHookeDamageT.i.h"

#endif // _FSNeoHookeDamageT_
