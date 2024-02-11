//
// $Id: FSNeoHookePZLinT.h,v 1.2 2008-12-12 18:58:15 amota Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.1  2008/09/03 18:40:50  beichuan
// Piezoelectricity
//
// Revision 1.2  2008/07/14 17:37:44  lxmota
// Various corrections related to initialization.
//
// Revision 1.1  2008/06/16 18:10:49  lxmota
// Piezoelectric material. Initial sources.
//
//

#if !defined(_FSNeoHookePZLinT_)
#define _FSNeoHookePZLinT_

#include <cassert>

#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "FSIsotropicMatT.h"
#include "FSPZMatSupportT.h"

namespace Tahoe {

  class FSNeoHookePZLinT: public FSIsotropicMatT
  {

    //
    // methods
    //

  public:

    //
    // constructors
    //
    FSNeoHookePZLinT();

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
    // information about subordinate parameter lists
    //
    virtual void DefineSubs(SubListT& sub_list) const;

    //
    // Interface required by Tahoe
    //
    double StrainEnergyDensity();

    //
    // material mechanical tangent modulus
    //
    virtual const dMatrixT& C_IJKL();
    virtual const dMatrixT& C_IJKL_Num();

    //
    // material piezoelectric tangent modulus
    //
    virtual const dMatrixT& H_IJK();

    //
    // material electric tangent modulus
    //
    virtual const dMatrixT& B_IJ();

    //
    // Second Piola-Kirchhoff stress
    //
    virtual const dSymMatrixT& S_IJ();
    virtual const dSymMatrixT& S_IJ_Num();
    virtual const dSymMatrixT& S_IJ(const dSymMatrixT& C);

    //
    // electric displacement
    //
    virtual const dArrayT& D_I();

    //
    // electric field
    //
    virtual const dArrayT& E_I();

    //
    // divergence vector potential
    //
    virtual double DivPhi();

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
    // @}
    //

    //
    // Helmholtz free energy density
    //
    double EnergyDensity(const dSymMatrixT& C, const dArrayT& D) const;

    //
    // 2nd Piola-Kirchhoff stress measures
    //
    const dSymMatrixT Stress(const dSymMatrixT& C, const dArrayT& D) const;

    //
    // Material electric field
    //
    const dArrayT ElectricField(const dSymMatrixT& C, const dArrayT& D) const;

    //
    // accessors and mutators for material constants
    //
    void SetElectricPermittivity(double epsilon);
    double GetElectricPermittivity() const;
    void SetPiezoelectricConstant(int i, int j, double gij);
    double GetPiezoelectricConstant(int i, int j) const;

    void SetPenaltyCoeffifient(double k);
    double GetPenaltyCoefficient() const;

    void SetFSPZMatSupport(const FSPZMatSupportT* support);

    const int ManifoldDim() const;
    const int StrainDim() const;
    const int ElectricalDim() const;

  protected:

    const FSPZMatSupportT* fFSPZMatSupport;

  private:

    void Initialize();

    double EnergyDensityMechanical(const dSymMatrixT& C) const;

    double EnergyDensityElectrical(const dSymMatrixT& C,
				   const dArrayT& D) const;

    double EnergyDensityPiezoelectrical(const dSymMatrixT& C,
					const dArrayT& D) const;

    double EnergyDensityElasticVol(const dSymMatrixT& C) const;
    double EnergyDensityElasticDev(const dSymMatrixT& C) const;

    const dSymMatrixT StressNum(const dSymMatrixT& C, const dArrayT& D) const;

    const dSymMatrixT StressMechanical(const dSymMatrixT& C) const;

    const dSymMatrixT StressElectrical(const dSymMatrixT& C,
				       const dArrayT& D) const;

    const dSymMatrixT StressPiezoelectrical(const dSymMatrixT& C,
					    const dArrayT& D) const;

    const dSymMatrixT StressElasticVol(const dSymMatrixT& C) const;
    const dSymMatrixT StressElasticDev(const dSymMatrixT& C) const;

    const dArrayT ElectricFieldElectrical(const dSymMatrixT& C,
					  const dArrayT& D) const;

    const dArrayT ElectricFieldPiezoelectrical(const dSymMatrixT& C,
					       const dArrayT& D) const;

    //
    // Tangent moduli
    //
    const dMatrixT TangentMechanical(const dSymMatrixT& C,
             const dArrayT& D) const;

    const dMatrixT TangentMechanicalNum(const dSymMatrixT& C,
             const dArrayT& D) const;

    const dMatrixT TangentElectrical(const dSymMatrixT& C,
             const dArrayT& D) const;

    const dMatrixT TangentPiezoelectrical(const dSymMatrixT& C,
            const dArrayT& D) const;

    const dMatrixT TangentMechanicalElasticVol(const dSymMatrixT& C) const;
    const dMatrixT TangentMechanicalElasticDev(const dSymMatrixT& C) const;
    const dMatrixT TangentMechanicalElectrical(const dSymMatrixT& C,
					       const dArrayT& D) const;

    const dSymMatrixT RightCauchyGreenDeformation();
    const dArrayT ElectricDisplacement();
    const dArrayT ElectricDisplacement(int ip);

    double DivergenceVectorPotential();
    double DivergenceVectorPotential(int ip);

    //
    // data
    //

  public:

    static const char* Name;

  protected:

  private:

    static const bool dependsOnC;

    double fElectricPermittivity;
    dMatrixT fPiezoelectricTensor;

    double fEnergyDensity;

    dArrayT fElectricField;
    dArrayT fElectricDisplacement;
    double  fDivergenceVectorPotential;

    dSymMatrixT fStress;
    dMatrixT fTangentMechanical;
    dMatrixT fTangentPiezoelectrical;
    dMatrixT fTangentElectrical;

    // To enforce divergence-free vector potential
    double fPenaltyCoefficient;

  };

} // namespace Tahoe

#include "FSNeoHookePZLinT.i.h"

#endif // _FSNeoHookePZLinT_
