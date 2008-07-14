//
// $Id: FSNeoHookePZLinT.h,v 1.2 2008-07-14 17:37:44 lxmota Exp $
//
// $Log: not supported by cvs2svn $
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
    const dMatrixT& C_ijkl();

    //
    // material piezoelectric tangent modulus
    //
    const dMatrixT& H_ijk();

    //
    // material electric tangent modulus
    //
    const dMatrixT& B_ij();

    //
    // Second Piola-Kirchhoff stress
    //
    const dSymMatrixT& S_ij();

    //
    // electric field
    //
    const dArrayT& E_i();

    //
    // spatial mechanical tangent modulus
    //
    const dMatrixT& c_ijkl();

    //
    // Cauchy stress
    //
    const dSymMatrixT& s_ij();

    //
    // pressure associated with the last computed stress
    //
    double Pressure() const;

    //
    // @}
    //

    //
    //
    //
    virtual bool Need_D() const;
    virtual bool Need_D_last() const;

    //
    // Helmholtz free energy density
    //
    double energyDensity(const dSymMatrixT& C, const dArrayT& D) const;

    //
    // 2nd Piola-Kirchhoff stress measures
    //
    const dSymMatrixT stress2PK(const dSymMatrixT& C, const dArrayT& D) const;

    //
    // Material electric field
    //
    const dArrayT electricField(const dSymMatrixT& C, const dArrayT& D) const;

    //
    // Tangent moduli
    //
    const dMatrixT tangentMechanical(const dSymMatrixT& C,
				     const dArrayT& D) const;

    const dMatrixT tangentElectrical(const dSymMatrixT& C,
				     const dArrayT& D) const;

    const dMatrixT tangentPiezoelectrical(const dSymMatrixT& C,
					  const dArrayT& D) const;

    //
    // accesors and mutators for material constants
    //
    void setShearModulus(double mu);
    double getShearModulus() const;
    void setBulkModulus(double kappa);
    double getBulkModulus() const;
    void setElectricPermittivity(double epsilon);
    double getElectricPermittivity() const;
    void setPiezoelectricConstant(int i, int j, double gij);
    double getPiezoelectricConstant(int i, int j) const;

    void setFSPZMatSupport(const FSPZMatSupportT* support);

    const int ManifoldDim() const;
    const int StrainDim() const;
    const int ElectricalDim() const;

  protected:

    const FSPZMatSupportT* fFSPZMatSupport;

  private:

    void initialize();

    double energyDensityMechanical(const dSymMatrixT& C) const;

    double energyDensityElectrical(const dSymMatrixT& C,
				   const dArrayT& D) const;

    double energyDensityPiezoelectrical(const dSymMatrixT& C,
					const dArrayT& D) const;

    double energyDensityElasticVol(const dSymMatrixT& C) const;
    double energyDensityElasticDev(const dSymMatrixT& C) const;

    const dSymMatrixT stressMechanical(const dSymMatrixT& C) const;

    const dSymMatrixT stressElectrical(const dSymMatrixT& C,
				       const dArrayT& D) const;

    const dSymMatrixT stressPiezoelectrical(const dSymMatrixT& C,
					    const dArrayT& D) const;

    const dSymMatrixT stressElasticVol(const dSymMatrixT& C) const;
    const dSymMatrixT stressElasticDev(const dSymMatrixT& C) const;

    const dArrayT electricFieldElectrical(const dSymMatrixT& C,
					  const dArrayT& D) const;

    const dArrayT electricFieldPiezoelectrical(const dSymMatrixT& C,
					       const dArrayT& D) const;

    const dMatrixT tangentMechanicalElasticVol(const dSymMatrixT& C) const;
    const dMatrixT tangentMechanicalElasticDev(const dSymMatrixT& C) const;
    const dMatrixT tangentMechanicalElectrical(const dSymMatrixT& C,
					       const dArrayT& D) const;

    const dSymMatrixT rightCauchyGreenDeformation();
    const dArrayT electricDisplacement();
    const dArrayT electricDisplacement(int ip);
    const dArrayT electricDisplacement_last();
    const dArrayT electricDisplacement_last(int ip);

    //
    // data
    //

  public:

    static const char* Name;

  protected:

  private:

    double fShearModulus;
    double fBulkModulus;
    double fElectricPermittivity;
    dMatrixT fPiezoelectricTensor;

    double fEnergyDensity;

    dArrayT fElectricField;
    dArrayT fElectricDisplacement;

    dSymMatrixT fStress;
    dMatrixT fTangentMechanical;
    dMatrixT fTangentPiezoelectrical;
    dMatrixT fTangentElectrical;

  };

} // namespace Tahoe

#include "FSNeoHookePZLinT.i.h"

#endif // _FSNeoHookePZLinT_
