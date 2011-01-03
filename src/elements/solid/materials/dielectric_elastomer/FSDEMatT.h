
#if !defined(_FSDEMatT_)
#define _FSDEMatT_

#include <cassert>

#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "NL_E_MatT.h"
#include "FSDEMatSupportT.h"

namespace Tahoe {

  class FSDEMatT: public NL_E_MatT
  {

    //
    // methods
    //

  public:

    //
    // constructors
    //
    FSDEMatT();

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

    //
    // material electromechanical tangent modulus
    //
    virtual const dMatrixT& E_IJK();

    //
    // material electric tangent modulus
    //
    virtual const dMatrixT& B_IJ();

    //
    // Second Piola-Kirchhoff stress
    //
    virtual const dSymMatrixT& S_IJ();

    //
    // electric displacement
    //
    virtual const dArrayT& D_I();

    //
    // electric field
    //
    virtual const dArrayT& E_I();

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
    // accessors and mutators for material constants
    //
    void SetElectricPermittivity(double epsilon);
    double GetElectricPermittivity() const;

    void SetFSDEMatSupport(const FSDEMatSupportT* support);

  protected:

	/* compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* compute the symetric 2nd Piola-Kirchhoff reduced index vector */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);

	/* returns the strain energy density for the specified strain */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);

    const FSDEMatSupportT* fFSDEMatSupport;

  private:

    void Initialize();

    const dMatrixT RightCauchyGreenDeformation();
    const dArrayT ElectricField();
    const dArrayT ElectricField(int ip);

    //
    // data
    //

  public:

    static const char* Name;

  protected:

  private:

    double fElectricPermittivity;
    double fEnergyDensity;
    double fMu;
    double fNrig;
    double fLambda;
	double fKappa;

    dArrayT fElectricField;
    dArrayT fElectricDisplacement;
	dArrayT fParams;
	
    dSymMatrixT fStress;
    dMatrixT fTangentMechanicalElec;
    dMatrixT fTangentMechanical;
    dMatrixT fTangentElectromechanical;
    dMatrixT fTangentElectrical;

  };

} // namespace Tahoe

#include "FSDEMatT.i.h"

#endif // _FSDEMatT_
