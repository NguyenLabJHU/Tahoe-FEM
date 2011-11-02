
#if !defined(_FSDEMatViscoT_)
#define _FSDEMatViscoT_

#include <cassert>

#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "FSSolidMatT.h"
#include "FSDEMatSupportViscoT.h"

namespace Tahoe {

  class FSDEMatViscoT: public FSSolidMatT
  {

    //
    // methods
    //

  public:

    //
    // constructors
    //
    FSDEMatViscoT();

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

    void SetFSDEMatSupportVisco(const FSDEMatSupportViscoT* support);

  protected:

    const FSDEMatSupportViscoT* fFSDEMatSupportVisco;

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

#include "FSDEMatViscoT.i.h"

#endif // _FSDEMatViscoT_
