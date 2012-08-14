
#if !defined(_FSDEMatQ1P0T_)
#define _FSDEMatQ1P0T_

#include <cassert>

#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "FSSolidMatT.h"
#include "FSDEMatSupportQ1P0T.h"

namespace Tahoe {

  class FSDEMatQ1P0T: public FSSolidMatT
  {

  public:

    // constructors
    FSDEMatQ1P0T();

    // set parameters
    void DefineParameters(ParameterListT& list) const;
    void TakeParameterList(const ParameterListT& list);

    // information about subordinate parameter lists
    virtual void DefineSubs(SubListT& sub_list) const;

    // Interface required by Tahoe
    double StrainEnergyDensity();

    // material mechanical tangent modulus
    virtual const dMatrixT& C_IJKL();

    // material electromechanical tangent modulus
    virtual const dMatrixT& E_IJK();

    // material electric tangent modulus
    virtual const dMatrixT& B_IJ();

    // Second Piola-Kirchhoff stress
    virtual const dSymMatrixT& S_IJ();

    // electric displacement
    virtual const dArrayT& D_I();

    // electric field
    virtual const dArrayT& E_I();

    // spatial mechanical tangent modulus
    virtual const dMatrixT& c_ijkl();

    // Cauchy stress
    virtual const dSymMatrixT& s_ij();

	// Q1P0 STUFF
	virtual const dMatrixT& b_ij();
	virtual const dArrayT& d_i();
	virtual const dMatrixT& e_ijk();

    // pressure associated with the last computed stress
    double Pressure() const;
	
    // accessors and mutators for material constants
    void SetElectricPermittivity(double epsilon);
    double GetElectricPermittivity() const;

    void SetFSDEMatSupportQ1P0(const FSDEMatSupportQ1P0T* support);

  protected:

    const FSDEMatSupportQ1P0T* fFSDEMatSupportQ1P0;

  private:

    void Initialize();

    const dMatrixT RightCauchyGreenDeformation();
    const dArrayT ElectricField();
    const dArrayT ElectricField(int ip);

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
    dMatrixT fTangentElectromechanicalSpatial;    
    dMatrixT fTangentElectrical;

  };

} // namespace Tahoe

#include "FSDEMatQ1P0T.i.h"

#endif // _FSDEMatQ1P0T_
