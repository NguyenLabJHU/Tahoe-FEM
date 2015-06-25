
#if !defined(_FSDEMatQ1P02DT_)
#define _FSDEMatQ1P02DT_

#include <cassert>

#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "FSSolidMatT.h"
#include "FSDEMatSupportQ1P02DT.h"

namespace Tahoe {

  class FSDEMatQ1P02DT: public FSSolidMatT
  {

  public:

    // constructors
    FSDEMatQ1P02DT();

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

    void SetFSDEMatSupportQ1P02D(const FSDEMatSupportQ1P02DT* support);

  protected:

    const FSDEMatSupportQ1P02DT* fFSDEMatSupportQ1P02D;

  private:

    void Initialize();

    /* ------------- */
    const dMatrixT DeformationMatrix();
    /* ------------- */

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
	double fGamma;
	double fT_0;

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

#include "FSDEMatQ1P02DT.i.h"

#endif // _FSDEMatQ1P02DT_
