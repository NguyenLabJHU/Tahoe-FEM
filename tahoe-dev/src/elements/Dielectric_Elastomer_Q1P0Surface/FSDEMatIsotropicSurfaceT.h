
#if !defined(_FSDEMatIsotropicSurfaceT_)
#define _FSDEMatIsotropicSurfaceT_

#include <cassert>

#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "FSSolidMatT.h"

namespace Tahoe {

  class FSDEMatIsotropicSurfaceT: public FSSolidMatT
  {

  public:

    // constructors
    FSDEMatIsotropicSurfaceT();

    // set parameters
    void DefineParameters(ParameterListT& list) const;
    void TakeParameterList(const ParameterListT& list);

    // information about subordinate parameter lists
    virtual void DefineSubs(SubListT& sub_list) const;

    // Interface required by Tahoe
    double StrainEnergyDensity();

    // material mechanical tangent modulus
    virtual const dMatrixT& C_IJKL();

    // Second Piola-Kirchhoff stress
    virtual const dSymMatrixT& S_IJ();
    
    // spatial mechanical tangent modulus
    virtual const dMatrixT& c_ijkl();
    
    // Cauchy stress
    virtual const dSymMatrixT& s_ij();

    // pressure associated with the last computed stress
    double Pressure() const;

  protected:

	/* FDKStV.h */
	/* set (material) tangent modulus */
	void SetModulus();

  private:

    void Initialize();

	/** return true if material implementation supports imposed thermal
	 * strains. This material does support multiplicative thermal
	 * strains. */
	virtual bool SupportsThermalStrain(void) const { return true; };

  public:

    static const char* Name;

  protected:

  private:

	double fE;
	double fNu;
	
	/* From FSDEMatQ1P0SurfaceT.h */
    double fElectricPermittivity;
    double fMu;
    double fNrig;
    double fLambda;
	double fKappa;	
	
    dSymMatrixT fIsotropicStress, fStress1;
	dMatrixT fModulusKStV;
	FrameT	fLastCall;
	dMatrixT fTM;

  };

} // namespace Tahoe

#endif // _FSDEMatIsotropicSurfaceT_
