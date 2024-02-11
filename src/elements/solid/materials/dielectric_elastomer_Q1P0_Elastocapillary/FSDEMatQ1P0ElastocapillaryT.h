
#if !defined(_FSDEMatQ1P0ElastocapillaryT_)
#define _FSDEMatQ1P0ElastocapillaryT_

#include <cassert>

#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "FSSolidMatT.h"

namespace Tahoe {

  class FSDEMatQ1P0ElastocapillaryT: public FSSolidMatT
  {

  public:

    // constructors
    FSDEMatQ1P0ElastocapillaryT();

    // set parameters
    void DefineParameters(ParameterListT& list) const;
    void TakeParameterList(const ParameterListT& list);

    // information about subordinate parameter lists
    virtual void DefineSubs(SubListT& sub_list) const;

    // Interface required by Tahoe
//    double StrainEnergyDensity();

    // spatial mechanical tangent modulus
    virtual const dMatrixT& c_ijkl();

    // Cauchy stress
    virtual const dSymMatrixT& s_ij();

  protected:

  private:

    void Initialize();

    const dMatrixT RightCauchyGreenDeformation();

  public:

    static const char* Name;

  protected:

  private:
	
	double fSurfTension;
    dSymMatrixT fCapillaryStress;
    dMatrixT fCapillaryTangent;

  };

} // namespace Tahoe

#endif // _FSDEMatQ1P0ElastocapillaryT_
