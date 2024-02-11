#include "FSDEMatSupportQ1P0T.h"
#include "FSDielectricElastomerQ1P0T.h"

namespace Tahoe{

  //
  //
  //
  void
  FSDEMatSupportQ1P0T::SetContinuumElement(const ContinuumElementT* p)
  {
    //
    // inherited
    //
    FSMatSupportT::SetContinuumElement(p);

    //
    // cast to finite deformation DE pointer
    //
    fFSDielectricElastomerQ1P0 =
      dynamic_cast<const FSDielectricElastomerQ1P0T*>(p);

  }

} //namespace Tahoe
