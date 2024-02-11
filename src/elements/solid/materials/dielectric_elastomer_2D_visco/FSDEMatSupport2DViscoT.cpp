#include "FSDEMatSupport2DViscoT.h"
#include "FSDielectricElastomer2DViscoT.h"

namespace Tahoe{

  //
  //
  //
  void
  FSDEMatSupport2DViscoT::SetContinuumElement(const ContinuumElementT* p)
  {
    //
    // inherited
    //
    FSMatSupportT::SetContinuumElement(p);

    //
    // cast to finite deformation DE pointer
    //
    fFSDielectricElastomer2DVisco =
      dynamic_cast<const FSDielectricElastomer2DViscoT*>(p);

  }

} //namespace Tahoe
