#include "FSDEMatSupport2DT.h"
#include "FSDielectricElastomer2DT.h"

namespace Tahoe{

  //
  //
  //
  void
  FSDEMatSupport2DT::SetContinuumElement(const ContinuumElementT* p)
  {
    //
    // inherited
    //
    FSMatSupport2DT::SetContinuumElement(p);

    //
    // cast to finite deformation DE pointer
    //
    fFSDielectricElastomer2D =
      dynamic_cast<const FSDielectricElastomer2DT*>(p);

  }

} //namespace Tahoe
