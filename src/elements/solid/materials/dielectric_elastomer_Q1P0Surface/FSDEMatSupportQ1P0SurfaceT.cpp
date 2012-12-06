#include "FSDEMatSupportQ1P0SurfaceT.h"
#include "FSDielectricElastomerQ1P0SurfaceT.h"

namespace Tahoe{

  //
  //
  //
  void
  FSDEMatSupportQ1P0SurfaceT::SetContinuumElement(const ContinuumElementT* p)
  {
    //
    // inherited
    //
    FSMatSupportT::SetContinuumElement(p);

    //
    // cast to finite deformation DE pointer
    //
    fFSDielectricElastomerQ1P0Surface =
      dynamic_cast<const FSDielectricElastomerQ1P0SurfaceT*>(p);

  }

} //namespace Tahoe
