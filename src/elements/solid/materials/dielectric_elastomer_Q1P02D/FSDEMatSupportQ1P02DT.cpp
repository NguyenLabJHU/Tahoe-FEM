#include "FSDEMatSupportQ1P02DT.h"
#include "FSDielectricElastomerQ1P02DT.h"

namespace Tahoe{

  //
  //
  //
  void
  FSDEMatSupportQ1P02DT::SetContinuumElement(const ContinuumElementT* p)
  {
    //
    // inherited
    //
    FSMatSupportT::SetContinuumElement(p);

    //
    // cast to finite deformation DE pointer
    //
    fFSDielectricElastomerQ1P02D =
      dynamic_cast<const FSDielectricElastomerQ1P02DT*>(p);

  }

} //namespace Tahoe
