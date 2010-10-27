#include "FSDEMatSupportT.h"
#include "FSDielectricElastomerT.h"

namespace Tahoe{

  //
  //
  //
  void
  FSDEMatSupportT::SetContinuumElement(const ContinuumElementT* p)
  {

    //
    // inherited
    //
    FSMatSupportT::SetContinuumElement(p);

    //
    // cast to finite deformation DE pointer
    //
    fFSDielectricElastomer =
      dynamic_cast<const FSDielectricElastomerT*>(p);

  }

} //namespace Tahoe
