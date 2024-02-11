#include "FSDEMatSupportQ1P0ViscoT.h"
#include "FSDielectricElastomerQ1P0ViscoT.h"

namespace Tahoe{

  //
  //
  //
  void
  FSDEMatSupportQ1P0ViscoT::SetContinuumElement(const ContinuumElementT* p)
  {
    //
    // inherited
    //
    FSMatSupportT::SetContinuumElement(p);

    //
    // cast to finite deformation DE pointer
    //
    fFSDielectricElastomerQ1P0Visco =
      dynamic_cast<const FSDielectricElastomerQ1P0ViscoT*>(p);

  }

} //namespace Tahoe
