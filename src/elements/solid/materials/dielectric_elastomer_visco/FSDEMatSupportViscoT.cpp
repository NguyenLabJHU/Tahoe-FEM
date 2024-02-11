#include "FSDEMatSupportViscoT.h"
#include "FSDielectricElastomerViscoT.h"

namespace Tahoe{

  //
  //
  //
  void
  FSDEMatSupportViscoT::SetContinuumElement(const ContinuumElementT* p)
  {
    //
    // inherited
    //
    FSMatSupportT::SetContinuumElement(p);

    //
    // cast to finite deformation DE pointer
    //
    fFSDielectricElastomerVisco =
      dynamic_cast<const FSDielectricElastomerViscoT*>(p);

  }

} //namespace Tahoe
