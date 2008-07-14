//
// $Id: FSHuWashizuUSCT.i.h,v 1.1 2008-07-14 16:13:25 lxmota Exp $
//
// $Log: not supported by cvs2svn $
//

namespace Tahoe {

  //
  //
  //
  inline
  FSHuWashizuUSCT::FSHuWashizuUSCT(const ElementSupportT& support):
    FiniteStrainT(support),
    fFSMatSupport(0)
  {

    SetName("Finite-deformation Hu-Washizu uSC");
    Initialize();

  }


} // namespace Tahoe
