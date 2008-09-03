//
// $Id: FSHuWashizuUSCT.i.h,v 1.1 2008-09-03 18:40:50 beichuan Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.1  2008/07/14 16:13:25  lxmota
// Initial sources (disabled for now)
//
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
