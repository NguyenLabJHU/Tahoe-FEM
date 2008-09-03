//
// $Id: FSPZMatSupportT.cpp,v 1.1 2008-09-03 18:40:50 beichuan Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.1  2008/06/16 18:21:41  lxmota
// Piezoelectric material support. Initial sources.
//
//

#include "FSPZMatSupportT.h"
#include "FSPiezoElectricSolidT.h"

namespace Tahoe{

  //
  //
  //
  void
  FSPZMatSupportT::SetContinuumElement(const ContinuumElementT* p)
  {

    //
    // inherited
    //
    FSMatSupportT::SetContinuumElement(p);

    //
    // cast to finite deformation piezoelectric pointer
    //
    fFSPiezoElectricSolid =
      dynamic_cast<const FSPiezoElectricSolidT*>(p);

  }

} //namespace Tahoe
