// DEVELOPMENT
/* $Id: MakeCrystalT.h,v 1.4 2005-02-04 00:59:27 rjones Exp $ */

#ifndef _MAKE_CRYSTAL_T_H_
#define _MAKE_CRYSTAL_T_H_

#include "StringT.h"

using namespace Tahoe;

    class MakeCrystalT {
   public:
       MakeCrystalT() { }
       ~MakeCrystalT() { }
   
      void Run(StringT& input_file);
   };

#endif
