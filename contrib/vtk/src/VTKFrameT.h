/* $Id: VTKFrameT.h,v 1.1 2001-10-23 00:23:08 recampb Exp $ */

#ifndef _VTK_FRAME_T_H_
#define _VTK_FRAME_T_H_

/* base class */
#include "iConsoleObjectT.h"
#include "StringT.h"

class VTKFrameT: public StringT, public iConsoleObjectT
{
 public:

  /* constructor */
  VTKFrameT(void);

 private:

  /* dummy variable */
  int fAge;
};

#endif
