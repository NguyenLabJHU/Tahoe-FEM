/* $Id: VTKFrameT.cpp,v 1.1 2001-10-23 00:23:08 recampb Exp $ */

#include "VTKFrameT.h"

/* constructor */
VTKFrameT::VTKFrameT(void): fAge(0)
{
  /* set up console modifiable variables */
  iAddVariable("age", fAge);
  iAddVariable("str_len", (const int) fLength);
  iAddVariable("string", (const char*) fArray);
}
