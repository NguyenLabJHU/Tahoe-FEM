/* $Id: VTKConsoleT.h,v 1.1 2001-09-19 01:12:35 recampb Exp $ */

#ifndef _VTK_CONSOLE_T_H_
#define _VTK_CONSOLE_T_H_

/* base class */
#include "iConsoleObjectT.h"

class VTKConsoleT: public iConsoleObjectT
{
 public:

  /* constructor */
  VTKConsoleT(void);

  /* execute given command - returns false on fail */
  virtual bool iDoCommand(const StringT& command, StringT& line);

 private:
  
  int fInteger;
  double fDouble;
  StringT fString;
};

#endif
