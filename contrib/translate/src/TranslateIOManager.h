/* $Id: TranslateIOManager.h,v 1.1 2001-09-01 00:05:36 paklein Exp $ */

#ifndef _TRANSLATE_IOMANAGER_H_
#define _TRANSLATE_IOMANAGER_H_

#include "IOManager.h"

class TranslateIOManager : public IOManager
{
 public:

  TranslateIOManager (ostream& outfile);

  virtual void Interactive (void);
};

#endif
