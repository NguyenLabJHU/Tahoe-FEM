/* $Id: TranslateIOManager.cpp,v 1.1 2001-09-01 00:05:36 paklein Exp $  */

#include "TranslateIOManager.h"

TranslateIOManager::TranslateIOManager (ostream& out) :
  IOManager (out)
{
}

void TranslateIOManager::Interactive (void)
{
  InteractiveIO ();
}

