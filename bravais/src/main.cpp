/* $Id: main.cpp,v 1.4 2002-07-24 01:14:59 saubry Exp $ */
#include <iostream>
#include "ExceptionCodes.h"
#include "PeriodicTableT.h"
#include "StringT.h"

#include "MakeCrystalT.h"

#include "MeshAtom.h"


int main()
{
  MakeCrystalT MC;
  MC.Run();
  
  return 0;
}	
