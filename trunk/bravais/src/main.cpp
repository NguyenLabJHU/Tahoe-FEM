// DEVELOPMENT
/* $Id: main.cpp,v 1.9 2004-08-18 19:53:00 bsun Exp $ */
#include <iostream.h>
#include "ExceptionCodes.h"
#include "PeriodicTableT.h"
#include "StringT.h"

#include "MakeCrystalT.h"

#include "MeshAtomT.h"


    int main()
   {
      MakeCrystalT MC;
      try
      {
         MC.Run();
      }
      
      
          catch (int ErrorCode) 
         {
            cout << "\n\n Exiting due to error . . . ";
            switch (ErrorCode) 
            {
               case eBadInputValue:
                  cout << " Bad Input Value\n";
                  break;
               case eOutOfRange:
                  cout << " Out of Range\n";
                  break;
               case eSizeMismatch:
                  cout << " Size Mismatch\n";
                  break;
               case eOutOfMemory:
                  cout << " Out of Memory\n";
                  break;
               case eDatabaseFail:
                  cout << " Error with database\n";
                  break;
            }
            cout << "\n\n Game Over\n\n";
         }
   
   
      return 0;
   }	
