// DEVELOPMENT
/* $Id: main.cpp,v 1.10 2005-02-04 00:59:27 rjones Exp $ */
#include <iostream.h>
#include "ExceptionCodes.h"
#include "PeriodicTableT.h"
#include "StringT.h"

#include "MakeCrystalT.h"

#include "MeshAtomT.h"


    int main(int argc, char* argv[])
   {
	StringT input_file;
	if (argc != 2)  {
//		cerr << "usage: " << argv[0] << " <infile> \n";
//		abort();
		cout << "Enter the name of the input date file:" << endl;
		cin >> input_file;
	} else {
		input_file = argv[1];
	}

	cout << " using input file : " << input_file << "\n"; 

      MakeCrystalT MC;
      try
      {
         MC.Run(input_file);
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
