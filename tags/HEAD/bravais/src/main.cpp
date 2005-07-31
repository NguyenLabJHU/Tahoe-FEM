/* $Id: main.cpp,v 1.1.1.1 2002-02-28 02:13:08 jzimmer Exp $ */
#include <iostream>
#include "ExceptionCodes.h"
#include "MakeCrystalT.h"

int main()
{
	cout << "Starting up Bravais program..." << endl;

	MakeCrystalT MC;
	try {
		MC.Run();
	}
	
	catch (int ErrorCode) {
		cout << "\n\n Exiting due to error . . . ";
		switch (ErrorCode) {
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

	cout << "Exiting Bravais program." << endl;

	return 0;
}	
