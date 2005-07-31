/* $Id: MakeCrystalT.cpp,v 1.1.1.1 2002-02-28 02:13:08 jzimmer Exp $ */
#include "MakeCrystalT.h"

#include <iostream>
#include <fstream>
#include "StringT.h"
#include "ifstreamT.h"
#include "ExceptionCodes.h"
#include "VolumeT.h"

void MakeCrystalT::Run() {
	cout << "Name of input file?" < "\n";
	StringT inputfile;
	cin >> inputfile;

	ifstreamT in('%');
	in.open(inputfile);

	int nsd; 
	StringT shape;

	in >> nsd >> shape;

	VolumeT* pv1;

	switch(nsd) {
	 case 2:
	  if (shape=="box") {
		pv1 = new BoxT(nsd); 
	  }
	  else { 
		throw eBadInputValue;
	  }
         break;
         case 3:
	  if (shape=="box") {
		pv1 = new BoxT(nsd); 
	  }
	  else { 
		throw eBadInputValue;
	  }
         break;
	 default :
	  throw eBadInputValue;
	}

	pv1->DefineBoundary(in);
	pv1->CalculateVolume();
        cout << "Volume equals " << pv1->GetVolume() << " cubic angstroms" << endl;
	pv1->FillVolume();
	pv1->WriteFile();


	in.close();



}

