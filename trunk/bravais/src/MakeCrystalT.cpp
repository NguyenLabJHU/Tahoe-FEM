/* $Id: MakeCrystalT.cpp,v 1.2 2002-03-06 01:55:43 jzimmer Exp $ */
#include "MakeCrystalT.h"

#include <iostream>
#include <fstream>
#include "StringT.h"
#include "ifstreamT.h"
#include "ExceptionCodes.h"
#include "VolumeT.h"
#include "BoxT.h"
#include "CrystalLatticeT.h"
#include "FCCT.h"

void MakeCrystalT::Run() {
	cout << "Name of input file?" << "\n";
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

	pv1->SetSize(in);

	StringT latticetype;

	in >> latticetype;

	CrystalLatticeT* pl1;

	if (latticetype=="FCC") {
		pl1 = new FCCT;
	}
	else {
		throw eBadInputValue;
	}	

	pl1->SetBasis();
	pl1->SetLatticeParameters(in);
	pl1->CalculateDensity();

	pv1->DefineBoundary(pl1);
	pv1->CalculateVolume();
        cout << "Volume equals " << pv1->GetVolume() << " cubic angstroms" << endl;

	pv1->FillVolume(pl1);
	pv1->WriteFile();


	in.close();



}

