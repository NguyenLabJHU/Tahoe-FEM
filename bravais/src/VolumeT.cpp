/* $Id: VolumeT.cpp,v 1.2 2002-03-06 01:55:43 jzimmer Exp $ */
#include "VolumeT.h"
#include <iostream>
#include <fstream>
#include "ifstreamT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "ArrayT.h"
#include "EnSightOutputT.h"
#include "OutputSetT.h"
#include "StringT.h"

VolumeT::VolumeT(int n) {
	nSD = n;
	nATOMS = 0;
}

VolumeT::~VolumeT() {
}

int VolumeT::GetDimensions() {
	return nSD;
}

double VolumeT::GetVolume() {
	return volume;
}

void VolumeT::WriteFile() {

	ArrayT<iArray2DT> fAtomSets(1);
	fAtomSets[0].Dimension(nATOMS,1);
	fAtomSets[0].SetValueToPosition();

	ArrayT<iArray2DT*> connects(1);
	connects[0] = fAtomSets.Pointer(0);
	ArrayT<StringT> ID(1);
	ID[0] = "nickel";

	


}

