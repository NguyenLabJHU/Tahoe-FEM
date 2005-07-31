/* $Id: VolumeT.cpp,v 1.1.1.1 2002-02-28 02:13:08 jzimmer Exp $ */
#include "VolumeT.h"
#include <iostream>
#include <fstream>
#include "ifstreamT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

VolumeT::VolumeT(int n) {
	nSD = n;
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
}

BoxT::BoxT(int n) : VolumeT(n) {
}

BoxT::~BoxT() {
}

void BoxT::DefineBoundary(ifstreamT& in) {
	length.Dimension(nSD);
	switch(nSD) {
	 case 2:
	  in >> length[0] >> length[1];
	  break;
	 case 3:
	  in >> length[0] >> length[1] >> length[2];
	  break;
	}

	surfaces.Dimension(nSD,2);
	for (int i=0; i<nSD; i++ ) {
		surfaces(i,0) = -length[i]/2;
		surfaces(i,1) =  length[i]/2;
	}

}

void BoxT::CalculateVolume() {
        switch(nSD) {
	 case 2:
	  volume = length[0]*length[1];
	  break;
	 case 3:
	  volume = length[0]*length[1]*length[2];
	  break;
	}
}

void BoxT::FillVolume() {
}
	
