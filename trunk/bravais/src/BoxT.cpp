/* $Id: BoxT.cpp,v 1.1 2002-03-06 01:55:43 jzimmer Exp $ */
#include "BoxT.h"

#include "VolumeT.h"
#include <iostream>
#include <fstream>
#include "ifstreamT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "CrystalLatticeT.h"

BoxT::BoxT(int n) : VolumeT(n) {
}

BoxT::~BoxT() {
}

void BoxT::SetSize(ifstreamT& in) {
	length.Dimension(nSD);
	ncells.Dimension(nSD);
	in >> dimunits;
	if (dimunits) {
		switch(nSD) {
		 case 2:
		  in >> length[0] >> length[1];
		  break;
		 case 3:
		  in >> length[0] >> length[1] >> length[2];
		  break;
		}
	}
	else {
		switch(nSD) {
		 case 2:
		  in >> ncells[0] >> ncells[1];
		  break;
		 case 3:
		  in >> ncells[0] >> ncells[1] >> ncells[2];
		  break;
		}
	}
}

void BoxT::DefineBoundary(CrystalLatticeT* pcl) {
	int nlsd = pcl->GetNLSD();
	if (nlsd!=nSD) throw eSizeMismatch;
	const dArrayT& vLP = pcl->GetLatticeParameters();
	if (dimunits) {
		for(int i=0;i<nSD;i++) {
			double d1 = length[i]/vLP[i];
			int i1 = (int) (2.0*d1+1.0)/2.0;
			ncells[i] = i1;
		}
	}
	else {
		for(int i=0;i<nSD;i++) {
			length[i] = vLP[i]*ncells[i];
		}
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

void BoxT::FillVolume(CrystalLatticeT* pcl) {
	int nuca = pcl->GetNUCA();
	const dArrayT& vLP = pcl->GetLatticeParameters();
	const dArray2DT& vB = pcl->GetBasis();
	int p,q,r,a = 0;
	int count = 1.1*(pcl->GetDensity())*volume;
	double offset = 0.0;
	atom_coord.Dimension(nSD,count);	
	if (nSD==2) {
		for (p=-(ncells[1]/2);p<(ncells[1]/2);p++) {
		for (q=-(ncells[0]/2);q<(ncells[0]/2);q++) {
			for (int m=0;m<nuca;m++) {
				atom_coord[0,a] = ((double) q + vB(0,m))*vLP[0] + offset;
				atom_coord[1,a] = ((double) p + vB(1,m))*vLP[1] + offset;
				a++;
			}
		}
		}
	}
	else if (nSD==3) {
		for (p=-(ncells[2]/2);p<(ncells[2]/2);p++) {
		for (q=-(ncells[1]/2);q<(ncells[1]/2);q++) {
		for (r=-(ncells[0]/2);r<(ncells[0]/2);r++) {
			for (int m=0;m<nuca;m++) {
				atom_coord[0,a] = ((double) r + vB(0,m))*vLP[0] + offset;
				atom_coord[1,a] = ((double) q + vB(1,m))*vLP[1] + offset;
				atom_coord[2,a] = ((double) p + vB(2,m))*vLP[2] + offset;
				a++;
			}
		}
		}
		}
	}

	cout << "The number of atoms anticipated was " << count << "\n";
	cout << "The number of atoms created is " << a << "\n";

	nATOMS = a;

}
		
