/* $Id: FCCT.cpp,v 1.5 2002-08-02 02:07:49 saubry Exp $ */
#include "FCCT.h"
#include "CrystalLatticeT.h"

#include <iostream>
#include <fstream>

#include "ifstreamT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"

FCCT::FCCT(int nlsd,int nuca,double alat,
	   dArray2DT mat_rot,double angle) : 
  CrystalLatticeT(nlsd,nuca,mat_rot,angle)
{
  if (nlsd==2)
    {
      if(nuca != 2) {cout << "Wrong nuca\n"; throw eSizeMismatch;}
      vBasis(0,0) = 0.0;
      vBasis(1,0) = 0.0;

      vBasis(0,1) = 0.5;
      vBasis(1,1) = 0.5;

      //if (fabs(angle_rotation) >=1.e-5) vBasis = BasisRotation(vBasis);

    }
  
  if (nlsd==3) 
    {
      if(nuca != 4) {cout << "Wrong nuca\n"; throw eSizeMismatch;}
      vBasis(0,0) = 0.0;
      vBasis(1,0) = 0.0;
      vBasis(2,0) = 0.0;

      vBasis(0,1) = 0.5;
      vBasis(1,1) = 0.5;
      vBasis(2,1) = 0.0;

      vBasis(0,2) = 0.5;
      vBasis(1,2) = 0.0;
      vBasis(2,2) = 0.5;

      vBasis(0,3) = 0.0;
      vBasis(1,3) = 0.5;
      vBasis(2,3) = 0.5;

      //double norm = norm_vec[0] + norm_vec[1] + norm_vec[2];
      //if (norm > 1.e-5) vBasis = BasisRotation(vBasis);
    }
  
  for (int i=0; i<nlsd; i++) 
    vLatticeParameters[i] = alat;  
}


FCCT::FCCT(const FCCT& source) : CrystalLatticeT(source.nLSD,source.nUCA,
						 source.matrix_rotation,
						 source.angle_rotation)
{
  
  for (int i=0; i<source.nLSD; i++) 
    for (int j=0; j<source.nUCA; j++) 
      vBasis(i,j) = source.vBasis(i,j);
 
  for (int i=0; i<source.nLSD; i++) 
    vLatticeParameters[i] = source.vLatticeParameters[i];
}

const dArray2DT& FCCT::GetBasis() {
	if (!&vBasis) {
		cout << "vBasis doesn't exist!" << "\n";
		throw eDatabaseFail;
	}
	else {
		return vBasis;		
	}

}

const dArrayT& FCCT::GetLatticeParameters() {
	if (!&vLatticeParameters) {
		cout << "vLatticeParameters doesn't exist!" << "\n";	
		throw eDatabaseFail;
	}
	else {
		return vLatticeParameters;
	}
}	










