#include "BCCT.h"
#include "CrystalLatticeT.h"
#include <iostream>
#include <fstream>
#include "StringT.h"
#include "ifstreamT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"

BCCT::BCCT(int nlsd,int nuca,double alat,dArrayT vec_rot) : CrystalLatticeT(nlsd,nuca,vec_rot)
{
  if (nlsd==2)
    {
      if(nuca != 1) {cout << "Wrong nuca\n"; throw eSizeMismatch;}
      vBasis(0,0) = 0.0;
      vBasis(1,0) = 0.0;
    }
  
  if (nlsd==3) 
    {
      if(nuca != 2) {cout << "Wrong nuca\n"; throw eSizeMismatch;}
      vBasis(0,0) = 0.0;
      vBasis(1,0) = 0.0;
      vBasis(2,0) = 0.0;
      vBasis(0,1) = 0.5;
      vBasis(1,1) = 0.5;
      vBasis(2,1) = 0.5;
    }
 
  for (int i=0; i<nlsd; i++) vLatticeParameters[i] = alat;
}


BCCT::BCCT(const BCCT& source) : CrystalLatticeT(source.nLSD,source.nUCA,source.vector_rotation)
{
  for (int i=0; i<source.nLSD; i++) 
    for (int j=0; j<source.nUCA; j++) 
      vBasis(i,j) = source.vBasis(i,j);
 
  for (int i=0; i<source.nLSD; i++) 
    vLatticeParameters[i] = source.vLatticeParameters[i];
}

const dArray2DT& BCCT::GetBasis() {
	if (!&vBasis) {
		cout << "vBasis doesn't exist!" << "\n";
		throw eDatabaseFail;
	}
	else {
		return vBasis;		
	}

}

const dArrayT& BCCT::GetLatticeParameters() {
	if (!&vLatticeParameters) {
		cout << "vLatticeParameters doesn't exist!" << "\n";	
		throw eDatabaseFail;
	}
	else {
		return vLatticeParameters;
	}
}	


