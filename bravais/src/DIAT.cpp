#include "DIAT.h"
#include "CrystalLatticeT.h"
#include <iostream>
#include <fstream>
#include "StringT.h"
#include "ifstreamT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"

DIAT::DIAT(int nlsd,int nuca,double alat) : CrystalLatticeT(nlsd,nuca)
{
  if (nlsd==2)
    {
      if(nuca != 2) {cout << "Wrong nuca\n"; throw eSizeMismatch;}
      vBasis(0,0) = 0.0;
      vBasis(1,0) = 0.0;
      vBasis(0,1) = 0.25;
      vBasis(1,1) = 0.25;
    }
  
  if (nlsd==3) 
    {
      if(nuca != 8) {cout << "Wrong nuca\n"; throw eSizeMismatch;}
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

      vBasis(0,4) = 0.25;
      vBasis(1,4) = 0.25;
      vBasis(2,4) = 0.25;

      vBasis(0,5) = 0.75;
      vBasis(1,5) = 0.75;
      vBasis(2,5) = 0.25;

      vBasis(0,6) = 0.75;
      vBasis(1,6) = 0.25;
      vBasis(2,6) = 0.75;

      vBasis(0,7) = 0.25;
      vBasis(1,7) = 0.75;
      vBasis(2,7) = 0.75;
    }
 
  for (int i=0; i<nlsd; i++) vLatticeParameters[i] = alat;
}


DIAT::DIAT(const DIAT& source) : CrystalLatticeT(source.nLSD,source.nUCA)
{
  if (source.nLSD==2)
    {
      if(source.nUCA != 2) {cout << "Wrong nuca\n"; throw eSizeMismatch;}
      vBasis(0,0) = 0.0;
      vBasis(1,0) = 0.0;
      vBasis(0,1) = 0.25;
      vBasis(1,1) = 0.25;
    }
  
  if (source.nLSD==3) 
    {
      if(source.nUCA != 8) {cout << "Wrong nuca\n"; throw eSizeMismatch;}
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

      vBasis(0,4) = 0.25;
      vBasis(1,4) = 0.25;
      vBasis(2,4) = 0.25;

      vBasis(0,5) = 0.75;
      vBasis(1,5) = 0.75;
      vBasis(2,5) = 0.25;

      vBasis(0,6) = 0.75;
      vBasis(1,6) = 0.25;
      vBasis(2,6) = 0.75;

      vBasis(0,7) = 0.25;
      vBasis(1,7) = 0.75;
      vBasis(2,7) = 0.75;
  }
 
  for (int i=0; i<source.nLSD; i++) 
    vLatticeParameters[i] = source.vLatticeParameters[i];
}

const dArray2DT& DIAT::GetBasis() {
	if (!&vBasis) {
		cout << "vBasis doesn't exist!" << "\n";
		throw eDatabaseFail;
	}
	else {
		return vBasis;		
	}

}

const dArrayT& DIAT::GetLatticeParameters() {
	if (!&vLatticeParameters) {
		cout << "vLatticeParameters doesn't exist!" << "\n";	
		throw eDatabaseFail;
	}
	else {
		return vLatticeParameters;
	}
}	


