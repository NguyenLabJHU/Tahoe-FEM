/* $Id: CrystalLatticeT.cpp,v 1.5 2002-07-29 19:16:26 saubry Exp $ */
#include "CrystalLatticeT.h"

#include <iostream>
#include <fstream>
#include "StringT.h"
#include "ifstreamT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"

#include "dMatrixT.h"
#include "Rotate2DT.h"
#include "Rotate3DT.h"

CrystalLatticeT::CrystalLatticeT(int nlsd, int nuca,
				 dArrayT vec_rot,double angle) 
{
	nLSD = nlsd;
	nUCA = nuca;
        vBasis.Dimension(nLSD,nUCA);
        vLatticeParameters.Dimension(nLSD);

	if(nLSD == 3)
	  {
	    // Direction in 3D only
	    if(vec_rot.Length() != nLSD) throw eSizeMismatch;
	    vector_rotation.Dimension(nLSD);

	    norm_vec = 0.0;
	    for (int i=0; i<nLSD; i++)
	      norm_vec += vec_rot[i]*vec_rot[i];
	    norm_vec = sqrt(norm_vec);
	    
	    if(norm_vec != 0) 
	      {
		for (int i=0; i<nLSD; i++) 
      		  vector_rotation[i] = vec_rot[i]/norm_vec;
	      }
	  }

	// Assumed here that the angle is in degree
	angle_rotation = angle;
}

CrystalLatticeT::CrystalLatticeT(const CrystalLatticeT& source) 
{

  nLSD = source.nLSD;
  nUCA = source.nUCA;

  vBasis.Dimension(nLSD,nUCA);
  vBasis = source.vBasis;

  vLatticeParameters.Dimension(nLSD);
  vLatticeParameters = source.vLatticeParameters;

  density = source.density;
}


void CrystalLatticeT::CalculateDensity() 
{
	double ucvolume = 1.0;
	for (int i=0; i<nLSD; i++) 
	  ucvolume *= vLatticeParameters[i];
	density = ((double) nUCA)/ucvolume;
}

double CrystalLatticeT::GetDensity() 
{
	return density;
}

dArray2DT  CrystalLatticeT::BasisRotation(dArray2DT A)
{
  int im = A.MajorDim();
  int jm = A.MinorDim();

  if(jm != nLSD) throw eSizeMismatch;
  if(fabs(angle_rotation) <= 1.e-5 && norm_vec <= 1.e-5) return A;

  dArrayT vectorIN(jm);
  dArrayT vectorOUT(jm);
  dArray2DT B(im,jm);

  if(jm==2)
    {
      Rotate2DT R(angle_rotation);
      
      for (int i=0; i<im; i++)
	{
	  for (int j=0; j<jm; j++)  vectorIN[j] = A(i,j);
	  vectorOUT = R.RotateVectorOut(vectorIN);
	  for (int j=0; j<jm; j++) 
	    B(i,j) = vectorOUT[j];
	}
    }
  else if(jm==3)
    {
      Rotate3DT R(vector_rotation,angle_rotation);
      for (int i=0; i<im; i++)
	{
	  for (int j=0; j<jm; j++)  vectorIN[j] = A(i,j);
	  vectorOUT = R.RotateVectorOut(vectorIN);
	  for (int j=0; j<jm; j++) 
	    B(i,j) = vectorOUT[j];
	}
    }

  return B;
}
