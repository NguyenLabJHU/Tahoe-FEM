/* $Id: CrystalLatticeT.cpp,v 1.4 2002-07-25 23:48:07 saubry Exp $ */
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

CrystalLatticeT::CrystalLatticeT(int nlsd, int nuca,dArrayT vec_rot) 
{
	nLSD = nlsd;
	nUCA = nuca;
        vBasis.Dimension(nLSD,nUCA);
        vLatticeParameters.Dimension(nLSD);

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
  double Pi = acos(-1.0);

  int im = A.MajorDim();
  int jm = A.MinorDim();

  cout << "im="<<im<<"\n";
  cout << "jm="<<jm<<"\n";

  if(jm != nLSD) throw eSizeMismatch;
  if(norm_vec <= 1.e-5) return A;

  dArrayT vectorIN(jm);
  dArrayT vectorOUT(jm);
  dArray2DT B(im,jm);

  if(jm==2)
    {
      double angle = atan2(vector_rotation[1],vector_rotation[0]);
      angle = 180.0*angle/Pi;

      Rotate2DT R(angle);
      
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
      double angle1=0.0,angle2=0.0;

      angle1 = acos(vector_rotation[2]);
      if( fabs(sin(angle1)) >= 1.e-5) 
	{
	  double sig = vector_rotation[1]/sin(angle1);
	  angle2 = acos(vector_rotation[0])/sin(angle1);
	  //if(sig <= 0.0) angle2 = 2*Pi - angle1;
	}

      angle1 = 180.0*angle1/Pi;
      angle2 = 180.0*angle2/Pi;

      cout << "angle1="<<angle1<<"\n";
      cout << "angle2="<<angle2<<"\n";

      Rotate3DT R(angle1,angle2);
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

/*
dArray2DT  CrystalLatticeT::BasisRotation(dArray2DT A)
{
  double Pi = acos(-1.0);

  int im = A.MajorDim();
  int jm = A.MinorDim();

  if(im != nLSD) throw eSizeMismatch;
  if(norm_vec <= 1.e-5) return A;

  dArrayT vectorIN(im);
  dArrayT vectorOUT(im);
  dArray2DT B(im,jm);

  if(im==2)
    {
      double angle = atan2(vector_rotation[1],vector_rotation[0]);
      angle = 180.0*angle/Pi;

      Rotate2DT R(angle);
      
      for (int j=0; j<jm; j++)
	{
	  for (int i=0; i<im; i++)  vectorIN[i] = A(i,j);
	  vectorOUT = R.RotateVectorOut(vectorIN);
	  for (int i=0; i<im; i++) 
	    B(i,j) = vectorOUT[i];
	}
    }
  else
    {
      double angle1=0.0,angle2=0.0;

      angle1 = acos(vector_rotation[2]);
      if( fabs(sin(angle1)) >= 1.e-5) 
	{
	  double sig = vector_rotation[1]/sin(angle1);
	  angle2 = acos(vector_rotation[0])/sin(angle1);
      }

      angle1 = 180.0*angle1/Pi;
      angle2 = 180.0*angle2/Pi;

      Rotate3DT R(angle1,angle2);
      for (int j=0; j<jm; j++)
	{
	  for (int i=0; i<im; i++)  vectorIN[i] = A(i,j);
	  vectorOUT = R.RotateVectorOut(vectorIN);
	  for (int i=0; i<im; i++) 
	    B(i,j) = vectorOUT[i];
	}
    }

  return B;
}
*/
