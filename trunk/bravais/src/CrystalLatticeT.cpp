/* $Id: CrystalLatticeT.cpp,v 1.6 2002-08-02 02:07:49 saubry Exp $ */
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
				 dArray2DT mat_rot,double angle) 
{
	nLSD = nlsd;
	nUCA = nuca;
        vBasis.Dimension(nLSD,nUCA);
        vLatticeParameters.Dimension(nLSD);

	if(nLSD == 3)
	  {
	    // Directions in 3D only
	    if(mat_rot.MinorDim() != nLSD) throw eSizeMismatch;
	    if(mat_rot.MajorDim() != nLSD) throw eSizeMismatch;
	    matrix_rotation.Dimension(nLSD,nLSD);
	    
	    norm_vec.Dimension(nLSD);

	    for (int j=0; j<nLSD; j++)
	      {
		norm_vec[j] = 0.0;
		for (int i=0; i<nLSD; i++)
		  norm_vec[j] += mat_rot(i,j)*mat_rot(i,j);
		norm_vec[j] = sqrt(norm_vec[j]);
	      }



	    for (int j=0; j<nLSD; j++)
	      {	
		if(norm_vec[j] != 0) 
		  {
		    for (int i=0; i<nLSD; i++) 
		      for (int j=0; j<nLSD; j++) 
			matrix_rotation(i,j) = mat_rot(i,j)/norm_vec[j];
		  }
	      }

	    // check input matrix
	    double yz =matrix_rotation(1,0)*matrix_rotation(2,1)-matrix_rotation(1,1)*matrix_rotation(2,0);
	    double xz =matrix_rotation(0,1)*matrix_rotation(2,0)-matrix_rotation(0,0)*matrix_rotation(2,1); 
	    double xy =matrix_rotation(0,0)*matrix_rotation(1,1)-matrix_rotation(0,1)*matrix_rotation(1,0);  
	    if ( fabs(yz-matrix_rotation(0,2)) <= 1.e-5 ||
		 fabs(xz-matrix_rotation(1,2)) <= 1.e-5 ||
		 fabs(xy-matrix_rotation(2,2)) <= 1.e-5 ) 

	      {
		/*
		  cout << "Matrix of rotation is not right.\n Changing it ...\n";
		  cout << "was " << matrix_rotation(0,2) << "  " 
	 	       << matrix_rotation(1,2) << "  " 
		       << matrix_rotation(2,2) << "\n";
		  cout << "is now " << yz << "  " 
		       << xz << "  " 
		       << xy << "\n";
		*/
		matrix_rotation(0,2) = yz;
		matrix_rotation(1,2) = xz;
		matrix_rotation(2,2) = xy;
	      }


	    angle_rotation = 0;
	  }
	else if(nLSD == 2)  
	  {
	    // Assumed here that the angle is in degree
	    angle_rotation = angle;
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

/*
dArray2DT  CrystalLatticeT::BasisRotation(dArray2DT A)
{
  int im = A.MajorDim();
  int jm = A.MinorDim();

  cout << "im=" << im << "\n";
  cout << "jm=" << jm << "\n";

  if(im != nLSD) throw eSizeMismatch;

  dArrayT vectorIN(jm);
  dArrayT vectorOUT(jm);
  dArray2DT B(im,jm);

  if(im==2)
    {
      if(fabs(angle_rotation) <= 1.e-5) return A;
      Rotate2DT R(angle_rotation);
      
      for (int j=0; j<jm; j++)
	{
	  for (int i=0; i<im; i++)  vectorIN[i] = A(i,j);
	  vectorOUT = R.RotateVectorIn(vectorIN);
	  for (int i=0; i<im; i++) 
	    B(i,j) = vectorOUT[i];
	}
    }
  else if(im==3)
    {
      double norm = norm_vec[0] + norm_vec[1] + norm_vec[2];
      if (norm < 1.e-5) return A;

      Rotate3DT R;
      R.GiveTransfoMatrix(matrix_rotation);
      for (int j=0; j<jm; j++)
	{
	  for (int i=0; i<im; i++)  vectorIN[i] = A(i,j);
	  vectorOUT = R.RotateVectorIn(vectorIN);
	  for (int i=0; i<im; i++) 
	    B(i,j) = vectorOUT[i];
	}
    }

  return B;
}
*/

dArray2DT  CrystalLatticeT::BasisRotation(dArray2DT A)
{
  int im = A.MajorDim();
  int jm = A.MinorDim();

  dMatrixT Q(jm,jm);

  if(jm != nLSD) throw eSizeMismatch;

  dArrayT vectorIN(jm);
  dArrayT vectorOUT(jm);
  dArray2DT B(im,jm);

  if(jm==2)
    {
      if(fabs(angle_rotation) <= 1.e-5) return A;
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
      double norm = norm_vec[0] + norm_vec[1] + norm_vec[2];
      if (norm < 1.e-5) return A;
      Rotate3DT R;
      R.GiveTransfoMatrix(matrix_rotation);

      /*
      Q = R.Q();
      cout << Q(0,0) << "  " << Q(1,0) << "  " << Q(2,0) << "\n";
      cout << Q(0,1) << "  " << Q(1,1) << "  " << Q(2,1) << "\n";
      cout << Q(0,2) << "  " << Q(1,2) << "  " << Q(2,2) << "\n\n";

      vectorIN[0]=vectorIN[1]=0;vectorIN[2] = 1;
      Rotate3DT R2(vectorIN,45);
      Q = R2.Q();
      cout << Q(0,0) << "  " << Q(1,0) << "  " << Q(2,0) << "\n";
      cout << Q(0,1) << "  " << Q(1,1) << "  " << Q(2,1) << "\n";
      cout << Q(0,2) << "  " << Q(1,2) << "  " << Q(2,2) << "\n";
      */
      


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
dArray2DT CrystalLatticeT::AxisRotation(dArrayT len)
{
  int jm = len.Length();
  dArrayT vectorIn(jm);
  dArrayT vectorOut(jm);

  dArrayT rotatelen(jm); 

  if(jm==2)
    {
      if(fabs(angle_rotation) <= 1.e-5) return len;
      Rotate2DT R(angle_rotation);    

      vectorIn[0] = len[0];
      vectorIn[1] = 0.0;
      vectorOut = R.RotateVectorOut(vectorIn);
      rotatelen[0] = fabs(vectorOut[1]);

      vectorIn[0] = 0.0;
      vectorIn[1] = len[1];
      vectorOut = R.RotateVectorOut(vectorIn);
      rotatelen[1] = fabs(vectorOut[1]);

      return rotatelen;
    }
  else if(jm==3)
    {
      double norm = norm_vec[0] + norm_vec[1] + norm_vec[2];
      if (norm < 1.e-5) return len;

      Rotate3DT R;
      R.GiveTransfoMatrix(matrix_rotation);

      vectorIn[0] = len[0];
      vectorIn[1] = vectorIn[2] = 0.0;
      vectorOut = R.RotateVectorOut(vectorIn);
      rotatelen[0] = vectorOut[0];

      vectorIn[1] = len[1];
      vectorIn[0] = vectorIn[2] = 0.0;
      vectorOut = R.RotateVectorOut(vectorIn);
      rotatelen[1] = vectorOut[1];

      vectorIn[2] = len[2];
      vectorIn[0] = vectorIn[1] = 0.0;
      vectorOut = R.RotateVectorOut(vectorIn);
      rotatelen[2] = vectorOut[2];

      return rotatelen;
    }
}
*/

