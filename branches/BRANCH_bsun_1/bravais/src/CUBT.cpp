// DEVELOPMENT
#include "CUBT.h"
#include "CrystalLatticeT.h"
#include "StringT.h"
#include "ifstreamT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"

// Lattice parameter: a = b <> c 


CUBT::CUBT(int nlsd,int nuca,dArrayT alat,
	   int which_rot,dArray2DT mat_rot,
	   double angle) : 
  CrystalLatticeT(nlsd,nuca,which_rot,mat_rot,angle)
{
  for (int i=0; i<nlsd; i++) 
    vLatticeParameters[i] = alat[i];

  // Define lattice type
  sLATTYPE = "CUB";

  if (nlsd==2)
    {
      if(nuca != 1) {cout << "Wrong nuca\n"; throw eSizeMismatch;}

      double Pi = 4.*atan(1.0);
      double Pio3 = Pi/3.;

      // Define basis vectors

      vBasis(0,0) = 0.0;
      vBasis(1,0) = 0.0;

      vAxis(0,0) = vLatticeParameters[0];
      vAxis(1,0) = 0.0;

      vAxis(0,1) = 0.0;
      vAxis(1,1) = vLatticeParameters[1];


      // Rotate axis if necessary
      if (fabs(angle_rotation) >=1.e-5) 
      	vAxis = AxisRotation(vAxis);
    }
  
  if (nlsd==3) 
    {

      if(nuca != 1) {cout << "Wrong nuca\n"; throw eSizeMismatch;}

      vBasis(0,0) = 0.0;
      vBasis(1,0) = 0.0;

      // Define basis vectors
      vAxis(0,0) = vLatticeParameters[0];
      vAxis(1,0) = 0.0;
      vAxis(2,0) = 0.0;

      vAxis(0,1) = 0.0;
      vAxis(1,1) = vLatticeParameters[1];
      vAxis(2,1) = 0.0;

      vAxis(0,2) = 0.0;
      vAxis(1,2) = 0.0;
      vAxis(2,2) = vLatticeParameters[2];

      // Rotate axis if necessary  (put a flag later...)
      double norm = sqrt(norm_vec[0] + norm_vec[1] + norm_vec[2]);
      if (norm > 1.e-5) vAxis = AxisRotation(vAxis);

      // Define primitive lattice vectors
      a[0] = 1.0;a[1]=1.0;a[2]=0.0;
      b[0] = 1.0;b[1]=0.0;b[2]=1.0;
      c[0] = 0.0;c[1]=1.0;c[2]=1.0;
      for (int i=0; i<nLSD; i++)
	{
	  a[i] *= vLatticeParameters[i]*0.5;
	  b[i] *= vLatticeParameters[i]*0.5;
	  c[i] *= vLatticeParameters[i]*0.5;
	}


      if(norm > 1.e-5) 
	{
	  a = VectorRotation(a);
	  b = VectorRotation(b);
	  c = VectorRotation(c);
	}
    }
}


CUBT::CUBT(const CUBT& source) : CrystalLatticeT(source.nLSD,source.nUCA,
						 source.WhichRot,source.matrix_rotation,
						 source.angle_rotation)
{
  nType = source.nType;

  for (int i=0; i<source.nLSD; i++) 
    for (int j=0; j<source.nUCA; j++) 
      vBasis(i,j) = source.vBasis(i,j);

  for (int i=0; i<source.nLSD; i++) 
    for (int j=0; j<source.nLSD; j++) 
      vAxis(i,j) = source.vAxis(i,j);
 
  for (int i=0; i<source.nLSD; i++) 
    vLatticeParameters[i] = source.vLatticeParameters[i];

  for (int i=0; i<source.nUCA; i++) 
    vType[i] = source.vType[i];

}

const dArray2DT& CUBT::GetBasis() 
{
  if (!&vBasis) 
    {
      cout << "vBasis doesn't exist!" << "\n";
      throw eDatabaseFail;
    }
  else 
    {
      return vBasis;		
    }
}

const dArray2DT& CUBT::GetAxis() 
{
  if (!&vAxis) 
    {
      cout << "vAxis doesn't exist!" << "\n";
      throw eDatabaseFail;
    }
  else 
    {
      return vAxis;		
    }
}

const dArrayT& CUBT::GetLatticeParameters() 
{
  if (!&vLatticeParameters) 
    {
      cout << "vLatticeParameters doesn't exist!" << "\n";	
      throw eDatabaseFail;
    }
  else 
    {
      return vLatticeParameters;
    }
}	


