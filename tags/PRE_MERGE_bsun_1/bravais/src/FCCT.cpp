// DEVELOPMENT
/* $Id: FCCT.cpp,v 1.16 2004-02-06 22:00:13 saubry Exp $ */
#include "FCCT.h"
#include "CrystalLatticeT.h"

#include "ifstreamT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"

FCCT::FCCT(int nlsd,int nuca,dArrayT alat,
	   int which_rot,dArray2DT mat_rot,
	   double angle) : 
  CrystalLatticeT(nlsd,nuca,which_rot,mat_rot,angle)
{

  for (int i=0; i<nlsd; i++) 
    vLatticeParameters[i] = alat[i]; 

  // Define lattice type
  sLATTYPE = "FCC";

  if (nlsd==2)
    {
      cout << "Cannot create a 2-dimensional FCC lattice!" << "\n";
      cout << "Try a 2-d HEX lattice instead." << "\n";
      throw eBadInputValue;
    }
  
  if (nlsd==3) 
    {
      if(nuca != 4) {cout << "Wrong nuca\n"; throw eSizeMismatch;}

      // Define atoms in cell
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

      // Rotate axis
      double norm = sqrt(norm_vec[0] + norm_vec[1] + norm_vec[2]);
      if (norm > 1.e-5) vAxis = AxisRotation(vAxis);

      // Define primitive lattice vectors
      a[0] = 1.0;a[1]=1.0;a[2]=0.0;
      b[0] = 0.0;b[1]=1.0;b[2]=1.0;
      c[0] = 1.0;c[1]=0.0;c[2]=1.0;

      if(norm > 1.e-5) 
	{
	  a = VectorRotation(a);
	  b = VectorRotation(b);
	  c = VectorRotation(c);
	}
    }
 
}


FCCT::FCCT(const FCCT& source) : CrystalLatticeT(source.nLSD,source.nUCA,
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

const dArray2DT& FCCT::GetBasis() 
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

const dArray2DT& FCCT::GetAxis() 
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

const dArrayT& FCCT::GetLatticeParameters() 
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










