// DEVELOPMENT
/* $Id: CrystalLatticeT.h,v 1.15 2004-08-18 19:53:00 bsun Exp $ */

#ifndef _CRYSTAL_LATTICE_T_H_
#define _CRYSTAL_LATTICE_T_H_

#include "StringT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "ifstreamT.h"

   using namespace Tahoe;

    class CrystalLatticeT {
   
   protected:
   
      int nLSD, nUCA;
      int nType;                   // total number of types
      StringT sLATTYPE;  
   
      dArrayT a,b,c;     // Primitive Lattice Vectors
      dArray2DT vBasis;            // atoms in cell
      dArrayT vLatticeParameters;  // lattice parameters
      dArray2DT vAxis;             // Bravais vectors
      iArrayT vType;               // atom type
   
      int WhichRot;
      dArray2DT matrix_rotation;
      double angle_rotation;
      dArrayT norm_vec;
   
      double density;
   
   public:
   
   // Constructor 
      CrystalLatticeT(int nlsd, int nuca,int which_rot,
      	dArray2DT mat_rot,double angle);
   // Copy Constructor 
      CrystalLatticeT(const CrystalLatticeT& source);
   // Destructor
       ~CrystalLatticeT() { }
   
       int GetNLSD() { 
         return nLSD; }
       int GetNUCA() { 
         return nUCA; }
       int GetNTYPE() { 
         return nType; }
       StringT GetSLATTYPE() { 
         return sLATTYPE; }
   
       int GetRotMeth() { 
         return WhichRot; };
       double GetAngleRotation() { 
         return angle_rotation; };
       dArray2DT GetMatrixRotation() { 
         return matrix_rotation;};
       iArrayT GetType() { 
         return vType; }
       int GetNumberOfType() { 
         return nType; }
   
   
      virtual const dArrayT& GetVector_a() = 0;
      virtual const dArrayT& GetVector_b() = 0;
      virtual const dArrayT& GetVector_c() = 0;
   
      virtual const dArrayT& GetLatticeParameters() = 0;
      virtual const dArray2DT& GetBasis() = 0;
      virtual const dArray2DT& GetAxis() = 0;
   
      void   CalculateDensity();
      double GetDensity();
   
      dArray2DT AxisRotation(dArray2DT A);
      dArrayT VectorRotation(dArrayT v);
   };

#endif
