// DEVELOPMENT
/* $Id: VolumeT.h,v 1.20 2004-08-18 19:53:00 bsun Exp $ */

#ifndef _VOLUME_T_H_
#define _VOLUME_T_H_

#include <iostream.h>

#include "iArrayT.h"
#include "iArray2DT.h"

#include "dArrayT.h"
#include "dArray2DT.h"
#include "CrystalLatticeT.h"
#include "ifstreamT.h"

   using namespace Tahoe;

    class VolumeT 
   {
   
   protected:
   
      int nSD;
      int nATOMS;
      double volume;
      StringT sLATTYPE;
      StringT VolType;
   
      StringT atom_names;
      iArrayT atom_number;
      ArrayT< StringT > atom_ID;
      dArray2DT atom_coord;
      ArrayT< const iArray2DT * >  atom_connect;
      dArray2DT atom_bounds;
      iArrayT atom_types;
      iArrayT atom_parts;
   
      iArrayT WhichSort;
      iArrayT Map;
   
   public:
   
      VolumeT(int n);
       ~VolumeT() {};
      VolumeT(const VolumeT& source);
   
      int GetDimensions();
      int GetNumberAtoms();
      double GetVolume();
   
      int GetNumberOfAtoms();
      StringT GetTypeOfVolume();
   
      virtual void CreateLattice(CrystalLatticeT* pcl) = 0;
      virtual void SortLattice(CrystalLatticeT* pcl) = 0;
   
      virtual void CalculateBounds(CrystalLatticeT* pcl) = 0;
   
      virtual iArrayT GetNCells() = 0;
      virtual dArray2DT GetLength() = 0; 
   
      StringT*   GetAtomNames();
      const ArrayT< StringT > *   GetAtomID();
      dArray2DT* GetAtomCoordinates();
      const ArrayT< const iArray2DT * > * GetAtomConnect();
      dArray2DT* GetAtomBounds();
      iArrayT*   GetAtomNumber();
      iArrayT*   GetAtomTypes();
      iArrayT*   GetAtomParts();
   
   };

    inline int VolumeT::GetNumberOfAtoms() {
      return nATOMS;};
    inline StringT VolumeT::GetTypeOfVolume() {
      return VolType;};

#endif

