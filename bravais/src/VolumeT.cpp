// DEVELOPMENT
/* $Id: VolumeT.cpp,v 1.13 2003-08-14 23:57:53 saubry Exp $ */
#include "VolumeT.h"

VolumeT::VolumeT(int n) 
{
  nSD = n;
  nATOMS = 0;
  VolType = "none";
}


VolumeT::VolumeT(const VolumeT& source)
{
  nSD = source.nSD;
  nATOMS = source.nATOMS;
}


int VolumeT::GetDimensions()
{
  return nSD;
}

double VolumeT::GetVolume() 
{
  return volume;
}

int VolumeT::GetNumberAtoms() 
{
  return nATOMS;
}

StringT* VolumeT::GetAtomNames() 
{
  return &atom_names;
}

dArray2DT* VolumeT::GetAtomCoordinates() 
{
  return &atom_coord;
}

dArray2DT* VolumeT::GetAtomBounds()
{
  return &atom_bounds;
}

iArrayT* VolumeT::GetAtomNumber()
{
  return &atom_number;
}

iArrayT* VolumeT::GetAtomTypes()
{
  return &atom_types;
}

iArrayT* VolumeT::GetAtomParts()
{
  return &atom_parts;
}

const ArrayT< const iArray2DT * > * VolumeT::GetAtomConnect()
{
  return &atom_connect;
}

const ArrayT< StringT > * VolumeT::GetAtomID()
{
  return &atom_ID;
}
