// DEVELOPMENT
/* $Id: VolumeT.cpp,v 1.11 2003-08-01 22:54:39 saubry Exp $ */
#include "VolumeT.h"

VolumeT::VolumeT(int n) 
{
  nSD = n;
  nATOMS = 0;
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

iArrayT* VolumeT::GetAtomID()
{
  return &atom_ID;
}

dArray2DT* VolumeT::GetAtomCoordinates() 
{
  return &atom_coord;
}

iArray2DT* VolumeT::GetAtomConnectivities()
{
  return &atom_connectivities;
}

dArray2DT* VolumeT::GetAtomBounds()
{
  return &atom_bounds;
}

iArrayT* VolumeT::GetAtomTypes()
{
  return &atom_types;
}

iArrayT* VolumeT::GetAtomParts()
{
  return &atom_parts;
}

const ArrayT< const iArray2DT * > * VolumeT::GetAtomArrayConnect()
{
  return &atom_array_connect;
}

const ArrayT< StringT > * VolumeT::GetAtomArrayID()
{
  return &atom_array_ID;
}
