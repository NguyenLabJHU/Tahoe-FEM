/* $Id: LocalCrystalPlast2D.cpp,v 1.6 2003-01-29 07:35:04 paklein Exp $ */
#include "LocalCrystalPlast2D.h"
#include "ElementCardT.h"
#include "ifstreamT.h"

using namespace Tahoe;

/* spatial dimensions of the problem */
const int kNSD = 2;

LocalCrystalPlast2D::LocalCrystalPlast2D(ifstreamT& in, const FSMatSupportT& support) :
  LocalCrystalPlast (in, support),  
  Material2DT       (in, Material2DT::kPlaneStrain),
  f2Dsavg_ij   (kNSD),
  f2Dcavg_ijkl (dSymMatrixT::NumValues(kNSD))
{
 
}

LocalCrystalPlast2D::~LocalCrystalPlast2D() {} 

const dSymMatrixT& LocalCrystalPlast2D::s_ij()
{
  // inherited
  const dSymMatrixT& savg_ij = LocalCrystalPlast::s_ij();

  // reduce savg_ij: 3D -> 2D
  f2Dsavg_ij.ReduceFrom3D(savg_ij);
  f2Dsavg_ij *= fThickness;

  return f2Dsavg_ij;
}

const dMatrixT& LocalCrystalPlast2D::c_ijkl()
{
  // inherited
  const dMatrixT& cavg_ijkl = LocalCrystalPlast::c_ijkl();

  // reduce cavg_ijkl: 3D -> 2D
  f2Dcavg_ijkl.Rank4ReduceFrom3D(cavg_ijkl);
  f2Dcavg_ijkl *= fThickness;

  return f2Dcavg_ijkl;
}

void LocalCrystalPlast2D::Print(ostream& out) const
{
  // inherited
  LocalCrystalPlast::Print(out);
  Material2DT::Print(out);
}

void LocalCrystalPlast2D::PrintName(ostream& out) const
{
  // inherited
  LocalCrystalPlast::PrintName(out);

  // output 2D case name
  out << "    Plane Strain\n";
}
