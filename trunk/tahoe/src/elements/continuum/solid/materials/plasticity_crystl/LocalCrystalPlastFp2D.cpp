/* $Id: LocalCrystalPlastFp2D.cpp,v 1.5 2003-01-29 07:35:04 paklein Exp $ */
#include "LocalCrystalPlastFp2D.h"
#include "ElementCardT.h"
#include "ifstreamT.h"

using namespace Tahoe;

/* spatial dimensions of the problem */
const int kNSD = 2;

LocalCrystalPlastFp2D::LocalCrystalPlastFp2D(ifstreamT& in, const FSMatSupportT& support) :
  LocalCrystalPlastFp (in, support),  
  Material2DT         (in, Material2DT::kPlaneStrain),
  f2Dsavg_ij   (kNSD),
  f2Dcavg_ijkl (dSymMatrixT::NumValues(kNSD))
{
 
}

LocalCrystalPlastFp2D::~LocalCrystalPlastFp2D() {} 

const dSymMatrixT& LocalCrystalPlastFp2D::s_ij()
{
  // inherited
  const dSymMatrixT& savg_ij = LocalCrystalPlastFp::s_ij();

  // reduce savg_ij: 3D -> 2D
  f2Dsavg_ij.ReduceFrom3D(savg_ij);
  f2Dsavg_ij *= fThickness;

  return f2Dsavg_ij;
}

const dMatrixT& LocalCrystalPlastFp2D::c_ijkl()
{
  // inherited
  const dMatrixT& cavg_ijkl = LocalCrystalPlastFp::c_ijkl();

  // reduce cavg_ijkl: 3D -> 2D
  f2Dcavg_ijkl.Rank4ReduceFrom3D(cavg_ijkl);
  f2Dcavg_ijkl *= fThickness;

  return f2Dcavg_ijkl;
}

void LocalCrystalPlastFp2D::Print(ostream& out) const
{
  // inherited
  LocalCrystalPlastFp::Print(out);
  Material2DT::Print(out);
}

void LocalCrystalPlastFp2D::PrintName(ostream& out) const
{
  // inherited
  LocalCrystalPlastFp::PrintName(out);

  // output 2D case name
  out << "    Plane Strain\n";
}
