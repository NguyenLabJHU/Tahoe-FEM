/* $Id: GradCrystalPlastFp2D.cpp,v 1.4 2003-01-29 07:35:04 paklein Exp $ */
#include "GradCrystalPlastFp2D.h"
#include "Utils.h"

#include "ifstreamT.h"

using namespace Tahoe;

/* spatial dimensions of the problem */
const int kNSD = 2;

GradCrystalPlastFp2D::GradCrystalPlastFp2D(ifstreamT& in, const FSMatSupportT& support) :
  GradCrystalPlastFp (in, support),  
  Material2DT        (in, Material2DT::kPlaneStrain),
  f2Ds_ij    (kNSD),
  f2Dc_ijkl  (dSymMatrixT::NumValues(kNSD))
{

}

GradCrystalPlastFp2D::~GradCrystalPlastFp2D() {} 

const dSymMatrixT& GradCrystalPlastFp2D::s_ij()
{
  // inherited
  const dSymMatrixT& s_ij = GradCrystalPlastFp::s_ij();

  // reduce savg_ij: 3D -> 2D
  f2Ds_ij.ReduceFrom3D(s_ij);
  f2Ds_ij *= fThickness;

  return f2Ds_ij;
}

const dMatrixT& GradCrystalPlastFp2D::c_ijkl()
{
  // inherited
  const dMatrixT& c_ijkl = GradCrystalPlastFp::c_ijkl();

  // reduce c_ijkl: 3D -> 2D
  f2Dc_ijkl.Rank4ReduceFrom3D(c_ijkl);
  f2Dc_ijkl *= fThickness;

  return f2Dc_ijkl;
}

void GradCrystalPlastFp2D::Print(ostream& out) const
{
  // inherited
  GradCrystalPlastFp::Print(out);
  Material2DT::Print(out);
}

void GradCrystalPlastFp2D::PrintName(ostream& out) const
{
  // inherited
  GradCrystalPlastFp::PrintName(out);

  // output 2D case name
  out << "    Plane Strain\n";
}
