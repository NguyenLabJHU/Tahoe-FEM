/* $Id: BCJHypoIsoDamageYC2D.cpp,v 1.4 2003-01-29 07:35:06 paklein Exp $ */
#include "BCJHypoIsoDamageYC2D.h"
#include "ifstreamT.h"
#include "Utils.h"

using namespace Tahoe;

/* spatial dimension of problem */
const int kNSD = 2;

BCJHypoIsoDamageYC2D::BCJHypoIsoDamageYC2D(ifstreamT& in, const FSMatSupportT& support) :
  BCJHypoIsoDamageYC3D   (in, support),  
  Material2DT (in, Material2DT::kPlaneStrain),
  f2Ds_ij   (kNSD),
  f2Dc_ijkl (dSymMatrixT::NumValues(kNSD))
{

}

BCJHypoIsoDamageYC2D::~BCJHypoIsoDamageYC2D() {} 

const dSymMatrixT& BCJHypoIsoDamageYC2D::s_ij()
{
  // inherited
  const dSymMatrixT& sij = BCJHypoIsoDamageYC3D::s_ij();

  // reduce stress: 3D -> 2D
  f2Ds_ij.ReduceFrom3D(sij);
  f2Ds_ij *= fThickness;

  return f2Ds_ij;
}

const dMatrixT& BCJHypoIsoDamageYC2D::c_ijkl()
{
  // inherited
  const dMatrixT& cijkl = BCJHypoIsoDamageYC3D::c_ijkl();

  // reduce cijkl: 3D -> 2D
  f2Dc_ijkl.Rank4ReduceFrom3D(cijkl);
  f2Dc_ijkl *= fThickness;

  return f2Dc_ijkl;
}

void BCJHypoIsoDamageYC2D::Print(ostream& out) const
{
  // inherited
  BCJHypoIsoDamageYC3D::Print(out);
  Material2DT::Print(out);
}

void BCJHypoIsoDamageYC2D::PrintName(ostream& out) const
{
  // inherited
  BCJHypoIsoDamageYC3D::PrintName(out);

  // output model name
  out << "    Plane Strain\n";
}
