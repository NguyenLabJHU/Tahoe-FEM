/*
  File: BCJHypo2D.cpp
*/

#include "BCJHypo2D.h"
#include "ifstreamT.h"
#include "Utils.h"

#include "ElasticT.h"

/* spatial dimension of problem */
const int kNSD = 2;

BCJHypo2D::BCJHypo2D(ifstreamT& in, const ElasticT& element) :
  BCJHypo3D   (in, element),  
  Material2DT (in, Material2DT::kPlaneStrain),
  f2Ds_ij   (kNSD),
  f2Dc_ijkl (dSymMatrixT::NumValues(kNSD))
{

}

BCJHypo2D::~BCJHypo2D() {} 

const dSymMatrixT& BCJHypo2D::s_ij()
{
  // inherited
  const dSymMatrixT& sij = BCJHypo3D::s_ij();

  // reduce stress: 3D -> 2D
  f2Ds_ij.ReduceFrom3D(sij);
  f2Ds_ij *= fThickness;

  return f2Ds_ij;
}

const dMatrixT& BCJHypo2D::c_ijkl()
{
  // inherited
  const dMatrixT& cijkl = BCJHypo3D::c_ijkl();

  // reduce cijkl: 3D -> 2D
  f2Dc_ijkl.Rank4ReduceFrom3D(cijkl);
  f2Dc_ijkl *= fThickness;

  return f2Dc_ijkl;
}

void BCJHypo2D::Print(ostream& out) const
{
  // inherited
  BCJHypo3D::Print(out);
  Material2DT::Print(out);
}

void BCJHypo2D::PrintName(ostream& out) const
{
  // inherited
  BCJHypo3D::PrintName(out);

  // output model name
  out << "    Plane Strain\n";
}

const dMatrixT& BCJHypo2D::DeformationGradient(const LocalArrayT& disp)
{
  // expand total deformation gradient: 2D -> 3D (plane strain)
  fmatx1.Rank2ExpandFrom2D(F(disp));    // fFtot or fFtot_n
  fmatx1(2, 2) = 1.;

  return fmatx1;
}
