/*
  File: LocalCrystalPlast2D.cpp
*/

#include "LocalCrystalPlast2D.h"
#include "ElementCardT.h"
#include "ifstreamT.h"

#include "ElasticT.h"
#include "FEManagerT.h"

/* spatial dimensions of the problem */
const int kNSD = 2;

LocalCrystalPlast2D::LocalCrystalPlast2D(ifstreamT& in, const ElasticT& element) :
  LocalCrystalPlast (in, element),  
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

const dMatrixT& LocalCrystalPlast2D::DeformationGradient(const LocalArrayT& disp)
{ 
  // expand total deformation gradient: 2D -> 3D (plane strain)
  fmatx1.Rank2ExpandFrom2D(F(disp));    // fFtot or fFtot_n
  fmatx1(2, 2) = 1.;

  return fmatx1;
}
