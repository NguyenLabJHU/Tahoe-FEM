/*
  File: GradCrystalPlast2D.cpp
*/

#include "GradCrystalPlast2D.h"
#include "Utils.h"

#include "ifstreamT.h"
//#include "ElasticT.h"
//DEV

//#include "ShapeFunctionT.h"
//DEV


/* spatial dimensions of the problem */
const int kNSD = 2;

GradCrystalPlast2D::GradCrystalPlast2D(ifstreamT& in, const FiniteStrainT& element) :
  GradCrystalPlast (in, element),  
  Material2DT      (in, Material2DT::kPlaneStrain),
  f2Ds_ij    (kNSD),
  f2Dc_ijkl  (dSymMatrixT::NumValues(kNSD))
{

}

GradCrystalPlast2D::~GradCrystalPlast2D() {} 

const dSymMatrixT& GradCrystalPlast2D::s_ij()
{
  // inherited
  const dSymMatrixT& s_ij = GradCrystalPlast::s_ij();

  // reduce savg_ij: 3D -> 2D
  f2Ds_ij.ReduceFrom3D(s_ij);
  f2Ds_ij *= fThickness;

  return f2Ds_ij;
}

const dMatrixT& GradCrystalPlast2D::c_ijkl()
{
  // inherited
  const dMatrixT& c_ijkl = GradCrystalPlast::c_ijkl();

  // reduce c_ijkl: 3D -> 2D
  f2Dc_ijkl.Rank4ReduceFrom3D(c_ijkl);
  f2Dc_ijkl *= fThickness;

  return f2Dc_ijkl;
}

void GradCrystalPlast2D::Print(ostream& out) const
{
  // inherited
  GradCrystalPlast::Print(out);
  Material2DT::Print(out);
}

void GradCrystalPlast2D::PrintName(ostream& out) const
{
  // inherited
  GradCrystalPlast::PrintName(out);

  // output 2D case name
  out << "    Plane Strain\n";
}

//DEV - deprecated
#if 0
const dMatrixT& GradCrystalPlast2D::DeformationGradient(const LocalArrayT& disp)
{
  // 2D displacement gradient at IP
  ShapeFunction().GradU(disp, fGradU, fIP);

  // expand total deformation gradient: 2D -> 3D (plane strain)
  fmatx1.Rank2ExpandFrom2D(FDContinuumT::F(fGradU));    // fFtot or fFtot_n
  fmatx1(2, 2) = 1.;

  return fmatx1;
}
#endif
