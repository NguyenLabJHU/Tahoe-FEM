/* $Id: GradCrystalPlast2D.cpp,v 1.7 2004-09-10 22:39:43 paklein Exp $ */
#include "GradCrystalPlast2D.h"
#include "Utils.h"

using namespace Tahoe;

/* spatial dimensions of the problem */
const int kNSD = 2;

GradCrystalPlast2D::GradCrystalPlast2D(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("gradient_crystal_plasticity_2D"),
  GradCrystalPlast (in, support),  
  f2Ds_ij    (kNSD),
  f2Dc_ijkl  (dSymMatrixT::NumValues(kNSD))
{
	/* reset default value */
	fConstraint = kPlaneStrain;
}

const dSymMatrixT& GradCrystalPlast2D::s_ij()
{
  // inherited
  const dSymMatrixT& s_ij = GradCrystalPlast::s_ij();

  // reduce savg_ij: 3D -> 2D
  f2Ds_ij.ReduceFrom3D(s_ij);

  return f2Ds_ij;
}

const dMatrixT& GradCrystalPlast2D::c_ijkl()
{
  // inherited
  const dMatrixT& c_ijkl = GradCrystalPlast::c_ijkl();

  // reduce c_ijkl: 3D -> 2D
  f2Dc_ijkl.Rank4ReduceFrom3D(c_ijkl);

  return f2Dc_ijkl;
}
