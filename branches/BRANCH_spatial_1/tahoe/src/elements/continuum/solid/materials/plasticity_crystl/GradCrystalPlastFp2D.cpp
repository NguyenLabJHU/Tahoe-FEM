/* $Id: GradCrystalPlastFp2D.cpp,v 1.6 2004-09-10 22:39:43 paklein Exp $ */
#include "GradCrystalPlastFp2D.h"
#include "Utils.h"

using namespace Tahoe;

/* spatial dimensions of the problem */
const int kNSD = 2;

GradCrystalPlastFp2D::GradCrystalPlastFp2D(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("gradient_crystal_plasticity_Fp_2D"),
  GradCrystalPlastFp (in, support),  
  f2Ds_ij    (kNSD),
  f2Dc_ijkl  (dSymMatrixT::NumValues(kNSD))
{
	/* reset default value */
	fConstraint = kPlaneStrain;
}

const dSymMatrixT& GradCrystalPlastFp2D::s_ij()
{
  // inherited
  const dSymMatrixT& s_ij = GradCrystalPlastFp::s_ij();

  // reduce savg_ij: 3D -> 2D
  f2Ds_ij.ReduceFrom3D(s_ij);

  return f2Ds_ij;
}

const dMatrixT& GradCrystalPlastFp2D::c_ijkl()
{
  // inherited
  const dMatrixT& c_ijkl = GradCrystalPlastFp::c_ijkl();

  // reduce c_ijkl: 3D -> 2D
  f2Dc_ijkl.Rank4ReduceFrom3D(c_ijkl);

  return f2Dc_ijkl;
}
