/* $Id: LocalCrystalPlast2D.cpp,v 1.8 2004-09-10 22:39:43 paklein Exp $ */
#include "LocalCrystalPlast2D.h"
#include "ElementCardT.h"

using namespace Tahoe;

/* spatial dimensions of the problem */
const int kNSD = 2;

LocalCrystalPlast2D::LocalCrystalPlast2D(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("local_crystal_plasticity_2D"),
  LocalCrystalPlast (in, support),  
  f2Dsavg_ij   (kNSD),
  f2Dcavg_ijkl (dSymMatrixT::NumValues(kNSD))
{
	/* reset default value */
	fConstraint = kPlaneStrain;
}

const dSymMatrixT& LocalCrystalPlast2D::s_ij()
{
  // inherited
  const dSymMatrixT& savg_ij = LocalCrystalPlast::s_ij();

  // reduce savg_ij: 3D -> 2D
  f2Dsavg_ij.ReduceFrom3D(savg_ij);

  return f2Dsavg_ij;
}

const dMatrixT& LocalCrystalPlast2D::c_ijkl()
{
  // inherited
  const dMatrixT& cavg_ijkl = LocalCrystalPlast::c_ijkl();

  // reduce cavg_ijkl: 3D -> 2D
  f2Dcavg_ijkl.Rank4ReduceFrom3D(cavg_ijkl);

  return f2Dcavg_ijkl;
}
