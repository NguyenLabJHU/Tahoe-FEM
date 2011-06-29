/* $Id: LocalCrystalPlast2D.cpp,v 1.9 2005-01-21 16:51:21 paklein Exp $ */
#include "LocalCrystalPlast2D.h"
#include "ElementCardT.h"

using namespace Tahoe;

/* spatial dimensions of the problem */
const int kNSD = 2;

LocalCrystalPlast2D::LocalCrystalPlast2D(void):
	ParameterInterfaceT("local_crystal_plasticity_2D")
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

/* accept parameter list */
void LocalCrystalPlast2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	LocalCrystalPlast::TakeParameterList(list);

	/* dimension work space */
	f2Dsavg_ij.Dimension(kNSD);
	f2Dcavg_ijkl.Dimension(dSymMatrixT::NumValues(kNSD));
}
