/* $Id: BCJHypoIsoDamageKE2D.cpp,v 1.6 2004-09-10 22:39:48 paklein Exp $ */
#include "BCJHypoIsoDamageKE2D.h"
#include "Utils.h"

using namespace Tahoe;

/* spatial dimension of problem */
const int kNSD = 2;

BCJHypoIsoDamageKE2D::BCJHypoIsoDamageKE2D(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("BCJHypoIsoDamageKE_2D"),
  BCJHypoIsoDamageKE3D   (in, support),  
  f2Ds_ij   (kNSD),
  f2Dc_ijkl (dSymMatrixT::NumValues(kNSD))
{
	/* reset default value */
	fConstraint = kPlaneStrain;
}

const dSymMatrixT& BCJHypoIsoDamageKE2D::s_ij()
{
  // inherited
  const dSymMatrixT& sij = BCJHypoIsoDamageKE3D::s_ij();

  // reduce stress: 3D -> 2D
  f2Ds_ij.ReduceFrom3D(sij);

  return f2Ds_ij;
}

const dMatrixT& BCJHypoIsoDamageKE2D::c_ijkl()
{
  // inherited
  const dMatrixT& cijkl = BCJHypoIsoDamageKE3D::c_ijkl();

  // reduce cijkl: 3D -> 2D
  f2Dc_ijkl.Rank4ReduceFrom3D(cijkl);

  return f2Dc_ijkl;
}
