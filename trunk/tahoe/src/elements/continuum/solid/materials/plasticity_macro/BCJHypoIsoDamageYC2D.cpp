/* $Id: BCJHypoIsoDamageYC2D.cpp,v 1.5 2004-07-15 08:29:14 paklein Exp $ */
#include "BCJHypoIsoDamageYC2D.h"

#include "Utils.h"

using namespace Tahoe;

/* spatial dimension of problem */
const int kNSD = 2;

BCJHypoIsoDamageYC2D::BCJHypoIsoDamageYC2D(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("BCJHypoIsoDamageYC_2D"),
  BCJHypoIsoDamageYC3D   (in, support),  
  f2Ds_ij   (kNSD),
  f2Dc_ijkl (dSymMatrixT::NumValues(kNSD))
{

}

const dSymMatrixT& BCJHypoIsoDamageYC2D::s_ij()
{
  // inherited
  const dSymMatrixT& sij = BCJHypoIsoDamageYC3D::s_ij();

  // reduce stress: 3D -> 2D
  f2Ds_ij.ReduceFrom3D(sij);

  return f2Ds_ij;
}

const dMatrixT& BCJHypoIsoDamageYC2D::c_ijkl()
{
  // inherited
  const dMatrixT& cijkl = BCJHypoIsoDamageYC3D::c_ijkl();

  // reduce cijkl: 3D -> 2D
  f2Dc_ijkl.Rank4ReduceFrom3D(cijkl);

  return f2Dc_ijkl;
}

/* describe the parameters needed by the interface */
void BCJHypoIsoDamageYC2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	BCJHypoIsoDamageYC3D::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);
}
