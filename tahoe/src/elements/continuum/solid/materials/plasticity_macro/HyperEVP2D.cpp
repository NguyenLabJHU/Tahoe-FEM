/* $Id: HyperEVP2D.cpp,v 1.5.30.2 2004-03-03 16:15:06 paklein Exp $ */
#include "HyperEVP2D.h"
#include "ifstreamT.h"
#include "Utils.h"

using namespace Tahoe;

/* spatial dimension of problem */
const int kNSD = 2;

HyperEVP2D::HyperEVP2D(ifstreamT& in, const FSMatSupportT& support) :
	ParameterInterfaceT("HyperEVP_2D"),
  HyperEVP3D  (in, support),  
  f2Ds_ij   (kNSD),
  f2Dc_ijkl (dSymMatrixT::NumValues(kNSD))
{

}

const dSymMatrixT& HyperEVP2D::s_ij()
{
  // inherited
  const dSymMatrixT& sij = HyperEVP3D::s_ij();

  // reduce stress: 3D -> 2D
  f2Ds_ij.ReduceFrom3D(sij);

  return f2Ds_ij;
}

const dMatrixT& HyperEVP2D::c_ijkl()
{
  // inherited
  const dMatrixT& cijkl = HyperEVP3D::c_ijkl();

  // reduce cijkl: 3D -> 2D
  f2Dc_ijkl.Rank4ReduceFrom3D(cijkl);

  return f2Dc_ijkl;
}

void HyperEVP2D::PrintName(ostream& out) const
{
  // inherited
  HyperEVP3D::PrintName(out);

  // output model name
  out << "    Plane Strain\n";
}

/* describe the parameters needed by the interface */
void HyperEVP2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	HyperEVP3D::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("2D_constraint");
	constraint.SetDefault(kPlaneStrain);
}
