/* $Id: J2Simo3D.cpp,v 1.1 2001-05-01 23:24:51 paklein Exp $ */
/* created: paklein (06/22/1997)                                          */

#include "J2Simo3D.h"
#include "ElasticT.h"
#include "ElementCardT.h"

/* constructor */
J2Simo3D::J2Simo3D(ifstreamT& in, const ElasticT& element):
	SimoIso3D(in, element),
	J2SimoLinHardT(in, NumIP(), Mu()),
	fLocLastDisp(element.LastDisplacements()),
	fRelDisp(LocalArrayT::kDisp, fLocLastDisp.NumberOfNodes(), fLocLastDisp.MinorDim()),
	fFtot(3),
	ffrel(3),
	fF_temp(3)
{
	/* check last displacements */
	if (!fLocLastDisp.IsRegistered() ||
		 fLocLastDisp.MinorDim() != NumDOF())
	{
		cout << "\n J2Simo3D::J2Simo3D: last local displacement vector is invalid" << endl;
		throw eGeneralFail;
	}
}

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT J2Simo3D::TangentType(void) const
{
	return GlobalT::kNonSymmetric;
}

/* update internal variables */
void J2Simo3D::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Update(element);
}

/* reset internal variables to last converged solution */
void J2Simo3D::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Reset(element);
}

/* print parameters */
void J2Simo3D::Print(ostream& out) const
{
	/* inherited */
	SimoIso3D::Print(out);
	J2SimoLinHardT::Print(out);
}

/* modulus */
const dMatrixT& J2Simo3D::c_ijkl(void)
{
	/* Compute F_total and f_relative 3D */
	ComputeGradients();

	int ip = CurrIP();
	ElementCardT& element = CurrentElement();

	/* compute isochoric elastic stretch */
	double J = fFtot.Det();
	const dSymMatrixT& b_els = ElasticStretch(fFtot, ffrel, element, ip);
	
	/* elastic tangent modulus */
	ComputeModuli(J, b_els, fModulus);
	
	/* elastoplastic correction */
	fModulus += ModuliCorrection(element, ip);
	
	return fModulus;
}

/* stress */
const dSymMatrixT& J2Simo3D::s_ij(void)
{
	/* Compute F_total and f_relative 3D */
	ComputeGradients();

	int ip = CurrIP();
	ElementCardT& element = CurrentElement();

	/* compute isochoric elastic stretch */
	double J = fFtot.Det();
	const dSymMatrixT& b_els = ElasticStretch(fFtot, ffrel, element, ip);
	
	/* elastic stress */
	ComputeCauchy(J, b_els, fStress);

	/* modify Cauchy stress (return mapping) */
	fStress += StressCorrection(fFtot, ffrel, element, ip);
	
	return fStress;
}

/* returns the strain energy density for the specified strain */
double J2Simo3D::StrainEnergyDensity(void)
{
	/* Compute F_total and f_relative */
	ComputeGradients();

	/* compute isochoric elastic stretch */
	double J = fFtot.Det();
	const dSymMatrixT& b_els = ElasticStretch(fFtot, ffrel,
		CurrentElement(), CurrIP());

	return ComputeEnergy(J, b_els);
}

/* required parameter flags */
bool J2Simo3D::NeedLastDisp(void) const { return true; }

/***********************************************************************
* Protected
***********************************************************************/

/* print name */
void J2Simo3D::PrintName(ostream& out) const
{
	/* inherited */
	SimoIso3D::PrintName(out);
	J2SimoLinHardT::PrintName(out);
}

/***********************************************************************
* Private
***********************************************************************/

/* compute F_total and f_relative */
void J2Simo3D::ComputeGradients(void)
{
	/* total deformation gradient */
	fFtot = F();

	/* relative deformation gradient */
	fF_temp.Inverse(F(fLocLastDisp));
	ffrel.MultAB(fFtot,fF_temp);
}
