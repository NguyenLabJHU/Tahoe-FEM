/* $Id: J2Simo2D.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/22/1997)                                          */

#include "J2Simo2D.h"
#include "ElasticT.h"
#include "ElementCardT.h"

/* constructor */
J2Simo2D::J2Simo2D(ifstreamT& in, const ElasticT& element):
	SimoIso2D(in, element),
	J2SimoLinHardT(in, NumIP(), Mu()),
	fLocLastDisp(element.LastDisplacements()),
	fRelDisp(LocalArrayT::kDisp, fLocLastDisp.NumberOfNodes(), fLocLastDisp.MinorDim()),
	fFtot(3),
	ffrel(3),
	fF_temp(2),
	fFtot_2D(2),
	ffrel_2D(2)	
{
	/* check last displacements */
	if (!fLocLastDisp.IsRegistered() ||
		 fLocLastDisp.MinorDim() != NumDOF())
	{
		cout << "\n J2Simo2D::J2Simo2D: last local displacement vector is invalid" << endl;
		throw eGeneralFail;
	}
}

/* form of tangent matrix (symmetric by default) */
GlobalT::SystemTypeT J2Simo2D::TangentType(void) const
{
	return GlobalT::kNonSymmetric;
}

/* update internal variables */
void J2Simo2D::UpdateHistory(void)
{
	/* update if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Update(element);
}

/* reset internal variables to last converged solution */
void J2Simo2D::ResetHistory(void)
{
	/* reset if plastic */
	ElementCardT& element = CurrentElement();
	if (element.IsAllocated()) Reset(element);
}

/* print parameters */
void J2Simo2D::Print(ostream& out) const
{
	/* inherited */
	SimoIso2D::Print(out);
	Material2DT::Print(out);
	J2SimoLinHardT::Print(out);
}

/* modulus */
const dMatrixT& J2Simo2D::c_ijkl(void)
{
	/* Compute F_total and f_relative 3D */
	ComputeGradients();

	int ip = CurrIP();
	ElementCardT& element = CurrentElement();

	/* compute isochoric elastic stretch */
	double J = fFtot.Det();
	const dSymMatrixT& b_els = ElasticStretch(fFtot, ffrel, element, ip);
	
	/* 3D elastic stress */
	ComputeModuli(J, b_els, fModulus);
	
	/* elastoplastic correction */
	fModulus += ModuliCorrection(element, ip);
	
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(fModulus);
	fModulus2D *= fThickness;

	return fModulus2D;
}

/* stress */
const dSymMatrixT& J2Simo2D::s_ij(void)
{
	/* Compute F_total and f_relative 3D */
	ComputeGradients();

	int ip = CurrIP();
	ElementCardT& element = CurrentElement();

	/* compute isochoric elastic stretch */
	double J = fFtot.Det();
	const dSymMatrixT& b_els = ElasticStretch(fFtot, ffrel, element, ip);
	
	/* 3D elastic stress */
	ComputeCauchy(J, b_els, fStress);

	/* modify Cauchy stress (return mapping) */
	fStress += StressCorrection(fFtot, ffrel, element, ip);
	
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(fStress);
	fStress2D *= fThickness;

	return fStress2D;
}

/* returns the strain energy density for the specified strain */
double J2Simo2D::StrainEnergyDensity(void)
{
	/* Compute F_total and f_relative 3D */
	ComputeGradients();

	/* compute isochoric elastic stretch */
	double J = fFtot.Det();
	const dSymMatrixT& b_els = ElasticStretch(fFtot, ffrel,
		CurrentElement(), CurrIP());

	return fThickness*ComputeEnergy(J, b_els);
}

/* required parameter flags */
bool J2Simo2D::NeedLastDisp(void) const { return true; }

/***********************************************************************
* Protected
***********************************************************************/

/* print name */
void J2Simo2D::PrintName(ostream& out) const
{
	/* inherited */
	SimoIso2D::PrintName(out);
	J2SimoLinHardT::PrintName(out);
}

/***********************************************************************
* Private
***********************************************************************/

/* compute F_total and f_relative */
void J2Simo2D::ComputeGradients(void)
{
	/* compute relative deformation gradient */
	fFtot_2D = F();
	fF_temp.Inverse( F(fLocLastDisp) );
	ffrel_2D.MultAB(fFtot_2D,fF_temp);

	/* 2D -> 3D */
	fFtot.Rank2ExpandFrom2D(fFtot_2D);
	fFtot(2,2) = 1.0;

	ffrel.Rank2ExpandFrom2D(ffrel_2D);
	ffrel(2,2) = 1.0;
	
}
