/* $Id: FiniteStrainT.cpp,v 1.1.2.2 2001-06-26 07:17:33 paklein Exp $ */

#include "FiniteStrainT.h"
#include "ShapeFunctionT.h"

/* constructor */
FiniteStrainT::FiniteStrainT(FEManagerT& fe_manager):
	ElasticT(fe_manager),
	fMatrixList(4),
	fIPSet(4)
{
	/* allocate return values */
	for (int i = 0; i < fMatrixList.Length(); i++)
		fMatrixList[i].Allocate(NumSD());

	/* initialize */
	fIPSet = -1;
}

/* total deformation gradient */
const dMatrixT& FiniteStrainT::DeformationGradient(void) const
{
	/* get reference */
	int dex = kF;
	int&  ip_flag = fIPSet[dex];
	dMatrixT& mat = fMatrixList[dex];

	/* recompute */
	if (ip_flag != CurrIP())
	{
		/* displacement gradient */
		fShapes->GradU(fLocDisp, mat);

		/* add identity */
		mat.PlusIdentity();

		/* flag */
		ip_flag = CurrIP();
	}
	return mat;
}

/* total deformation gradient */
const dMatrixT& FiniteStrainT::DeformationGradient(int ip) const
{
	/* get reference */
	int dex = kF_ip;
	int&  ip_flag = fIPSet[dex];
	dMatrixT& mat = fMatrixList[dex];

	/* recompute */
	if (ip_flag != CurrIP())
	{
		/* displacement gradient */
		fShapes->GradU(fLocLastDisp, mat, ip);

		/* add identity */
		mat.PlusIdentity();

		/* flag */
		ip_flag = ip;
	}
	return mat;
}

/* total strain from the end of the previous time step */
const dMatrixT& FiniteStrainT::DeformationGradient_last(void) const
{
	/* get reference */
	int dex = kF_last;
	int&  ip_flag = fIPSet[dex];
	dMatrixT& mat = fMatrixList[dex];

	/* recompute */
	if (ip_flag != CurrIP())
	{
		/* displacement gradient */
		fShapes->GradU(fLocLastDisp, mat);

		/* add identity */
		mat.PlusIdentity();

		/* flag */
		ip_flag = CurrIP();
	}
	return mat;
}

/* total strain from the end of the previous time step */
const dMatrixT& FiniteStrainT::DeformationGradient_last(int ip) const
{
	/* get reference */
	int dex = KF_last_ip;
	int&  ip_flag = fIPSet[dex];
	dMatrixT& mat = fMatrixList[dex];

	/* recompute */
	if (ip_flag != CurrIP())
	{
		/* displacement gradient */
		fShapes->GradU(fLocLastDisp, mat, ip);

		/* add identity */
		mat.PlusIdentity();

		/* flag */
		ip_flag = CurrIP();
	}
	return mat;
}

/* compute field gradients with respect to current coordinates */
void FiniteStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const
{
	/* field gradient */
	fShapes->GradU(u, grad_u);
}

/* compute field gradients with respect to current coordinates */
void FiniteStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, 
	int ip) const
{
	/* field gradient */
	fShapes->GradU(u, grad_u, ip);
}

/* compute field gradients with respect to reference coordinates */
void FiniteStrainT::ComputeGradient_reference(const LocalArrayT& u, 
	dMatrixT& grad_u) const
{
#pragma unused(u)
#pragma unused(grad_u)
	cout << "\n FiniteStrainT::ComputeGradient_reference: not implemented" << endl;
	throw;
}

void FiniteStrainT::ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u, 
	int ip) const
{
#pragma unused(u)
#pragma unused(grad_u)
#pragma unused(ip)
	cout << "\n FiniteStrainT::ComputeGradient_reference: not implemented" << endl;
	throw;
}

/***********************************************************************
* Protected
***********************************************************************/

/* increment current element */
bool FiniteStrainT::NextElement(void)
{
	/* sets flags */
	fIPSet = -1;
	
	/* inherited */
	return ElasticT::NextElement();
}
