/* $Id: FiniteStrainT.cpp,v 1.1.2.3 2001-06-28 01:24:11 paklein Exp $ */

#include "FiniteStrainT.h"
#include "ShapeFunctionT.h"

/* constructor */
FiniteStrainT::FiniteStrainT(FEManagerT& fe_manager):
	ElasticT(fe_manager)
{

}

/* called immediately after constructor */
void FiniteStrainT::Initialize(void)
{
	/* inherited */
	ElasticT::Initialize();

	/* allocate deformation gradient list */
	fF_List.Allocate(NumIP());
	for (int i = 0; i < NumIP(); i++)
		fF_List.Allocate(NumSD());

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
