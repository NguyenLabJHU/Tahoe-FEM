/* $Id: SmallStrainT.cpp,v 1.1.2.2 2001-06-26 07:17:33 paklein Exp $ */

#include "SmallStrainT.h"
#include "ShapeFunctionT.h"

/* constructor */
SmallStrainT::SmallStrainT(FEManagerT& fe_manager):
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

/* total strain */
const dSymMatrixT& SmallStrainT::LinearStrain(void) const
{
bug

//finish implementation of strain
}

/* total strain */
const dSymMatrixT& SmallStrainT::LinearStrain(int ip) const
{

}

/* total strain from the end of the previous time step */
const dSymMatrixT& SmallStrainT::LinearStrain_last(void) const
{

}

/* total strain from the end of the previous time step */
const dSymMatrixT& SmallStrainT::LinearStrain_last(int ip) const
{

}

/* compute field gradients */
void SmallStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const
{

}

/* compute field gradients */
void SmallStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, int ip) const
{

}

/* compute field gradients from the end of the previous time step */
void SmallStrainT::ComputeGradient_last(const LocalArrayT& u, dMatrixT& grad_u) const
{

}

/* compute field gradients from the end of the previous time step */
void SmallStrainT::ComputeGradient_last(const LocalArrayT& u, dMatrixT& grad_u, 
	int ip) const
{

}
