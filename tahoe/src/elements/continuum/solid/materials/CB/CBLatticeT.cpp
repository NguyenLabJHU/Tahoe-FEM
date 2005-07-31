/* $Id: CBLatticeT.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (12/02/1996)                                          */
/* CBLatticeT.cpp                                                         */

#include "CBLatticeT.h"

/* Constructor */
CBLatticeT::CBLatticeT(int numlatticedim, int numspatialdim,
	int numbonds): BondLatticeT(numlatticedim, numspatialdim, numbonds)
{

}

/*
* The Q matrix passed into this constructor is used to rotate the
* bond vectors into the orientation prescribed by Q.  No check is
* performed on the orthogonality of Q, only its dimensions.  Q is
* deep copied.  Q is defined as:
*
*			Q = d x_global / d x_natural
*
* So that the vectors are transformed by:
*
*			r_new = Q.r_natural
*
*/
CBLatticeT::CBLatticeT(const dMatrixT& Q, int numspatialdim,
	int numbonds): BondLatticeT(Q, numspatialdim, numbonds)
{

}
	
/*
* Fetch bond component tensor (R_I R_J R_K R_L) in reduced index
* form.
*/
void CBLatticeT::BondComponentTensor4(int numbond, dMatrixT& matrix) const
{
	/* temp */
	dArrayT temp(fNumSpatialDim, fBonds(numbond));

	if (fNumSpatialDim == 2)
		BondTensor4_2D(temp, matrix);
	else
		BondTensor4_3D(temp, matrix);
}

/*
* Fetch bond component tensor (R_I R_J) in reduced index form, ie
* this function expects a vector.
*/
void CBLatticeT::BondComponentTensor2(int numbond, dArrayT& vector) const
{
	/* wrap */
	dArrayT temp(fNumSpatialDim, fBonds(numbond));

	if (fNumSpatialDim == 2)
		BondTensor2_2D(temp, vector);
	else
		BondTensor2_3D(temp, vector);
}

void CBLatticeT::BatchBondComponentTensor2(dArray2DT& comptable) const
{
	if (fNumSpatialDim == 2)
		BatchBondTensor2_2D(comptable);
	else
		BatchBondTensor2_3D(comptable);
}

/**********************************************************************
* Private
**********************************************************************/

/*
* Building the bond component tensors.
*/
void CBLatticeT::BondTensor4_2D(const dArrayT& comps, dMatrixT& matrix) const
{
	/* dimension check */
	if (matrix.Rows() != 3 || matrix.Cols() != 3) throw eGeneralFail;

	double R0 = comps[0];
	double R1 = comps[1];

	double m[] = {R0*R0, R1*R1, R0*R1};

	for (int j = 0; j < 3; j++)
		for (int i = 0; i <= j; i++)
			matrix(i,j) = m[i]*m[j];

	matrix.CopySymmetric();
}	

void CBLatticeT::BondTensor4_3D(const dArrayT& comps, dMatrixT& matrix) const
{
	/* dimension check */
	if (matrix.Rows() != 6 || matrix.Cols() != 6) throw eGeneralFail;

	double R0 = comps[0];
	double R1 = comps[1];
	double R2 = comps[2];

	double m[] = {R0*R0, R1*R1, R2*R2, R1*R2, R0*R2, R0*R1};

	for (int j = 0; j < 6; j++)
		for (int i = 0; i <= j; i++)
			matrix(i,j) = m[i]*m[j];

	matrix.CopySymmetric();
}	

void CBLatticeT::BondTensor2_2D(const dArrayT& comps, dArrayT& vector) const
{	
	/* dimension check */
	if (vector.Length() != 3) throw eGeneralFail;
	
	double R0 = comps[0];
	double R1 = comps[1];

	vector[0] = R0*R0;
	vector[1] = R1*R1;
	vector[2] = R0*R1;
}
	
void CBLatticeT::BondTensor2_3D(const dArrayT& comps, dArrayT& vector) const
{
	/* dimension check */
	if (vector.Length() != 6) throw eGeneralFail;
	
	double R0 = comps[0];
	double R1 = comps[1];
	double R2 = comps[2];

	vector[0] = R0*R0;
	vector[1] = R1*R1;
	vector[2] = R2*R2;
	vector[3] = R1*R2;
	vector[4] = R0*R2;
	vector[5] = R0*R1;	
}	

/* batched versions */	
void CBLatticeT::BatchBondTensor2_2D(dArray2DT& comptable) const
{
	/* dimension check */
	if (comptable.MinorDim() != 3) throw eGeneralFail;

	for (int i = 0; i < fNumBonds; i++)
	{
		double* pbond = fBonds(i);
		double* pcomp = comptable(i);
	
		double R0 = pbond[0];
		double R1 = pbond[1];

		pcomp[0] = R0*R0;
		pcomp[1] = R1*R1;
		pcomp[2] = R0*R1;
	}
}	

void CBLatticeT::BatchBondTensor2_3D(dArray2DT& comptable) const
{
	/* dimension check */
	if (comptable.MinorDim() != 6) throw eGeneralFail;

	for (int i = 0; i < fNumBonds; i++)
	{
		double* pbond = fBonds(i);
		double* pcomp = comptable(i);
	
		double R0 = pbond[0];
		double R1 = pbond[1];
		double R2 = pbond[2];
	
		pcomp[0] = R0*R0;
		pcomp[1] = R1*R1;
		pcomp[2] = R2*R2;
		pcomp[3] = R1*R2;
		pcomp[4] = R0*R2;
		pcomp[5] = R0*R1;	
	}
}
