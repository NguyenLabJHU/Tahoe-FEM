/* $Id: CBLatticeT.cpp,v 1.7 2004-07-15 08:26:42 paklein Exp $ */
/* created: paklein (12/02/1996) */
#include "CBLatticeT.h"

using namespace Tahoe;

/* constructor */
CBLatticeT::CBLatticeT(void) { }
	
/* fetch bond component tensor (R_I R_J R_K R_L) in reduced index form */
void CBLatticeT::BondComponentTensor4(int numbond, dMatrixT& matrix) const
{
	/* temp */
	dArrayT temp;
	fBonds.RowAlias(numbond, temp);

	if (matrix.Rows() == 3) /* 3 stress components in 2D */
		BondTensor4_2D(temp, matrix);
	else if (matrix.Rows() == 6) /* 6 stress components in 2D */
		BondTensor4_3D(temp, matrix);
	else if (matrix.Rows() == 1) /* 1D */ {
		double R = fBonds[numbond];
		matrix[0] = R*R*R*R;
	}
	else
		ExceptionT::GeneralFail("CBLatticeT::BondComponentTensor2");
}

/*
* Fetch bond component tensor (R_I R_J) in reduced index form, ie
* this function expects a vector.
*/
void CBLatticeT::BondComponentTensor2(int numbond, dArrayT& vector) const
{
	/* temp */
	dArrayT temp;
	fBonds.RowAlias(numbond, temp);

	if (vector.Length() == 3) /* 3 components in 2D */
		BondTensor2_2D(temp, vector);
	else if (vector.Length() == 6) /* 6 stress components in 2D */
		BondTensor2_3D(temp, vector);
	else if (vector.Length() == 1) /* 1D */{
		double R = fBonds[numbond];
		vector[0] = R*R;
	}
	else
		ExceptionT::GeneralFail("CBLatticeT::BondComponentTensor2");
}

void CBLatticeT::BatchBondComponentTensor2(dArray2DT& comptable) const
{
	if (comptable.MinorDim() == 3) /* 3 components in 2D */
		BatchBondTensor2_2D(comptable);
	else
		BatchBondTensor2_3D(comptable);
}

/**********************************************************************
 * Private
 **********************************************************************/

/* building the bond component tensors */
void CBLatticeT::BondTensor4_2D(const dArrayT& comps, dMatrixT& matrix) const
{
	/* dimension check */
	if (matrix.Rows() != 3 || matrix.Cols() != 3) ExceptionT::GeneralFail("CBLatticeT::BondTensor4_2D");

	double R0 = comps[0];
	double R1 = comps[1];

	double m[] = {R0*R0, R1*R1, R0*R1};

	for (int j = 0; j < 3; j++)
		for (int i = 0; i <= j; i++)
			matrix(i,j) = m[i]*m[j];

	for (int j = 0; j < 3; j++) {
		double* col = matrix(j);
		*col++ = m[0]*m[j];
		*col++ = m[1]*m[j];
		*col   = m[2]*m[j];			
	}
}	

void CBLatticeT::BondTensor4_3D(const dArrayT& comps, dMatrixT& matrix) const
{
	/* dimension check */
	if (matrix.Rows() != 6 || matrix.Cols() != 6) ExceptionT::GeneralFail("CBLatticeT::BondTensor4_3D");

	double R0 = comps[0];
	double R1 = comps[1];
	double R2 = comps[2];

	double m[] = {R0*R0, R1*R1, R2*R2, R1*R2, R0*R2, R0*R1};

	for (int j = 0; j < 6; j++)
		for (int i = 0; i <= j; i++)
			matrix(i,j) = m[i]*m[j];

	for (int j = 0; j < 6; j++) {
		double* col = matrix(j);
		*col++ = m[0]*m[j];
		*col++ = m[1]*m[j];
		*col++ = m[2]*m[j];
		*col++ = m[3]*m[j];
		*col++ = m[4]*m[j];
		*col   = m[5]*m[j];
	}
}	

void CBLatticeT::BondTensor2_2D(const dArrayT& comps, dArrayT& vector) const
{	
	/* dimension check */
	if (vector.Length() != 3) ExceptionT::GeneralFail("CBLatticeT::BondTensor2_2D");
	
	double R0 = comps[0];
	double R1 = comps[1];

	vector[0] = R0*R0;
	vector[1] = R1*R1;
	vector[2] = R0*R1;
}
	
void CBLatticeT::BondTensor2_3D(const dArrayT& comps, dArrayT& vector) const
{
	/* dimension check */
	if (vector.Length() != 6) ExceptionT::GeneralFail("CBLatticeT::BondTensor2_3D");
	
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
	if (comptable.MinorDim() != 3) ExceptionT::GeneralFail("CBLatticeT::BatchBondTensor2_2D");

	for (int i = 0; i < fBonds.MajorDim(); i++)
	{
		const double* pbond = fBonds(i);
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
	if (comptable.MinorDim() != 6) ExceptionT::GeneralFail("CBLatticeT::BatchBondTensor2_3D");

	for (int i = 0; i < fBonds.MajorDim(); i++)
	{
		const double* pbond = fBonds(i);
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
