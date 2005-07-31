/* $Id: MLSSolverT.cpp,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (12/08/1999)                                          */

#include "MLSSolverT.h"

#include "ExceptionCodes.h"
#include "dSymMatrixT.h"

/* basis functions */
#include "PolyBasis1DT.h"
#include "PolyBasis2DT.h"
#include "PolyBasis3DT.h"

/* C1 window functions */
#include "C1FunctionT.h"

/* constants */
const double sqrtPi = sqrt(acos(-1.0));

/* constructor */
MLSSolverT::MLSSolverT(int nsd, int complete):
	fNumSD(nsd),
	fComplete(complete),
	fNumNeighbors(0),
	fOrder(0),
	fDb(fNumSD),
	fDDb(dSymMatrixT::NumValues(fNumSD)),
	fDM(fNumSD),
	fDDM(dSymMatrixT::NumValues(fNumSD)),
	/* variable memory managers */
	fArrayGroup(0),
	fArray2DGroup2(0, 0),
	fArray2DGroup3(0, 0),
	fLocCoords_man(0, fLocCoords, fNumSD),
	
	/* work space */
	fNSDsym(fNumSD)
{
	/* basis functions */
	switch (fNumSD)
	{
		case 1:
			fBasis = new PolyBasis1DT(fComplete);
			break;
		case 2:
			fBasis = new PolyBasis2DT(fComplete);
			break;
		case 3:
			fBasis = new PolyBasis3DT(fComplete);
			break;
		default:
			cout << "\n MLSSolverT::MLSSolverT: unsupported spatial dimensions "
			     << fNumSD << endl;
			throw eBadInputValue;
	}
	if (!fBasis) throw eOutOfMemory;

	/* construct window function */
	fWindow = NULL;
	
	/* dimension arrays */
	int m = fBasis->BasisDimension();
	fb.Allocate(m);
	for (int j = 0; j < fDb.Length(); j++)
		fDb[j].Allocate(m);
	for (int jj = 0; jj < fDDb.Length(); jj++)
		fDDb[jj].Allocate(m);

	fMinv.Allocate(m);
	for (int i = 0; i < fDM.Length(); i++)
		fDM[i].Allocate(m);
	for (int ii = 0; ii < fDDM.Length(); ii++)
		fDDM[ii].Allocate(m);

	/* work space */
	fMtemp.Allocate(m);
	fbtemp1.Allocate(m);
	fbtemp2.Allocate(m);
	fbtemp3.Allocate(m);
	fbtemp4.Allocate(m);
}
	
/* destructor */
MLSSolverT::~MLSSolverT(void)
{
	delete fBasis;
	delete fWindow;
}

/* class dependent initializations */
void MLSSolverT::Initialize(void)
{
	/* dimensions */
//	int   m = fBasis->BasisDimension();
//	int dim = NumValues(fNumSD);

	/* window function */
	fArrayGroup.Register(fw);
	fArray2DGroup2.Register(fDw);		
	fArray2DGroup3.Register(fDDw);		

	/* correction function */
	fArrayGroup.Register(fC);
	fArray2DGroup2.Register(fDC);		
	fArray2DGroup3.Register(fDDC);		

	/* shape function */
	fArrayGroup.Register(fphi);
	fArray2DGroup2.Register(fDphi);		
	fArray2DGroup3.Register(fDDphi);		
}

/* set MLS at coords given sampling points */
int MLSSolverT::SetField(const dArray2DT& coords, const dArrayT& dmax,
		const dArrayT& volume, const dArrayT& fieldpt, int order)
{
#if __option(extended_errorcheck)
	if (dmax.Length()    != coords.MajorDim()) throw eSizeMismatch;
	if (volume.Length()  != coords.MajorDim()) throw eSizeMismatch;
	if (fieldpt.Length() != fNumSD) throw eSizeMismatch;
	if (order < 0 || order > 2) throw eOutOfRange;
#endif
	
	/* dimension work space */
	fOrder = order;
	fNumNeighbors = coords.MajorDim();
	Dimension();

	/* convolution coordinates */
	fLocCoords.SetToScaled(-1.0, coords);
	fLocCoords.AddToRowsScaled(1.0, fieldpt);

	/* window functions */
	int numactive = SetWindow(dmax);
	if (numactive < fBasis->BasisDimension())
	{
		cout << "\n MLSSolverT::SetField: not enough nodes for fit: ";
		cout << numactive << '/' << fBasis->BasisDimension() << '\n';
		return 0;
	}

	/* evaluate basis functions at coords */
	fBasis->SetBasis(fLocCoords, fOrder);
	
	/* set moment matrix, inverse, and derivatives */
	if (!SetMomentMartrix(volume))
	{
		cout << "\n MLSSolverT::SetField: error in momentum matrix: ";
		return 0;
	}	
	
	/* set correction function coefficients */
	SetCorrectionCoefficient();
	
	/* set correction function */
	SetCorrection();

	/* set shape functions */
	SetShapeFunctions(volume);

	/* OK */
	return 1;
}

/* basis dimension */
int MLSSolverT::BasisDimension(void) const
{
	return fBasis->BasisDimension();
}

/***********************************************************************
* Protected
***********************************************************************/

/***********************************************************************
* Private
***********************************************************************/

/* configure solver for current number of neighbors */
void MLSSolverT::Dimension(void)
{
	/* set variable memory */
	fArrayGroup.Dimension(fNumNeighbors, false);
	fLocCoords_man.SetMajorDimension(fNumNeighbors, false);
	if (fOrder > 0) fArray2DGroup2.Dimension(fNumSD, fNumNeighbors);
	if (fOrder > 1) fArray2DGroup3.Dimension(dSymMatrixT::NumValues(fNumSD),
		fNumNeighbors);
}

/* set window functions */
int MLSSolverT::SetWindow(const dArrayT& dmax)
{
	dArrayT dx;
	dArrayT Dw(fNumSD); //make into work space variables
	int count = 0;
	for (int i = 0; i < fNumNeighbors; i++)
	{
		/* fetch local coords */
		fLocCoords.RowAlias(i, dx);
		
		double dm = dmax[i];
		double di = dx.Magnitude();
	
		/* out of influence range (or inactive) */
		if (di > 4.0*dm)
		{
			 fw[i] = 0.0;
			 if (fOrder > 0)
			 {
			 	fDw.SetColumn(i, 0.0);
			 	
			 	if (fOrder > 1)
			 		fDDw.SetColumn(i, 0.0);
			 }
		}
		/* Gaussian window function and derivatives */
		else
		{
			double    a = 0.4;
//			double    a = 1.0;
			double  adm = a*dm;
			double adm2 = adm*adm;
			double    q = di/adm;
			double&   w = fw[i];
			
			/* window */
	w = exp(-q*q)/(sqrtPi*adm);

			if (fOrder > 0)
			{
				/* 1st derivative */
				Dw.SetToScaled(-2.0*w/adm2, dx);
				fDw.SetColumn(i, Dw);

				if (fOrder > 1)
				{
					/* 2nd derivative */
					fNSDsym.Outer(dx);
					fNSDsym *= 4.0*w/(adm2*adm2);
					fNSDsym.PlusIdentity(-2.0*w/adm2);
					fDDw.SetColumn(i, fNSDsym);
				}
			}

			count++;
		}
	}
	
	return count;
}

/* set moment matrix, inverse, and derivatives */
int MLSSolverT::SetMomentMartrix(const dArrayT& volume)
{
	/* moment matrix */
	ComputeM(volume);
	if (fOrder > 0)
	{
		ComputeDM(volume);
		if (fOrder > 1)
			ComputeDDM(volume);
	}

	/* compute inverse */
	int OK = 1;
	switch (fMinv.Rows())
	{
		case 1:
			fMinv(0,0) = 1.0/fMinv(0,0);
			break;
		case 2:
		{
			double m00 = fMinv(0,0);
			double m11 = fMinv(1,1);
			double m01 = fMinv(0,1);

			double det = m00*m11 - m01*m01;
			if (fabs(det) < kSmall) OK = 0;

			fMinv(0,0) = m11/det;
			fMinv(1,1) = m00/det;
			fMinv(0,1) = fMinv(1,0) = -m01/det;
			break;
		}
		case 3:
			OK = SymmetricInverse3x3(fMinv);
			break;
		case 4:
			OK = SymmetricInverse4x4(fMinv);
			break;
		default:
			throw eGeneralFail;
	}
	return OK;
}

void MLSSolverT::ComputeM(const dArrayT& volume)
{
	int dim = fBasis->BasisDimension();

	//TEMP - use general routine for now
	if (1)
	{
		const dArray2DT& basis = fBasis->P();
		for (int i = 0; i < dim; i++)
			for (int j = i; j < dim; j++)
			{
				double*   w = fw.Pointer();
				double*  pi = basis(i);
				double*  pj = basis(j);
				double* vol = volume.Pointer();
				double  mij = 0.0;
				for (int k = 0; k < fNumNeighbors; k++)
					mij += (*pi++)*(*w++)*(*pj++)*(*vol++);
				fMinv(i,j) = fMinv(j,i) = mij;
			}
	}

	return;

	if (dim == 2)
	{
		const dArray2DT& basis = fBasis->P();
		double*  p0 = basis(0);
		double*  p1 = basis(1);
		double*   w = fw.Pointer();
		double* vol = volume.Pointer();

		/* integrate */
		double m00 = 0.0;
		double m11 = 0.0;
		double m01 = 0.0;	
		for (int i = 0; i < fNumNeighbors; i++)
		{
			double wxvol = (*w++)*(*vol++);
			m00 += (*p0)*wxvol*(*p0);
			m11 += (*p1)*wxvol*(*p1);
			m01 += (*p0)*wxvol*(*p1);
			
			p0++; p1++;		
		}
		
		fMinv(0,0) = m00;
		fMinv(1,1) = m11;
		fMinv(1,0) = fMinv(0,1) = m01;
	}
	if (dim == 3)
	{
		const dArray2DT& basis = fBasis->P();
		double*  p0 = basis(0);
		double*  p1 = basis(1);
		double*  p2 = basis(2);
		double*   w = fw.Pointer();
		double* vol = volume.Pointer();

		/* integrate */
		double m00 = 0.0;
		double m11 = 0.0;
		double m22 = 0.0;
		double m12 = 0.0;	
		double m02 = 0.0;	
		double m01 = 0.0;	
		for (int i = 0; i < fNumNeighbors; i++)
		{
			double wxvol = (*w++)*(*vol++);
			m00 += (*p0)*wxvol*(*p0);
			m11 += (*p1)*wxvol*(*p1);
			m22 += (*p2)*wxvol*(*p2);
			m12 += (*p1)*wxvol*(*p2);
			m02 += (*p0)*wxvol*(*p2);
			m01 += (*p0)*wxvol*(*p1);
			
			p0++; p1++; p2++;		
		}
		
		fMinv(0,0) = m00;
		fMinv(1,1) = m11;
		fMinv(2,2) = m22;
		fMinv(1,2) = fMinv(2,1) = m12;
		fMinv(2,0) = fMinv(0,2) = m02;
		fMinv(1,0) = fMinv(0,1) = m01;
	}
	if (dim == 4)
	{
		const dArray2DT& basis = fBasis->P();
		double*  p0 = basis(0);
		double*  p1 = basis(1);
		double*  p2 = basis(2);
		double*  p3 = basis(3);
		double*  w = fw.Pointer();
		double* vol = volume.Pointer();

		/* integrate */
		double m00 = 0.0;
		double m11 = 0.0;
		double m22 = 0.0;
		double m33 = 0.0;
		double m23 = 0.0;	
		double m13 = 0.0;	
		double m12 = 0.0;	
		double m03 = 0.0;	
		double m02 = 0.0;	
		double m01 = 0.0;	
		for (int i = 0; i < fNumNeighbors; i++)
		{
			double wxvol = (*w++)*(*vol++);
			m00 += (*p0)*wxvol*(*p0);
			m11 += (*p1)*wxvol*(*p1);
			m22 += (*p2)*wxvol*(*p2);
			m33 += (*p3)*wxvol*(*p3);
			m23 += (*p2)*wxvol*(*p3);
			m12 += (*p1)*wxvol*(*p3);
			m12 += (*p1)*wxvol*(*p2);
			m03 += (*p0)*wxvol*(*p3);
			m02 += (*p0)*wxvol*(*p2);
			m01 += (*p0)*wxvol*(*p1);
			
			p0++; p1++; p2++; p3++;
		}
		
		fMinv(0,0) = m00;
		fMinv(1,1) = m11;
		fMinv(2,2) = m22;
		fMinv(3,3) = m33;
		fMinv(2,3) = fMinv(3,2) = m23;
		fMinv(1,3) = fMinv(3,1) = m13;
		fMinv(1,2) = fMinv(2,1) = m12;
		fMinv(3,0) = fMinv(0,3) = m03;
		fMinv(2,0) = fMinv(0,2) = m02;
		fMinv(1,0) = fMinv(0,1) = m01;
	}
	else
	{
		const dArray2DT& basis = fBasis->P();
		for (int i = 0; i < dim; i++)
			for (int j = i; j < dim; j++)
			{
				double*   w = fw.Pointer();
				double*  pi = basis(i);
				double*  pj = basis(j);
				double* vol = volume.Pointer();
				double  mij = 0.0;
				for (int k = 0; k < fNumNeighbors; k++)
					mij += (*pi++)*(*w++)*(*pj++)*(*vol++);
				fMinv(i,j) = fMinv(j,i) = mij;
			}
	}
}

void MLSSolverT::ComputeDM(const dArrayT& volume)
{
	int dim = fBasis->BasisDimension();
	const dArray2DT& basis = fBasis->P();
	for (int m = 0; m < fNumSD; m++)
	{
		dMatrixT& DM = fDM[m];
		const dArray2DT& Dbasis = fBasis->DP(m);
		for (int i = 0; i < dim; i++)
			for (int j = i; j < dim; j++)
			{
				double*   w = fw.Pointer();
				double*  Dw = fDw(m);
				double*  pi = basis(i);
				double*  pj = basis(j);
				double* Dpi = Dbasis(i);
				double* Dpj = Dbasis(j);
				double* vol = volume.Pointer();
				double Dmij = 0.0;
				for (int k = 0; k < fNumNeighbors; k++)
				{
					Dmij += ((*Dpi)*(*w)*(*pj)
					       + (*pi)*(*Dw)*(*pj)
					       + (*pi)*(*w)*(*Dpj))*(*vol);
					
					w++; Dw++; pi++; pj++; Dpi++; Dpj++; vol++;
				}
				DM(i,j) = DM(j,i) = Dmij;
			}
	}
}

void MLSSolverT::ComputeDDM(const dArrayT& volume)
{
	int dim = fBasis->BasisDimension();
	const dArray2DT& basis = fBasis->P();
	for (int rs = 0; rs < dSymMatrixT::NumValues(fNumSD); rs++)
	{
		dMatrixT& DDM = fDDM[rs];

		/* resolve components */
		int r, s;
		dSymMatrixT::ExpandIndex(fNumSD, rs, r, s);
		
		const dArray2DT& Dbasis_r = fBasis->DP(r);
		const dArray2DT& Dbasis_s = fBasis->DP(s);
		const dArray2DT& DDbasis  = fBasis->DDP(rs);

		for (int i = 0; i < dim; i++)
			for (int j = i; j < dim; j++)
			{
				double*   w = fw.Pointer();
				double*  Dw_r = fDw(r);
				double*  Dw_s = fDw(s);
				double*  DDw  = fDDw(rs);

				double*    pi = basis(i);
				double*    pj = basis(j);
				double* Dpi_r = Dbasis_r(i);
				double* Dpi_s = Dbasis_s(i);
				double* Dpj_r = Dbasis_r(j);
				double* Dpj_s = Dbasis_s(j);
				double*  DDpi = DDbasis(i);
				double*  DDpj = DDbasis(j);

				double* vol = volume.Pointer();
				double DDmij = 0.0;
				for (int k = 0; k < fNumNeighbors; k++)
				{
					DDmij += ((*DDpi)*(*w)*(*pj)
					        + (*pi)*(*DDw)*(*pj)
					        + (*pi)*(*w)*(*DDpj)
					        + (*pi)*((*Dpj_r)*(*Dw_s) + (*Dpj_s)*(*Dw_r))
					        + (*pj)*((*Dpi_r)*(*Dw_s) + (*Dpi_s)*(*Dw_r))
					        + (*w)*((*Dpi_r)*(*Dpj_s) + (*Dpi_s)*(*Dpj_r)))*(*vol);
					
					w++; Dw_r++; Dw_s++; DDw++;
					pi++; pj++;
					Dpi_r++; Dpi_s++;
					Dpj_r++; Dpj_s++;
					DDpi++; DDpj++;
					vol++;
				}
				DDM(i,j) = DDM(j,i) = DDmij;
			}
	}
}

/* set correction function coefficients */
void MLSSolverT::SetCorrectionCoefficient(void)
{
	/* coefficients */
	fMinv.CopyColumn(0, fb);
	if (fOrder > 0)
	{
		/* first derivative */
		for (int i = 0; i < fNumSD; i++)
		{
			fDM[i].Multx(fb, fbtemp1);
			fbtemp1 *= -1.0;
			fMinv.Multx(fbtemp1, fDb[i]);
		}
		
		if (fOrder > 1)
		{
			/* second derivative */
			int nstr = dSymMatrixT::NumValues(fNumSD);
			for (int ij = 0; ij < nstr; ij++)
			{
				/* resolve indices */
				int i, j;
				dSymMatrixT::ExpandIndex(fNumSD, ij, i, j);
		
				fDDM[ij].Multx(fb, fbtemp1);
				fDM[i].Multx(fDb[j], fbtemp2);
				fDM[j].Multx(fDb[i], fbtemp3);
				fbtemp4.SetToCombination(-1.0, fbtemp1,
				                         -1.0, fbtemp2,
				                         -1.0, fbtemp3);
				fMinv.Multx(fbtemp4, fDDb[ij]);
			}
		}
	}
}

/* set correction function */
void MLSSolverT::SetCorrection(void)
{
	//NOTE: the loops for [nnd] and [nbasis] are reversed in here
	//      since [nnd] > [nbasis]

	const dArray2DT& basis = fBasis->P();

	/* correction function */
	fC = 0.0;
	for (int j = 0; j < fb.Length(); j++)
	{
		double* C = fC.Pointer();
		double* P = basis(j);
		double  b = fb[j];
		for (int n = 0; n < fNumNeighbors; n++)
			*C++ += (*P++)*b;
	}
	
	if (fOrder > 0)
	{
		/* 1st derivative */
		fDC = 0.0;
		for (int i = 0; i < fNumSD; i++)
		{
			const dArray2DT& Dbasis = fBasis->DP(i);
			for (int j = 0; j < fb.Length(); j++)
			{
				double* DC = fDC(i);
				double*  P = basis(j);				
				double* DP = Dbasis(j);				
				double   b = fb[j];
				double  Db = fDb[i][j];
				for (int n = 0; n < fNumNeighbors; n++)
					*DC++ += (*DP++)*b + (*P++)*Db;
			}
		}
		
		if (fOrder > 1)
		{
			/* 2nd derivative */
			fDDC = 0.0;
			for (int ij = 0; ij < dSymMatrixT::NumValues(fNumSD); ij++)
			{
				/* expand index */
				int i, j;
				dSymMatrixT::ExpandIndex(fNumSD, ij, i, j);
				const dArray2DT& Dbasis_i = fBasis->DP(i);
				const dArray2DT& Dbasis_j = fBasis->DP(j);
				const dArray2DT&  DDbasis = fBasis->DDP(ij);
				for (int m = 0; m < fb.Length(); m++)
				{
					double* DDC = fDDC(ij);
					double*   P = basis(m);
					double* DPi = Dbasis_i(m);
					double* DPj = Dbasis_j(m);
					double* DDP = DDbasis(m);
					double    b = fb[m];
					double  Dbi = fDb[i][m];
					double  Dbj = fDb[j][m];
					double  DDb = fDDb[ij][m];
					for (int n = 0; n < fNumNeighbors; n++)
						*DDC++ += (*DDP++)*b
						        + (*DPi++)*Dbj
						        + (*DPj++)*Dbi
						        + (*P++)*DDb;
					
				}
			}
		}
	}
}

/* set shape functions */
void MLSSolverT::SetShapeFunctions(const dArrayT& volume)
{
	/* shape function */
	double* phi = fphi.Pointer();
	double*   C = fC.Pointer();
	double*   w = fw.Pointer();
	double* vol = volume.Pointer();
	for (int i = 0; i < fNumNeighbors; i++)
		*phi++ = (*C++)*(*w++)*(*vol++);

	if (fOrder > 0)
	{
		for (int j = 0; j < fNumSD; j++)
		{
			double* Dphi = fDphi(j);
			double*    C = fC.Pointer();
			double*   DC = fDC(j);
			double*    w = fw.Pointer();
			double*   Dw = fDw(j);
			double*  vol = volume.Pointer();
			for (int i = 0; i < fNumNeighbors; i++)
				*Dphi++ = ((*DC++)*(*w++)
				         + (*C++)*(*Dw++))*(*vol++);
		}
		
		if (fOrder > 1)
		{
			for (int ij = 0; ij < dSymMatrixT::NumValues(fNumSD); ij++)
			{
				/* resolve index */
				int i, j;
				dSymMatrixT::ExpandIndex(fNumSD, ij, i, j);
			
				double* DDphi = fDDphi(ij);
				double*     C = fC.Pointer();
				double*  DC_i = fDC(i);
				double*  DC_j = fDC(j);
				double*   DDC = fDDC(ij);
				double*     w = fw.Pointer();
				double*  Dw_i = fDw(i);
				double*  Dw_j = fDw(j);
				double*   DDw = fDDw(ij);
				double*   vol = volume.Pointer();
				for (int n = 0; n < fNumNeighbors; n++)
					*DDphi++ = ((*DDC++)*(*w++)
					          + (*DC_i++)*(*Dw_j++)
					          + (*DC_j++)*(*Dw_i++)
					          + (*C++)*(*DDw++))*(*vol++);
			}
		}
	}
}

/* replace M with its inverse */
int MLSSolverT::SymmetricInverse3x3(dMatrixT& M)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16;
	
	double m00 = M(0,0);
	double m11 = M(1,1);
	double m22 = M(2,2);
	double m12 = M(1,2);
	double m02 = M(0,2);
	double m01 = M(0,1);
	
	z1 = m01*m01;
	z2 = m01*m02;
	z3 = m02*m02;
	z4 = m00*m11;
	z5 = -m02*m11;
	z6 = -m00*m12;
	z7 = m01*m12;
	z8 = m02*m12;
	z9 = m12*m12;
	z10 = m00*m22;
	z11 = -m01*m22;
	z12 = m11*m22;
	z13 = 2.*m12*z2;
	z14 = m22*z4;
	z1 = -z1;
	z15 = m22*z1;
	z3 = -z3;
	z16 = m11*z3;
	z2 = z2 + z6;
	z5 = z5 + z7;
	z6 = z11 + z8;
	z7 = -z9;
	z8 = m00*z7;
	z3 = z10 + z3;
	z7 = z12 + z7;
	z8 = z13 + z14 + z15 + z16 + z8;

	// check the determinant
	if (fabs(z8) < kSmall)
	{
		cout << "\n MLSSolverT::SymmetricInverse3x3: matrix is singular" << endl;
		return 0;
	}
	
	z1 = z1 + z4;
	z4 = 1./z8;
	z2 = z2*z4;
	z5 = z4*z5;
	z6 = z4*z6;
	z3 = z3*z4;
	z7 = z4*z7;
	z1 = z1*z4;
	
	// return value
	// {{z7, z6, z5},
	//  {z6, z3, z2},
	//  {z5, z2, z1}}

	M(0,0) = z7;
	M(1,1) = z3;
	M(2,2) = z1;
	M(1,2) = M(2,1) = z2;
	M(0,2) = M(2,0) = z5;
	M(0,1) = M(1,0) = z6;
	
	return 1;
}

int MLSSolverT::SymmetricInverse4x4(dMatrixT& M)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39, z40, z41, z42, z43, z44, z45, z46, z47, z48;
	double z49, z50, z51, z52, z53, z54, z55, z56, z57, z58, z59, z60;
	double z61, z62, z63, z64, z65, z66, z67, z68, z69, z70, z71, z72;
	double z73;

	double m00 = M(0,0);
	double m11 = M(1,1);
	double m22 = M(2,2);
	double m33 = M(3,3);
	double m23 = M(2,3);
	double m13 = M(1,3);
	double m12 = M(1,2);
	double m03 = M(0,3);
	double m02 = M(0,2);
	double m01 = M(0,1);

	z1 = m01*m01;
	z2 = m02*m02;
	z3 = m02*m03;
	z4 = m03*m03;
	z5 = 2.*m01*m02*m12;
	z6 = -m01*m03*m12;
	z7 = m12*m12;
	z8 = -m01*m02*m13;
	z9 = 2.*m01*m03*m13;
	z10 = m00*m12*m13;
	z11 = -m02*m12*m13;
	z12 = -m03*m12*m13;
	z13 = m13*m13;
	z14 = m01*m03*m22;
	z15 = m00*m11*m22;
	z16 = -m03*m11*m22;
	z17 = -m00*m13*m22;
	z18 = m01*m13*m22;
	z19 = m03*m13*m22;
	z20 = -m01*m02*m23;
	z21 = -m01*m03*m23;
	z22 = -m00*m11*m23;
	z23 = m02*m11*m23;
	z24 = m03*m11*m23;
	z25 = m00*m12*m23;
	z26 = -m01*m12*m23;
	z27 = -m03*m12*m23;
	z28 = -2.*m01*m03*m12*m23;
	z29 = m00*m13*m23;
	z30 = -m01*m13*m23;
	z31 = -m02*m13*m23;
	z32 = -2.*m01*m02*m13*m23;
	z33 = 2.*m12*m13*m23;
	z34 = m23*m23;
	z35 = m01*m02*m33;
	z36 = m00*m11*m33;
	z37 = -m02*m11*m33;
	z38 = -m00*m12*m33;
	z39 = m01*m12*m33;
	z40 = m02*m12*m33;
	z41 = m00*m22*m33;
	z42 = -m01*m22*m33;
	z43 = m11*m22*m33;
	z44 = 2.*m23*z10;
	z45 = m33*z15;
	z46 = m11*z3;
	z47 = -m12*z3;
	z48 = -m13*z3;
	z49 = -2.*m12*m13*z3;
	z3 = 2.*m23*z3;
	z50 = 2.*m23*z46;
	z51 = m33*z5;
	z52 = m22*z9;
	z53 = -m22*z1;
	z54 = m23*z1;
	z55 = -m33*z1;
	z56 = -m00*z13;
	z57 = m02*z13;
	z58 = -m22*z13;
	z59 = -m11*z2;
	z60 = m13*z2;
	z61 = -m33*z2;
	z2 = z13*z2;
	z13 = -m00*z34;
	z62 = m01*z34;
	z63 = -m11*z34;
	z1 = z1*z34;
	z34 = -m11*z4;
	z64 = m12*z4;
	z65 = -m22*z4;
	z66 = m33*z53;
	z67 = m22*z56;
	z68 = m33*z59;
	z69 = m11*z13;
	z70 = m22*z34;
	z71 = -m00*z7;
	z72 = m03*z7;
	z73 = -m33*z7;
	z4 = z4*z7;
	z7 = m33*z71;
	z12 = z12 + z24 + z30 + z37 + z39 + z57;
	z14 = z14 + z17 + z20 + z25 + z47 + z60;
	z17 = z19 + z27 + z31 + z40 + z42 + z62;
	z19 = z21 + z29 + z35 + z38 + z48 + z64;
	z3 = z13 + z3 + z41 + z61 + z65;
	z5 = z15 + z5 + z53 + z59 + z71;
	z11 = z11 + z16 + z18 + z23 + z26 + z72;
	z13 = z33 + z43 + z58 + z63 + z73;
	z6 = z10 + z22 + z46 + z54 + z6 + z8;
	z1 = z1 + z2 + z28 + z32 + z44 + z45 + z49 + z50 + z51 + z52;
	z2 = z34 + z36 + z55 + z56 + z9;
	z1 = z1 + z4 + z66 + z67 + z68 + z69 + z7 + z70;

	// check the determinant
	if (fabs(z1) < kSmall)
	{
		cout << "\n MLSSolverT::SymmetricInverse4x4: matrix is singular" << endl;
		return 0;
	}

	z1 = 1./z1;
	z4 = z1*z12;
	z7 = z1*z14;
	z8 = z1*z17;
	z9 = z1*z19;
	z3 = z1*z3;
	z5 = z1*z5;
	z10 = z1*z11;
	z11 = z1*z13;
	z6 = z1*z6;
	z1 = z1*z2;

	// return value
	// {{z11, z8, z4, z10},
	//  {z8, z3, z9, z7},
	//  {z4, z9, z1, z6},
	//  {z10, z7, z6, z5}}

	M(0,0) = z11;
	M(1,1) = z3;
	M(2,2) = z1;
	M(3,3) = z5;
	M(2,3) = M(3,2) = z6;
	M(1,3) = M(3,1) = z7;
	M(1,2) = M(2,1) = z9;
	M(0,3) = M(3,0) = z10;
	M(0,2) = M(2,0) = z4;
	M(0,1) = M(1,0) = z8;

	return 1;
}
