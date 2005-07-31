/* $Id: D2OrthoMLSSolverT.cpp,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (10/17/1999)                                          */

#include "D2OrthoMLSSolverT.h"

#include "ExceptionCodes.h"
#include "dSymMatrixT.h"

/* constants */
const double sqrtPi = sqrt(acos(-1.0));

/* constructor */
D2OrthoMLSSolverT::D2OrthoMLSSolverT(int nsd, int complete):
	OrthoMLSSolverT(nsd, complete),
	fDDqJ(dSymMatrixT::NumValues(fNumSD)),

	/* variable memory managers */
	fArray2DGroup3(0, 0),
	
	/* work space */
	fNSDsym(fNumSD)
{

}

/* class dependent initializations */
void D2OrthoMLSSolverT::Initialize(void)
{
	/* inherited */
	OrthoMLSSolverT::Initialize();

	/* dimensions */
	int   m = NumberOfMonomials(fComplete);
	int dim = dSymMatrixT::NumValues(fNumSD);

	/* b (2.11 b) derivatives */
	fDDb.Allocate(dim, m);

	/* monomial derivatives */
	fDDp.Allocate(dim, m);

	/* (orthogonal) basis functions derivatives */
	fDDq.Allocate(dim, m);

	/* rows of alpha derivatives */
	fDDa.Allocate(dim, m-1);
	
	/* register variable length arrays */
	for (int i = 0; i < dim; i++)
		fArray2DGroup1.Register(fDDqJ[i]);
		
	fArray2DGroup3.Register(fDDw);		
	fArray2DGroup3.Register(fDDphi);		
}

/* set MLS at coords given sampling points */
int D2OrthoMLSSolverT::SetField(const dArray2DT& nodalcoords,
	const dArrayT& dmax, const dArrayT& samplept)
{
#if __option(extended_errorcheck)
	if (dmax.Length() != nodalcoords.MajorDim()) throw eSizeMismatch;
	if (samplept.Length() != fNumSD) throw eSizeMismatch;
#endif

	/* set size of current working set */
	Dimension(nodalcoords.MajorDim());

	/* compute local coordinates */
	fLocCoords = nodalcoords;
	dArrayT x_sh;
	for (int k = 0; k < fNumNeighbors; k++)
	{
		/* wrt local origin */
		fLocCoords.RowAlias(k,x_sh);
		x_sh -= samplept;
	}

	/* compute nodal point weights and derivatives */
	int numactive = SetWeight(fLocCoords, dmax);
	if (numactive < NumberOfMonomials(fComplete))
	{
		//TEMP - will handle externally later
		cout << "\n D2OrthoMLSSolverT::SetField: not enough nodes for fit: ";
		cout << numactive << '/' << NumberOfMonomials(fComplete) << endl;
		throw eGeneralFail;
		
		return 0;
	}
	
	/* initialize b */
	fb = 0.0;
	fDb = 0.0;
	fDDb = 0.0;
	
	/* set {q,Dq}, {b,Dq}, and {q,Dq}_nodal */
	dArrayT origin(fNumSD);
	origin = 0.0;	
	SetBasis(fLocCoords, origin);

	/* phi */
	fphi = 0.0;
	double* b  = fb.Pointer();
	double* q  = fq.Pointer();
	for (int j = 0; j < fq.Length(); j++)
	{
		double cj = (*q++)/(*b++);
	
		double* phi = fphi.Pointer();
		double* w   = fw.Pointer();
		double* qJ  = fqJ(j);
		for (int J = 0; J < fNumNeighbors; J++)
			*phi++ += cj*(*qJ++)*(*w++);
			//testing
			//*phi++ += cj*(*qJ++);
	}

	/* 1st derivatives of phi */
	fDphi = 0.0;
	for (int i = 0; i < fNumSD; i++)
	{
		double* b = fb.Pointer();
		double* q = fq.Pointer();
		double* Db = fDb(i);
		double* Dq = fDq(i);
		
		for (int j = 0; j < fq.Length(); j++)
		{
			double* w    = fw.Pointer();
			double* qJ   = fqJ(j);
			double* DqJ  = (fDqJ[i])(j);
			double* Dw   = fDw(i); 		
			double* Dphi = fDphi(i);

			for (int J = 0; J < fNumNeighbors; J++)
			{
				double CjI  = (*q)*(*qJ)/(*b);
				double dCjI = ((*Dq)*(*qJ) + (*q)*(*DqJ) - (*Db)*CjI)/(*b);
			
				*Dphi++ += (*Dw)*CjI + (*w)*dCjI;
				//testing
				//*Dphi++ += dCjI;
				
				w++; qJ++;
				Dw++; DqJ++;
			}
			
			b++; q++; Db++; Dq++;
		}
	}

	/* 2nd derivatives of phi */
	fDDphi = 0.0;
	int stressdim = dSymMatrixT::NumValues(fNumSD);
	for (int l = 0; l < stressdim; l++)
	{
		/* expand stress index */
		int dex_1, dex_2;
		dSymMatrixT::ExpandIndex(fNumSD, l, dex_1, dex_2);

		double* b = fb.Pointer();
		double* q = fq.Pointer();

		double* Db_1 = fDb(dex_1);
		double* Db_2 = fDb(dex_2);
		double* Dq_1 = fDq(dex_1);
		double* Dq_2 = fDq(dex_2);

		double* DDb = fDDb(l);
		double* DDq = fDDq(l);
		
		for (int j = 0; j < fq.Length(); j++)
		{
			double* w  = fw.Pointer();
			double* qJ = fqJ(j);

			double* DqJ_1  = (fDqJ[dex_1])(j);
			double* DqJ_2  = (fDqJ[dex_2])(j);
			double* Dw_1   = fDw(dex_1); 		
			double* Dw_2   = fDw(dex_2); 		

			double* DDqJ  = (fDDqJ[l])(j);
			double* DDw   = fDDw(l); 		

			double* DDphi = fDDphi(l);
			for (int J = 0; J < fNumNeighbors; J++)
			{
				double CjI = (*q)*(*qJ)/(*b);

				double dCjI_1 = ((*Dq_1)*(*qJ) + (*q)*(*DqJ_1) - (*Db_1)*CjI)/(*b);
				double dCjI_2 = ((*Dq_2)*(*qJ) + (*q)*(*DqJ_2) - (*Db_2)*CjI)/(*b);

				double ddCjI = ((*DDq)*(*qJ) + (*q)*(*DDqJ) +
				                (*Dq_1)*(*DqJ_2) + (*Dq_2)*(*DqJ_1) -
				               ((*DDb) + (*Db_1)*dCjI_2 + (*Db_2)*dCjI_1))/(*b);
			
				*DDphi++ += (*DDw)*CjI +
				            (*Dw_1)*dCjI_2 + (*Dw_2)*dCjI_1 +
			               (*w)*ddCjI;

				//testing
				//*DDphi++ += ddCjI;

				w++; qJ++;
				Dw_1++; Dw_2++; DqJ_1++; DqJ_2++;
				DDqJ++; DDw++;
			}
			
			b++; q++;
			Db_1++; Db_2++; Dq_1++; Dq_2++;
			DDb++; DDq++;
		}
	}

	/* OK */
	return 1;
}

/***********************************************************************
* Protected
***********************************************************************/

/* configure solver for current number of neighbors */
void D2OrthoMLSSolverT::Dimension(int numneighbors)
{
	/* inherited */
	OrthoMLSSolverT::Dimension(numneighbors);

	/* set variable memory */
	fArray2DGroup3.Dimension(dSymMatrixT::NumValues(fNumSD),
		fNumNeighbors);
}

/***********************************************************************
* Private
***********************************************************************/

/* set weight functions and derivatives */
int D2OrthoMLSSolverT::SetWeight(const dArray2DT& localcoords,
	const ArrayT<double>& dmax)
{
	dArrayT dx;
	dArrayT Dw(fNumSD); //make into work space variables
	int count = 0;
	for (int i = 0; i < fNumNeighbors; i++)
	{
		/* fetch local coords */
		localcoords.RowAlias(i, dx);
		
		/* weight function */
		double& w = fw[i];

		double dm = dmax[i];
		double di = dx.Magnitude();
	
		/* out of influence range (or inactive) */
		if (di > 3.0*dm)
		{
			 w = 0.0;
			 fDw.SetColumn(i, 0.0);
			 fDDw.SetColumn(i, 0.0);
		}
		/* Gaussian weight function and derivatives */
		else
		{
			double    a = 0.4;
			double  adm = a*dm;
			double adm2 = adm*adm;
			double    q = di/adm;
			
			/* weight */
	w = exp(-q*q)/(sqrtPi*adm);

			/* 1st derivative */
			Dw.SetToScaled(2.0*w/adm2, dx);

			/* store */
			fDw.SetColumn(i, Dw);

			/* 2nd derivative */
			fNSDsym.Outer(dx);
			fNSDsym *= 4.0*w/(adm2*adm2);
			fNSDsym.PlusIdentity(-2.0*w/adm2);

			/* store */
			fDDw.SetColumn(i, fNSDsym);

			count++;
		}
	}
	
	return count;
}

/* set basis functions and derivatives */
void D2OrthoMLSSolverT::SetBasis(const dArray2DT& localcoords,
	const dArrayT& samplept)
{
	/* get p at all neighbors */
	dArrayT dx;
	for (int i = 0; i < fNumNeighbors; i++)
	{
		/* fetch local coords */
		localcoords.RowAlias(i, dx);
	
		/* mononials and derivatives */
		_SetMonomials(dx, fp, fDp, fDDp);
	
		/* store nodal p */
		fpJ.SetColumn(i,fp);
		
		/* write 1st monomial of nodal p/Dp/DDp to nodal q/Dq/DDq */
		fqJ(0,i) = fp[0];

		for (int k = 0; k < fNumSD; k++)
		{
			dArray2DT& DqJ = fDqJ[k];
			DqJ(0,i) = fDp(k,0);
		}

		int stressdim = dSymMatrixT::NumValues(fNumSD);
		for (int j = 0; j < stressdim; j++)
		{
			dArray2DT& DDqJ = fDDqJ[j];
			DDqJ(0,i) = fDDp(j,0);
		}
	}
	
	/* initialize b */
	Setb(0);
	
	/* initialize q */
	_SetMonomials(samplept, fq, fDq, fDDq);
	
	/* Schmidt orthogonalization */
	int stressdim = dSymMatrixT::NumValues(fNumSD);
	for (int j = 1; j < fq.Length(); j++)
	{
		/* set row of alpha (up to col j-1) */
		Setalpha(j, j-1);
		
		/* int pt basis */
		double sum = 0.0;
		double* a  = fa.Pointer();
		double* q  = fq.Pointer();
		for (int k = 0; k < j; k++)
			sum += (*a++)*(*q++);			

		fq[j] -= sum;	

		/* int pt basis 1st derivatives */
		for (int i = 0; i < fNumSD; i++)
		{
			double sum = 0.0;
			double* a  = fa.Pointer();
			double* q  = fq.Pointer();
			double* Da = fDa(i);
			double* Dq = fDq(i);
			for (int k = 0; k < j; k++)
				sum += (*Da++)*(*q++) + (*a++)*(*Dq++);			

			fDq(i,j) -= sum;
		}

		/* int pt basis 2nd derivatives */
		for (int l = 0; l < stressdim; l++)
		{
			double sum = 0.0;
			double*  a = fa.Pointer();
			double*  q = fq.Pointer();

			/* expand stress index */
			int dex_1, dex_2;
			dSymMatrixT::ExpandIndex(fNumSD, l, dex_1, dex_2);
			double* Da_1 = fDa(dex_1);
			double* Da_2 = fDa(dex_2);
			double* Dq_1 = fDq(dex_1);
			double* Dq_2 = fDq(dex_2);

			double* DDa = fDDa(l);
			double* DDq = fDDq(l);
			for (int k = 0; k < j; k++)
				sum += (*DDa++)*(*q++) +
				       (*Da_1++)*(*Dq_2++) +
				       (*Da_2++)*(*Dq_1++) +
				       (*a++)*(*DDq++);			

			fDDq(l,j) -= sum;
			//TEMP - seems like this should always be zero for
			//       linear basis functions - hard wire for 1D tests
			//if (fNumSD == 1) fDDq(l,j) = 0.0;
		}
		
		/* nodal basis functions */
		for (int J = 0; J < fNumNeighbors; J++)
		{
			/* nodal basis */
			double sum = 0.0;
			double*  a = fa.Pointer();
			double*  q = fqJ(0) + J;
			for (int k = 0; k < j; k++)
			{
				sum += (*a++)*(*q);
				 q  += fNumNeighbors;	
			}

			fqJ(j,J) = fpJ(j,J) - sum;
		
			/* nodal basis 1st derivatives */
			for (int i = 0; i < fNumSD; i++)
			{
				double sum = 0.0;
		
				dArray2DT& DqJ = fDqJ[i];
				double*  a = fa.Pointer();
				double*  q = fqJ(0) + J;
				double* Da = fDa(i);
				double* Dq = DqJ(0) + J;
				for (int k = 0; k < j; k++)
				{
					sum += (*Da++)*(*q) + (*a++)*(*Dq);
					q   += fNumNeighbors;
					Dq  += fNumNeighbors; 	
				}
	
				DqJ(j,J) = 0.0 - sum; // why is this zero here???
			}

			/* nodal basis 2nd derivatives */
			for (int l = 0; l < stressdim; l++)
			{			
				double sum = 0.0;
		
				double* a = fa.Pointer();
				double* q = fqJ(0) + J;

				/* expand stress index */
				int dex_1, dex_2;
				dSymMatrixT::ExpandIndex(fNumSD, l, dex_1, dex_2);
				dArray2DT& DqJ_1 = fDqJ[dex_1];
				dArray2DT& DqJ_2 = fDqJ[dex_2];
				double* Da_1 = fDa(dex_1);
				double* Da_2 = fDa(dex_2);
				double* Dq_1 = DqJ_1(0) + J;
				double* Dq_2 = DqJ_2(0) + J;

				dArray2DT& DDqJ = fDDqJ[l];
				double* DDa = fDDa(l);
				double* DDq = DDqJ(0) + J;
				for (int k = 0; k < j; k++)
				{
					sum += (*DDa++)*(*q) +
					       (*Da_1++)*(*Dq_2) +
					       (*Da_2++)*(*Dq_1) +
					       (*a++)*(*DDq);
					
					q    += fNumNeighbors;
					Dq_1 += fNumNeighbors;
					Dq_2 += fNumNeighbors;
					DDq  += fNumNeighbors; 	
				}
	
				DDqJ(j,J) = 0.0 - sum;
				//TEMP - seems like this should always be zero for
				//       linear basis functions - hard wire for 1D tests
				//if (fNumSD == 1) DDqJ(j,J) = 0.0;
			}
		}
	
		/* next term in b */
		Setb(j);
	}
}

/* set b and derivatives */
void D2OrthoMLSSolverT::Setb(int pterm)
{
	/* inherited */
	OrthoMLSSolverT::Setb(pterm);

	/* derivatives of b */
	int stressdim = dSymMatrixT::NumValues(fNumSD);
	for (int k = 0; k < stressdim; k++)
	{
		double& DDb = fDDb(k, pterm);

		double* q = fqJ(pterm);
		double* w = fw.Pointer();

		/* expand stress index */
		int dex_1, dex_2;
		dSymMatrixT::ExpandIndex(fNumSD, k, dex_1, dex_2);

		double* Dw_1 = fDw(dex_1);
		double* Dw_2 = fDw(dex_2);
		double* Dq_1 = (fDqJ[dex_1])(pterm);
		double* Dq_2 = (fDqJ[dex_2])(pterm);

		double* DDw = fDDw(k);
		double* DDq = (fDDqJ[k])(pterm);
		for (int J = 0; J < fNumNeighbors; J++)
		{
			DDb += (*DDw)*(*q)*(*q) +
			   2.0*((*q)*((*Dw_1)*(*Dq_2) + (*Dw_2)*(*Dq_1)) +
			        (*w)*((*Dq_1)*(*Dq_2) + (*q)*(*DDq)));
			
			w++; q++;
			Dw_1++; Dw_2++; Dq_1++; Dq_2++;
			DDw++; DDq++;
		}
	}
}

/* set Schmidt orthogonalization */
void D2OrthoMLSSolverT::Setalpha(int row, int maxcol)
{
	/* inherited */
	OrthoMLSSolverT::Setalpha(row, maxcol);

	int stressdim = dSymMatrixT::NumValues(fNumSD);
	for (int i = 0; i <= maxcol; i++)
	{
		/* derivatives of alpha - from (2.5) */		
		for (int k = 0; k < stressdim; k++)
		{
			double DDa = 0.0;

			double* w  = fw.Pointer();
			double* p  = fpJ(row);
			double* q  = fqJ(i);

			/* expand stress index */
			int dex_1, dex_2;
			dSymMatrixT::ExpandIndex(fNumSD, k, dex_1, dex_2);

			double* Dw_1 = fDw(dex_1);
			double* Dw_2 = fDw(dex_2);
			double* Dq_1 = (fDqJ[dex_1])(i);
			double* Dq_2 = (fDqJ[dex_2])(i);

			double* DDw = fDDw(k);
			double* DDq = (fDDqJ[k])(i);
			for (int J = 0; J < fNumNeighbors; J++)
			{
				DDa += (*p)*((*DDw)*(*q) +
				             (*w)*(*DDq) +
				             (*Dw_1)*(*Dq_2) +
				             (*Dw_2)*(*Dq_1));
				
				w++; p++; q++;
				Dw_1++; Dw_2++; Dq_1++; Dq_2++;
				DDw++; DDq++;
			}

			fDDa(k,i) = (DDa - fa[i]*fDDb(k,i) -
			             fDa(dex_1, i)*fDb(dex_2, i) -
			             fDa(dex_2, i)*fDb(dex_1, i))/fb[i];
		}
	}
}
