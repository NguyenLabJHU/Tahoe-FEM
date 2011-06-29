/* $Id: OrthoMLSSolverT.cpp,v 1.8 2004-06-26 06:11:09 paklein Exp $ */
/* created: paklein (07/03/1998) */
#include "OrthoMLSSolverT.h"

using namespace Tahoe;

/* constants */
const double sqrtPi = sqrt(acos(-1.0));

/* constructor */
OrthoMLSSolverT::OrthoMLSSolverT(int nsd, int complete):
	fNumSD(nsd),
	fComplete(complete), //1 => linear, which has 2 terms
	fNumNeighbors(0),
	fDqJ(fNumSD),

	/* variable memory managers */
	fArrayGroup(0, true),
	fArray2DGroup1(0, 0),
	fArray2DGroup2(0, 0),
	fLocCoords_man(0, fLocCoords, fNumSD)
{
	/* error checking */
	if (fNumSD < 1 && fNumSD > 3) throw ExceptionT::kBadInputValue;
}
	
/* destructor */
OrthoMLSSolverT::~OrthoMLSSolverT(void) { }

/* class dependent initializations */
void OrthoMLSSolverT::Initialize(void)
{
	/* number of monomials */
	int m = NumberOfMonomials(fComplete);

	/* b (2.11 b) and derivatives */
	fb.Dimension(m);
	fDb.Dimension(fNumSD, m);

	/* monomials and derivatives */
	fp.Dimension(m);
	fDp.Dimension(fNumSD, m);

	/* (orthogonal) basis functions (2.4) and derivatives */
	fq.Dimension(m);
	fDq.Dimension(fNumSD, m);

	/* rows of alpha and derivatives */
	fa.Dimension(m-1);
	fDa.Dimension(fNumSD, m-1);
	
	/* register variable length arrays */
	fArrayGroup.Register(fw);
	fArrayGroup.Register(fphi);
	
	fArray2DGroup1.Register(fpJ);
	fArray2DGroup1.Register(fqJ);
	for (int i = 0; i < fNumSD; i++)
		fArray2DGroup1.Register(fDqJ[i]);
		
	fArray2DGroup2.Register(fDw);		
	fArray2DGroup2.Register(fDphi);		
}

/* set MLS at coords given sampling points */
int OrthoMLSSolverT::SetField(const dArray2DT& nodalcoords,
	const nArrayT<double>& dmax, const dArrayT& samplept)
{
#if __option(extended_errorcheck)
	if (dmax.Length() != nodalcoords.MajorDim()) throw ExceptionT::kSizeMismatch;
	if (samplept.Length() != fNumSD) throw ExceptionT::kSizeMismatch;
#endif

	/* set size of current working set */
	Dimension(nodalcoords.MajorDim());

	/* compute local coordinates */
	fLocCoords = nodalcoords;
	dArrayT x_sh;
	for (int k = 0; k < fNumNeighbors; k++)
	{
		/* wrt local origin */
		fLocCoords.RowAlias(k, x_sh);
		x_sh -= samplept;
	}

	/* compute nodal point weights and derivatives */
	int numactive = SetWeight(fLocCoords, dmax);
	if (numactive < NumberOfMonomials(fComplete))
	{
		//TEMP - will handle externally later
		cout << "\n OrthoMLSSolverT::SetField: not enough nodes for fit: ";
		cout << numactive << '/' << NumberOfMonomials(fComplete) << '\n';
		return 0;
	}
	
	/* initialize b */
	fb  = 0.0;
	fDb = 0.0;
	
	/* set {q,Dq}, {b,Dq}, and {q,Dq}_nodal */
	dArrayT origin(fNumSD);
	origin = 0.0;	
	SetBasis(fLocCoords,origin);

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
			(*phi++) += cj*(*qJ++)*(*w++);
	}

	/* derivatives of phi */
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
			
				(*Dphi++) += (*Dw)*CjI + (*w)*dCjI;
				
				w++; qJ++; Dw++; DqJ++;
			}
			
			b++; q++; Db++; Dq++;
		}
	}

	/* OK */
	return 1;
}

/* orthogonality checking (dimensions mat) */
void OrthoMLSSolverT::ComputeOrtho(dMatrixT& mat) const
{
	/* dimension output matrix */
	int m = fq.Length();
	mat.Dimension(m);

	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
		{
			double sum = 0.0;
			
			const double* w = fw.Pointer();
			const double* qi = fqJ(i);
			const double* qj = fqJ(j);
			for (int J = 0; J < fNumNeighbors; J++)		
				sum += (*w++)*(*qi++)*(*qj++);
		
			mat(i,j) = sum;
		}
}

void OrthoMLSSolverT::ComputeDOrtho(int deriv, dMatrixT& mat) const
{
	/* dimension output matrix */
	int m = fq.Length();
	mat.Dimension(m);

	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
		{
			double sum = 0.0;
			
			const double* w = fw.Pointer();
			const double* qi = fqJ(i);
			const double* qj = fqJ(j);

			const double* Dw  = fDw(deriv);
			const double* Dqi = (fDqJ[deriv])(i);
			const double* Dqj = (fDqJ[deriv])(j);
			for (int J = 0; J < fNumNeighbors; J++)
			{
				sum += (*Dw++)*(*qi)*(*qj) +
				       (*w)*(*Dqi++)*(*qj) +
				       (*w)*(*qi)*(*Dqj++);
				
				w++; qi++; qj++;
			}
		
			mat(i,j) = sum;
		}
}

/* set weight functions and derivatives */
int OrthoMLSSolverT::SetWeight(const dArray2DT& localcoords,
	const ArrayT<double>& dmax)
{
	dArrayT dx, Dw(fNumSD); //make into work space variables?
	int count = 0;
	for (int i = 0; i < fNumNeighbors; i++)
	{
		/* fetch local coords */
		localcoords.RowAlias(i,dx);
		
		/* weight function */
		double& w = fw[i];

		double dm = dmax[i];
		double di = dx.Magnitude();
	
		/* out of influence range (or inactive) */
		if (di > dm)
		{
			 w = 0.0;
			
			 fDw.SetColumn(i,0.0);
		}
		/* Gaussian weight function and derivatives */
		else
		{
			/* influence scaling */
			double a = 0.4;
			double q = di/(a*dm);
			
			/* weight */
	w=exp(-q*q)/(sqrtPi*a*dm);

			/* derivatives */
			Dw.SetToScaled(2.0*w/pow(a*dm,2.0),dx);

			/* store */
			fDw.SetColumn(i,Dw);
			
			count++;
		}
	}
	
	return count;
}

/***********************************************************************
* Protected
***********************************************************************/

/* set basis functions and derivatives */
void OrthoMLSSolverT::SetBasis(const dArray2DT& localcoords,
	const dArrayT& samplept)
{
	/* get p at all neighbors */
	dArrayT dx;
	for (int i = 0; i < fNumNeighbors; i++)
	{
		/* fetch local coords */
		localcoords.RowAlias(i,dx);
	
		/* mononials and derivatives */
		SetMonomials(dx,fp,fDp);
	
		/* store nodal p */
		fpJ.SetColumn(i,fp);
		
		/* write 1st monomial of nodal p/Dp to nodal q/Dq */
		fqJ(0,i) = fp[0];
		for (int k = 0; k < fNumSD; k++) //not all derivatives of 1st
		{                                //term are just zero!
			dArray2DT& DqJ = fDqJ[k];
			DqJ(0,i) = fDp(k,0);
		}
	}
	
	/* initialize b */
	Setb(0);
	
	/* initialize q */
	SetMonomials(samplept, fq, fDq);
	
	/* Schmidt orthogonalization */
	for (int j = 1; j < fq.Length(); j++)
	{
		/* set row of alpha (up to col j-1) */
		Setalpha(j,j-1);
		
		/* int pt basis */
		double sum = 0.0;
		double* a  = fa.Pointer();
		double* q  = fq.Pointer();
		for (int k = 0; k < j; k++)
			sum += (*a++)*(*q++);			

		fq[j] -= sum;	

		/* int pt basis derivatives */
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
		
			/* nodal basis derivatives */
			for (int i = 0; i < fNumSD; i++)
			{
				dArray2DT& DqJ = fDqJ[i];
			
				double sum = 0.0;
		
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
	
				DqJ(j,J) = 0.0 - sum;
			}
		}
	
		/* next term in b */
		Setb(j);
	}
}

/* set b and derivatives */
void OrthoMLSSolverT::Setb(int pterm)
{
	/* set b */
	double& b = fb[pterm];
	double* q = fqJ(pterm);
	double* w = fw.Pointer();
	for (int J = 0; J < fNumNeighbors; J++)
	{
		b += (*w++)*(*q)*(*q); // (2.11 b)
		q++;
	}
	
	/* derivatives of b */
	for (int k = 0; k < fNumSD; k++)
	{
		double& Db = fDb(k, pterm);
	
		double* q = fqJ(pterm);
		double* w = fw.Pointer();

		double* Dw = fDw(k);
		double* Dq = (fDqJ[k])(pterm);

		for (int J = 0; J < fNumNeighbors; J++)
		{
			Db += (*Dw++)*(*q)*(*q) + 2.0*(*w++)*(*q)*(*Dq++); // (2.16 b)
			q++;
		}
	}
}

/* set Schmidt orthogonalization */
void OrthoMLSSolverT::Setalpha(int row, int maxcol)
{
	for (int i = 0; i <= maxcol; i++)
	{
		/* alpha (2.5) */
		double a = 0.0;

		double* w = fw.Pointer();
		double* p = fpJ(row);
		double* q = fqJ(i);

		for (int J = 0; J < fNumNeighbors; J++)
			a += (*p++)*(*w++)*(*q++);
			
		fa[i] = a/fb[i];	

		/* derivatives of alpha - from (2.5) */		
		for (int k = 0; k < fNumSD; k++)
		{
			double Da = 0.0;

			double* w  = fw.Pointer();
			double* p  = fpJ(row);
			double* q  = fqJ(i);
			double* Dw = fDw(k);
			double* Dq = (fDqJ[k])(i);

			for (int J = 0; J < fNumNeighbors; J++)
				Da += (*p++)*((*Dw++)*(*q++) + (*w++)*(*Dq++));

			fDa(k,i) = (Da - fa[i]*fDb(k,i))/fb[i];
		}
	}
}

/* configure solver for current number of neighbors */
void OrthoMLSSolverT::Dimension(int numneighbors)
{
	fNumNeighbors = numneighbors;

	/* set variable memory */
	fArrayGroup.Dimension(fNumNeighbors, false);
	fArray2DGroup1.Dimension(NumberOfMonomials(fComplete), fNumNeighbors);
	fArray2DGroup2.Dimension(fNumSD, fNumNeighbors);
	fLocCoords_man.SetMajorDimension(fNumNeighbors, false);
}
