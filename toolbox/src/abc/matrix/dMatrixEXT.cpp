/* $Id: dMatrixEXT.cpp,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (03/06/1998)                                          */

#include "dMatrixEXT.h"
#include "Constants.h"
#include "LAdMatrixT.h"
#include "dArray2DT.h"

/* Numerical Recipies macros */
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* constructor */
dMatrixEXT::dMatrixEXT(void): v1(NULL), v2(NULL) { }

dMatrixEXT::dMatrixEXT(int squaredim):
	dMatrixT(squaredim),
	fworkspace(2*squaredim)
{
	/* set pointers */
	v1 = fworkspace.Pointer();
	v2 = v1 + squaredim;
}

dMatrixEXT::dMatrixEXT(int squaredim, double* p):
	dMatrixT(squaredim,squaredim,p)
{
	fworkspace.Allocate(2*squaredim);

	/* set pointers */
	v1 = fworkspace.Pointer();
	v2 = v1 + squaredim;
}

/* post constructor (re-)dimensioning */
void dMatrixEXT::Allocate(int squaredim)
{
	/* inherited */
	dMatrixT::Allocate(squaredim);
	
	/* temp space */
	fworkspace.Allocate(2*squaredim);

	/* set pointers */
	v1 = fworkspace.Pointer();
	v2 = v1 + squaredim;
}

/* diagonalize (using symmetric QR algorithm) */
int dMatrixEXT::Diagonalize(dArrayT& eigs)
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (fRows != fCols) throw eGeneralFail;
	if (fRows != eigs.Length()) throw eSizeMismatch;
#endif

	/* make tri-diagonal */
	TriDiagonalForm();

	int    length = eigs.Length();
	double* diags = eigs.Pointer() - 1;
	double* subs  = Pointer() - 1;
	
	/* copy in */
	diags[1] = (*this)(0,0);
	for (int i = 2; i <= length; i++)
	{
		diags[i] = (*this)(i-1,i-1);
		subs[i] = (*this)(i-1,i-2);
	}
	
	/* Numerical Recipies function */
	int count = tqli(diags, subs, eigs.Length());

	/* clear all values and sort/fill with eigenvalues */
	(*this) = 0.0;
	double* peig  = eigs.Pointer();
	double* pdiag = Pointer();
	for (int j = 0; j < fRows; j++)
	{
		*pdiag = *peig++;
		pdiag += (fRows+1);
	}
	
	return count;
}

/* return the eigenvector corresponding to the
* approximate eigenvalue that is passed in.
*/
void dMatrixEXT::Eigenvector(double& eig_guess, dArrayT& eigenvector) const
{
#if __option(extended_errorcheck)
	if (fRows != fCols) throw eGeneralFail;
	if (fRows != eigenvector.Length()) throw eSizeMismatch;
#endif
	
	/* Rayleigh quotient iteration (8.2.3) */

	/* perturb eigenvalue guess (if "exact" is passed in) */
	eig_guess = float(eig_guess);

	/* initial guess */
	//eigenvector = 1.0/sqrt(eigenvector.Length()); //not a godd guess
	eigenvector = 0.0;                            //for element stiffness
	eigenvector[0] = 1.0;

	/* work space */
	dArrayT temp(fRows,v1);
	
	/* error */
	Multx(eigenvector,temp);
	temp.AddScaled(-eig_guess,eigenvector);
	double error0 = temp.Magnitude();
	double error  = error0;
	
	/* iteration */
	LAdMatrixT solver(fRows);
	dMatrixT shsolver;
	shsolver.Alias(solver);
	double tol = 1.0e-7;
	while (error > tol && (error/error0) > tol)
	{
		shsolver = (*this);
		shsolver.PlusIdentity(-eig_guess);

		/* next eigenvector */
		solver.LinearSolve(eigenvector);
		eigenvector.UnitVector();

		/* next eigenvalue (Rayleigh quotient) */
		//if (--hold < 0) eig_guess = MultmBn(eigenvector,eigenvector);
		eig_guess = MultmBn(eigenvector,eigenvector);
		
		/* error */
		Multx(eigenvector,temp);
		temp.AddScaled(-eig_guess,eigenvector);
		error = temp.Magnitude();
	}
}

/* compute tridiagonal decomposition */
void dMatrixEXT::TriDiagonalForm(void)
{
	/* aliases */
	double* v = v1;
	double* p = v2;

	/* algorithm (8.3.1) */
	for (int k = 0; k < fRows - 2; k++)
	{
		int length = fRows - (k + 1);
	
		double beta;
		double* pA = (*this)(k) + (k+1);
		HouseholderVector(pA, v, beta, length);

		double* pAi = pA + fRows;
		for (int i = 0; i < length; i++)
		{
			p[i] = beta*Dot(pAi,v,length);
			pAi += fRows;
		}

		double factor = beta*Dot(p,v,length)/2.0;
		double* pv = v;
		double* pw = p;
		for (int j = 0; j < length; j++)
			(*pw++) -= factor*(*pv++);
			
		*pA = *(pA + fRows - 1) = sqrt(Dot(pA,pA,length));
		
		double* pAn = pA + fRows;
		for (int n = 0; n < length; n++)
		{
			double* wn = p+n;
			double* vn = v+n;
		
			/* diagonal value */
			*pAn -= 2.0*(*wn)*(*vn);
		
			/* row and column together */
			double* pcA = pAn + 1;
			double* prA = pAn + fRows;
			double* wm = wn + 1;
			double* vm = vn + 1;
			for (int m = n+1; m < length; m++)
			{
				double factor = (*wn)*(*vm++) + (*wm++)*(*vn);
				*pcA -= factor;
				*prA -= factor;
			
				pcA++;
				prA += fRows;
			}
			
			/* next diagonal */
			pAn += (fRows+1);
		}
	}
}

/* returns the Householder vector (5.1.2-3), where the
* Householder reflection is given by:
*
*           P = I - 2/(v.v) v (x) v
*
* where beta is then defined as 2/(v.v)
*
* such that P.x = sqrt(x.x) e_1
*/
void dMatrixEXT::HouseholderVector(const double* x, double* v, double& beta,
	int length) const
{
	/* algorithm (5.1.1) */

	double sigma = Dot(x+1,x+1,length-1);

	if (fabs(sigma) < kSmall)
		beta = 0.0;
	else
	{
		double mu = sqrt(x[0]*x[0] + sigma);
		double v0 = (x[0] > 0.0) ? -sigma/(x[0] + mu) : x[0] - mu;
		beta   = 2*v0*v0/(sigma + v0*v0);
		
		v[0] = 1.0;
		double*       pv = v+1;
		const double* px = x+1;
		for (int i = 1; i < length; i++)
			(*pv++) = (*px++)/v0;
	}
}

double dMatrixEXT::pythag(double a, double b)
{
	double absa=fabs(a);
	double absb=fabs(b);

	if (absa > absb)
		return (absa*sqrt(1.0+SQR(absb/absa)));
	else
		return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

int dMatrixEXT::tqli(double d[], double e[], int n)
{
	int m,l,iter,i;
	double s,r,p,g,f,dd,c,b;

	for (i=2;i<=n;i++) e[i-1]=e[i];
	e[n]=0.0;
	for (l=1;l<=n;l++) {
		iter=0;
		do {
			for (m=l;m<=n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if ((double)(fabs(e[m])+dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == 3*n)
				{
					cout << "\n dMatrixEXT::tqli: no convergence after " << 3*n;
					cout << " iterations." << endl;
					throw eGeneralFail;
				}
					
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
	
	return iter;
}

