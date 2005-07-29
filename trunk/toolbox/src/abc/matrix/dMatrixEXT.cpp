/* $Id: dMatrixEXT.cpp,v 1.16 2005-07-29 03:09:33 paklein Exp $ */
/* created: paklein (03/06/1998)                                          */
#include "dMatrixEXT.h"
#include "toolboxConstants.h"
#include "LAdMatrixT.h"
#include "dArray2DT.h"
#include <math.h>
#include "dTensor4DT.h"
//#include "NRUTIL.H"

using namespace Tahoe;
const char caller[] = "dMatrixEXT";

/* Numerical Recipies macros */
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SWAP(g,h) {y=(g);(g)=(h);(h)=y;}

static inline double FMAX(double a, double b) { return (a > b) ? a : b; };
static inline int IMIN(int a, int b) { return (a > b) ? b : a; };

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

dMatrixEXT::dMatrixEXT(int squaredim, const double* p):
	dMatrixT(squaredim,squaredim,p)
{
	fworkspace.Dimension(2*squaredim);

	/* set pointers */
	v1 = fworkspace.Pointer();
	v2 = v1 + squaredim;
}

/* post constructor (re-)dimensioning */
void dMatrixEXT::Dimension(int squaredim)
{
	/* inherited */
	dMatrixT::Dimension(squaredim);
	
	/* temp space */
	fworkspace.Dimension(2*squaredim);

	/* set pointers */
	v1 = fworkspace.Pointer();
	v2 = v1 + squaredim;
}

/* diagonalize (using symmetric QR algorithm) */
int dMatrixEXT::Diagonalize(dArrayT& eigs)
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (fRows != fCols) ExceptionT::GeneralFail(caller);
	if (fRows != eigs.Length()) ExceptionT::SizeMismatch(caller);
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
	if (fRows != fCols) ExceptionT::GeneralFail(caller);
	if (fRows != eigenvector.Length()) ExceptionT::SizeMismatch(caller);
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

/* generate singular value decomposition of *this = U*W*V^T */
void dMatrixEXT::Compute_SVD(dMatrixT& U, dArrayT& W, dMatrixT& V, double threshold, 
	int max_its) const
{
#if __option(extended_errorcheck)
	if (fRows != U.Rows()  || 
	    fCols != U.Cols()  ||
	    fCols != W.Length() ||
	    fCols != V.Rows()  ||
	    fCols != V.Cols()) ExceptionT::SizeMismatch(caller);
	if (threshold > 1 || threshold < 0) ExceptionT::GeneralFail(caller);
#endif

	/* copy */
	U = *this;
	
	/* clear */
	W = 0.0;
	V = 0.0;

	/* call driver */
	if (!svdcmp(U.Pointer(), fRows, fCols, W.Pointer(), V.Pointer(), v1, max_its))
		ExceptionT::GeneralFail("dMatrixEXT::Compute_SVD", "failed to converge");
	
	/* find max singular value */
	double w_max = W.Max();
	
	/* drop smaller values */
	W.Chop(w_max*threshold);
}

/* back substitute given a decomposition computed with dMatrixEXT::Compute_SVD */
void dMatrixEXT::BackSubstitute_SVD(const dMatrixT& U, const dArrayT& W, 
	const dMatrixT& V, dArrayT& RHS) const
{
#if __option(extended_errorcheck)
	if (fRows != U.Rows()  || 
	    fCols != U.Cols()  ||
	    fCols != W.Length() ||
	    fRows != RHS.Length()) ExceptionT::SizeMismatch(caller);
#endif

	/* call driver */
	svbksb(U.Pointer(), W.Pointer(), V.Pointer(), fRows, fCols, RHS.Pointer(), v1, v2);

	/* copy in */
	RHS.Copy(v1);
}

/*************************************************************************
* Private
*************************************************************************/

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

double dMatrixEXT::pythag(double a, double b) const
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
	int m,l,iter = -1,i;
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
					ExceptionT::GeneralFail("dMatrixEXT::tqli", "no convergence after %d iterations", 3*n);
					
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

/* finds eigenvalues of a general matrix by first reducing to Hessian form*/
/* Was having problems with this so not in use, along with numerical recipes
 *elmes and hqr*/
void dMatrixEXT::eigvalfinder (dMatrixEXT& matrix, dArrayT& realev, dArrayT& imev)
{

#if __option(extended_errorcheck)
	if (fRows != fCols) ExceptionT::GeneralFail(caller);
	if (fRows != realev.Length()) ExceptionT::SizeMismatch(caller);
        if (fRows != imev.Length()) ExceptionT::SizeMismatch(caller);
#endif

	dMatrixEXT matrixCopy(fRows);
	matrixCopy = matrix;

elmhes(matrixCopy, fRows);

//matrixCopy(0,2) = 0.0;

hqr(matrixCopy, fRows, realev, imev);

}

void dMatrixEXT::elmhes(dMatrixEXT& a,int n)
  //float **a;
  //int n;
{
	int m,j,i;
	double y,x;

	// for (m=2;m<n<m++)
	for (m=1;m<n-1;m++) {
		x=0.0;
		i=m;
		//	for (j=m;j<=n;j++) {
		for (j=m;j<=n-1;j++) {
			if (fabs(a(j,m-1)) > fabs(x)) {
				x=a(j,m-1);
				i=j;
			}
		}
		if (i != m) {
		  //for (j=m-1;j<=n;j++) SWAP(a(i,j),a(m,j))
		  //for (j=1;j<=n;j++) SWAP(a(j,i),a(j,m))
		  for (j=m-1;j<=n-1;j++) SWAP(a(i,j),a(m,j))
		  for (j=0;j<=n-1;j++) SWAP(a(j,i),a(j,m))
		}
		if (x) {
		  //for (i=m+1;i<=n;i++) {
		  for (i=m+1;i<=n-1;i++) {
				if ((y=a(i,m-1)) != 0.0) {
					y /= x;
					a(i,m-1)=y;
					//for (j=m;j<=n;j++)
					for (j=m;j<=n-1;j++)
						a(i,j) -= y*a(m,j);
					//for (j=1;j<=n;j++)
					for (j=0;j<=n-1;j++)
						a(j,m) += y*a(j,i);
				}
			}
		}
	}
}


void dMatrixEXT::hqr(dMatrixEXT& a, int n, dArrayT& wr, dArrayT& wi)
  //float **a,wi[],wr[];
  //int n;
{
	int nn,m,l,k,j,its,i,mmin;
	double z,y,x,w,v,u,t,s,r = 0.,q = 0.,p = 0.,anorm;

	//anorm=fabs(a(1,1));
	anorm=fabs(a(0,0));	
	//for (i=2;i<=n;i++)
	//	for (j=(i-1);j<=n;j++)
	for (i=1;i<n;i++)
		for (j=(i-1);j<n;j++)
			anorm += fabs(a(i,j));
	nn=n-1;
	t=0.0;

	//while (nn >= 1) {
	while (nn >= 0) {
		its=0;
		do {
		      //for (l=nn;l>=2;l--) {
			for (l=nn;l>=1;l--) {
				s=fabs(a(l-1,l-1))+fabs(a(l,l));
				if (s == 0.0) s=anorm;
				if ((double)(fabs(a(l,l-1)) + s) == s) break;
			}
			x=a(nn,nn);
			if (l == nn) {
				wr[nn]=x+t;
				wi[nn--]=0.0;
			} else {
				y=a(nn-1,nn-1);
				w=a(nn,nn-1)*a(nn-1,nn);
				if (l == (nn-1)) {
					p=0.5*(y-x);
					q=p*p+w;
					z=sqrt(fabs(q));
					x += t;
					if (q >= 0.0) {
						z=p+SIGN(z,p);
						wr[nn-1]=wr[nn]=x+z;
						if (z) wr[nn]=x-w/z;
						wi[nn-1]=wi[nn]=0.0;
					} else {
						wr[nn-1]=wr[nn]=x+p;
						wi[nn-1]= -(wi[nn]=z);
					}
					nn -= 2;
				} else {
					if (its == 30)
					  ExceptionT::GeneralFail("dMatrixEXT::hqr", "%d iterations", its);

					if (its == 10 || its == 20) {
						t += x;
						//for (i=1;i<=nn;i++) a(i,i) -= x;
						for (i=0;i<nn;i++) a(i,i) -= x;	
						s=fabs(a(nn,nn-1))+fabs(a(nn-1,nn-2));
						y=x=0.75*s;
						w = -0.4375*s*s;
					}
					++its;
					// should be okay
					for (m=(nn-2);m>=l;m--) {
						z=a(m,m);
						r=x-z;
						s=y-z;
						p=(r*s-w)/a(m+1,m)+a(m,m+1);
						q=a(m+1,m+1)-z-r-s;
						r=a(m+2,m+1);
						s=fabs(p)+fabs(q)+fabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l) break;
						u=fabs(a(m,m-1))*(fabs(q)+fabs(r));
						v=fabs(p)*(fabs(a(m-1,m-1))+fabs(z)+fabs(a(m+1,m+1)));
						if ((float)(u+v) == v) break;
					}
					// should be okay  for (i=m+2;i<=nn;i++) {
					for (i=m+2;i<=nn;i++) {
						a(i,i-2)=0.0;
						if (i != (m+2)) a(i,i-3)=0.0;
					}
					//  should be okay for (k=m;k<=nn-1;k++) {
					for (k=m;k<=nn-1;k++) {
						if (k != m) {
							p=a(k,k-1);
							q=a(k+1,k-1);
							r=0.0;
							if (k != (nn-1)) r=a(k+2,k-1);
							if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0) {
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s=SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
							if (k == m) {
								if (l != m)
								a(k,k-1) = -a(k,k-1);
							} else
								a(k,k-1) = -s*x;
							p += s;
							x=p/s;
							y=q/s;
							z=r/s;
							q /= p;
							r /= p;
							// should be okay
							for (j=k;j<=nn;j++) {
								p=a(k,j)+q*a(k+1,j);
								if (k != (nn-1)) {
									p += r*a(k+2,j);
									a(k+2,j) -= p*z;
								}
								a(k+1,j) -= p*y;
								a(k,j) -= p*x;
							}
							mmin = nn<k+3 ? nn : k+3;
							// should be okay l not 1
							for (i=l;i<=mmin;i++) {
								p=x*a(i,k)+y*a(i,k+1);
								if (k != (nn-1)) {
									p += z*a(i,k+2);
									a(i,k+2) -= p*r;
								}
								a(i,k+1) -= p*q;
								a(i,k) -= p;
							}
						}
					}
				}
			}
		} while (l < nn-1);
	}
}


/* void svdcmp(double **a, int m, int n, double w[], double **v) */
int dMatrixEXT::svdcmp(double* a, int m, int n, double* w, double* v, double* rv1, int max_its) const
{
//	double pythag(double a, double b);
	int flag,i,its,j,jj,k,l=1,nm=0;
	double anorm,c,f,g,h,s,scale,x,y,z;

//	rv1=vector(1,n);
	g=scale=anorm=0.0;
//	for (i=1;i<=n;i++) {
	for (i=0;i<n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
//		if (i <= m) {
		if (i < m) {
//			for (k=i;k<=m;k++) scale += fabs(a[i*m + k]);
			for (k=i;k<m;k++) scale += fabs(a[i*m + k]);
			if (scale) {
//				for (k=i;k<=m;k++) {
				for (k=i;k<m;k++) {
					a[i*m + k] /= scale;
					s += a[i*m + k]*a[i*m + k];
				}
				f=a[i*m + i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i*m + i]=f-g;
//				for (j=l;j<=n;j++) {
				for (j=l;j<n;j++) {
//					for (s=0.0,k=i;k<=m;k++) s += a[i*m + k]*a[j*m + k];
					for (s=0.0,k=i;k<m;k++) s += a[i*m + k]*a[j*m + k];
					f=s/h;
//					for (k=i;k<=m;k++) a[j*m + k] += f*a[i*m + k];
					for (k=i;k<m;k++) a[j*m + k] += f*a[i*m + k];
				}
//				for (k=i;k<=m;k++) a[i*m + k] *= scale;
				for (k=i;k<m;k++) a[i*m + k] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
//		if (i <= m && i != n) {
		if (i < m && i != n) {
//			for (k=l;k<=n;k++) scale += fabs(a[k*m + i]);
			for (k=l;k<n;k++) scale += fabs(a[k*m + i]);
			if (scale) {
//				for (k=l;k<=n;k++) {
				for (k=l;k<n;k++) {
					a[k*m + i] /= scale;
					s += a[k*m + i]*a[k*m + i];
				}
				f=a[l*m + i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[l*m + i]=f-g;
//				for (k=l;k<=n;k++) rv1[k]=a[k*m + i]/h;
				for (k=l;k<n;k++) rv1[k]=a[k*m + i]/h;
//				for (j=l;j<=m;j++) {
				for (j=l;j<m;j++) {
//					for (s=0.0,k=l;k<=n;k++) s += a[k*m + j]*a[k*m + i];
					for (s=0.0,k=l;k<n;k++) s += a[k*m + j]*a[k*m + i];
//					for (k=l;k<=n;k++) a[k*m + j] += s*rv1[k];
					for (k=l;k<n;k++) a[k*m + j] += s*rv1[k];
				}
//				for (k=l;k<=n;k++) a[k*m + i] *= scale;
				for (k=l;k<n;k++) a[k*m + i] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
//	for (i=n;i>=1;i--) {
	for (i=n-1;i>-1;i--) {
//		if (i < n) {
		if (i < n-1) {
			if (g) {
//				for (j=l;j<=n;j++)
				for (j=l;j<n;j++)
					v[i*n + j]=(a[j*m + i]/a[l*m + i])/g;
//				for (j=l;j<=n;j++) {
				for (j=l;j<n;j++) {
//					for (s=0.0,k=l;k<=n;k++) s += a[k*m + i]*v[j*n + k];
					for (s=0.0,k=l;k<n;k++) s += a[k*m + i]*v[j*n + k];
//					for (k=l;k<=n;k++) v[j*n + k] += s*v[i*n + k];
					for (k=l;k<n;k++) v[j*n + k] += s*v[i*n + k];
				}
			}
//			for (j=l;j<=n;j++) v[j*n + i]=v[i*n + j]=0.0;
			for (j=l;j<n;j++) v[j*n + i]=v[i*n + j]=0.0;
		}
		v[i*n + i]=1.0;
		g=rv1[i];
		l=i;
	}
//	for (i=IMIN(m,n);i>=1;i--) {
	for (i=IMIN(m,n) - 1; i> -1;i--) {
		l=i+1;
		g=w[i];
//		for (j=l;j<=n;j++) a[j*m + i]=0.0;
		for (j=l;j<n;j++) a[j*m + i]=0.0;
		if (g) {
			g=1.0/g;
//			for (j=l;j<=n;j++) {
			for (j=l;j<n;j++) {
//				for (s=0.0,k=l;k<=m;k++) s += a[i*m + k]*a[j*m + k];
				for (s=0.0,k=l;k<m;k++) s += a[i*m + k]*a[j*m + k];
				f=(s/a[i*m + i])*g;
//				for (k=i;k<=m;k++) a[j*m + k] += f*a[i*m + k];
				for (k=i;k<m;k++) a[j*m + k] += f*a[i*m + k];
			}
//			for (j=i;j<=m;j++) a[i*m + j] *= g;
			for (j=i;j<m;j++) a[i*m + j] *= g;
//		} else for (j=i;j<=m;j++) a[i*m + j]=0.0;
		} else for (j=i;j<m;j++) a[i*m + j]=0.0;
		++a[i*m + i];
	}
//	for (k=n;k>=1;k--) {
	for (k=n-1;k > -1;k--) {
		for (its=1;its<=max_its;its++) {
			flag=1;
//			for (l=k;l>=1;l--) {
			for (l=k;l > -1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) { //keep
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
//					for (j=1;j<=m;j++) {
					for (j=0;j<m;j++) {
						y=a[nm*m + j];
						z=a[i*m + j];
						a[nm*m + j]=y*c+z*s;
						a[i*m + j]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
//					for (j=1;j<=n;j++) v[k*n + j] = -v[k*n + j];
					for (j=0;j<n;j++) v[k*n + j] = -v[k*n + j];
				}
				break;
			}
			if (its == max_its)
			{
				cout << " dMatrixEXT:svdcmp: no convergence in " << max_its << " svdcmp iterations" << endl;
				return 0;
			}
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) { //keep
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
//				for (jj=1;jj<=n;jj++) {
				for (jj=0;jj<n;jj++) {
					x=v[j*n + jj];
					z=v[i*n + jj];
					v[j*n + jj]=x*c+z*s;
					v[i*n + jj]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
//				for (jj=1;jj<=m;jj++) {
				for (jj=0;jj<m;jj++) {
					y=a[j*m + jj];
					z=a[i*m + jj];
					a[j*m + jj]=y*c+z*s;
					a[i*m + jj]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
//	free_vector(rv1,1,n);
	return 1;
}

/* void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]) */
void dMatrixEXT::svbksb(const double* u, const double* w, const double* v, int m, int n, 
	const double* b, double* x, double* tmp) const
{
	int jj,j,i;
	double s;

	//tmp=vector(1,n);
//	for (j=1;j<=n;j++) {
	for (j=0;j<n;j++) {
		s=0.0;
		if (w[j]) {
//			for (i=1;i<=m;i++) s += u[j*m + i]*b[i];
			for (i=0;i<m;i++) s += u[j*m + i]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
//	for (j=1;j<=n;j++) {
	for (j=0;j<n;j++) {
		s=0.0;
//		for (jj=1;jj<=n;jj++) s += v[jj*n + j]*tmp[jj];
		for (jj=0;jj<n;jj++) s += v[jj*n + j]*tmp[jj];
		x[j]=s;
	}
	//free_vector(tmp,1,n);
}



/* find eigenvalues of 3x3 matrix explicitly */
void dMatrixEXT::eigenvalue3x3(dMatrixEXT& J, dArrayT& reroot, dArrayT& imroot)
{
	if (fRows != 3) ExceptionT::GeneralFail("dMatrixEXT::eigenvalue3x3", "matrix must be dimension 3 not %d", fRows);

dMatrixEXT matrix(3);
  int n, flag;  // counter
double p,q,r; //coeffs of cubic
double a,b;   //coeffs of depressed cubic
double s1,s2,s3,s4,s5; //calculations numbers
 double radius, cangle; //polar form for cube root
double reA, reB, imA, imB; //real and imag parts of coeffs
double tol=10e-10;

double sqrt3=sqrt(3.0);
double pi = acos(0.0);

matrix=J;
double det=J.Det();

// calculate coefficients of determinant cubic y^3+p*y^2+q*y+r=0
//see Beyer, Standard Mathematic Tables for more info on solution of cubic

p=-1*(J(0,0)+J(1,1)+J(2,2));
q=J(1,1)*J(2,2)+J(2,2)*J(0,0)+J(0,0)*J(1,1)-J(1,2)*J(2,1)-J(0,1)*J(1,0)-J(0,2)*J(2,0);
r=-1*J.Det();



// coeffs of depressed cubic x^3 + a*x +b = 0

a=(3*q-p*p)/3;
b=(2*p*p*p-9*p*q+27*r)/27;


s1=-0.5*b;
s3=b*b/4+a*a*a/27;
 if (s3>=0){

   // s3 > 0 => two complex conjugate roots, one real root
  // if (s3>0.0) 
    // cout << "Warning: complex roots detected \n";
   //else
     //cout << "Repeated eigenvalue detected \n";

s2=sqrt(s3);

s4=s1+s2;
s5=s1-s2;

 if (s4 != 0)
   reA=fabs(s4)/s4*pow(fabs(s4),1.0/3.0);
 else
   reA=0;

if (s5!=0)
  reB=fabs(s5)/s5*pow(fabs(s5),1.0/3.0);
else
  reB=0;

reroot[0]=reA+reB;
imroot[0]=0;
reroot[1]=-.5*(reA+reB);
imroot[1]=.5*(reA-reB)*sqrt3;
reroot[2]=reroot[1];
imroot[2]=-1*imroot[1];

reroot[0]=reroot[0]-p/3;
reroot[1]=reroot[1]-p/3;
reroot[2]=reroot[2]-p/3;

 }
else
{

  //case of 3 real roots

s2=sqrt(-1*s3);

radius=sqrt(s1*s1-s3);
radius=pow(radius,1.0/3.0);

cangle=atan(s2/s1);
cangle=cangle/3;

n=0;
flag=0;

 while ( flag==0 && n<3){

   //must cycle through all possible complex cube roots of equation

reA=radius*cos(cangle);
imA=radius*sin(cangle);

 reB=reA;            //=-1*radius*cos(-1*cangle);
 imB=-1*imA;                    //radius*sin(-1*cangle);

 reroot[0]=reA+reB;              
 imroot[0]=0;
 reroot[1]=.5*((imB-imA)*sqrt3-reA-reB);
 imroot[1]=0;
 reroot[2]=.5*((imA-imB)*sqrt3-reA-reB);
 imroot[2]=0;

reroot[0]=reroot[0]-p/3;
reroot[1]=reroot[1]-p/3;
reroot[2]=reroot[2]-p/3;

 matrix(0,0)=J(0,0)-reroot[0];
 matrix(1,1)=J(1,1)-reroot[0];
 matrix(2,2)=J(2,2)-reroot[0];


 // check to see if proper roots found
 
 if (fabs(matrix.Det()/det) < tol || (det==0  && fabs(double(matrix.Det()<tol)) ))
   {
     matrix(0,0)=J(0,0)-reroot[1];
     matrix(1,1)=J(1,1)-reroot[1];
     matrix(2,2)=J(2,2)-reroot[1];
     
     if (fabs(matrix.Det()/det) < tol || (det==0  && fabs(double(matrix.Det()<tol)) ))
       {
	 
	 matrix(0,0)=J(0,0)-reroot[2];
	 matrix(1,1)=J(1,1)-reroot[2];
	 matrix(2,2)=J(2,2)-reroot[2];

	 if (fabs(matrix.Det()/det) < tol)
	   flag=1;
       }
   }
 
 cangle=cangle+2*pi/3;
n++;
//if (n==3 && flag==0)
 // cout << "eigenvalues not found \n";

 }

}



}

// finds eigenvector of 3x3 matrix given its eigenvalue
// may expnad to return multiple vector for same value if close

void dMatrixEXT::eigenvector3x3(dMatrixEXT& J, double value, int numvector,  dArrayT& vector, dArrayT& vector2, dArrayT& vector3)
{

double tol =1.0e-10;
double det;
double Jdet=J.Det();
dMatrixEXT matrix(3);
//int rnk;

matrix=J;


matrix(0,0)=matrix(0,0)-value;
matrix(1,1)=matrix(1,1)-value;
matrix(2,2)=matrix(2,2)-value;


if (fabs(matrix.Det()/Jdet)>tol && fabs(matrix.Det())>tol)
 {
   //cout << "Value entered not an eigenvalue \n";
   numvector=0;

   //dummy- i.e. not real eigenvalues
   vector[0]=1;
   vector[1]=0;
   vector[2]=0;
   
   vector2[0]=0;
   vector2[1]=1;
   vector2[2]=0;
	 
   vector3[0]=0;
   vector3[1]=0;
   vector3[2]=1;
       }
else
  {
    numvector =1; // finding one eigenvector
 Jdet = pow(fabs(Jdet), 2.0/3.0); //puts on same order of magnitude as 2x2

 //solve for eigenvector
if (fabs((matrix(0,0)*matrix(1,1)-matrix(1,0)*matrix(0,1))/Jdet)>tol && fabs(matrix(0,0)*matrix(1,1)-matrix(1,0)*matrix(0,1))>tol)
  {
vector[2]=1;

det=matrix(0,0)*matrix(1,1)-matrix(1,0)*matrix(0,1);

vector[0]=1/det*(matrix(0,1)*matrix(1,2)-matrix(1,1)*matrix(0,2));
vector[1]=1/det*(matrix(1,0)*matrix(0,2)-matrix(0,0)*matrix(1,2));

/*dummy*/
vector2[0]=1;
vector2[1]=0;
vector2[2]=0;

vector3[0]=1;
vector3[1]=0;
vector3[2]=0;
  }
 else

   // used if determinant of 1st submatrix is close to 0
   if (fabs((matrix(0,0)*matrix(2,2)-matrix(2,0)*matrix(0,2))/Jdet) > tol && (fabs(matrix(0,0)*matrix(2,2)-matrix(2,0)*matrix(0,2))>tol ))
     {
       //cout << " stage 2 \n";
       
       vector[1]=1;
       
       det=matrix(0,0)*matrix(2,2)-matrix(2,0)*matrix(0,2);
       
       vector[0]=1/det*(matrix(0,2)*matrix(2,1)-matrix(2,2)*matrix(0,1));
       vector[2]=1/det*(matrix(2,0)*matrix(0,1)-matrix(0,0)*matrix(2,1));
       
       /*dummy*/
       vector2[0]=1;
       vector2[1]=0;
       vector2[2]=0;
       
       vector3[0]=1;
       vector3[1]=0;
       vector3[2]=0;

       //cout << "vector = \n";
       //cout << vector << endl;

     }
   else
     //used if 2nd submatrix also close to 0
     if (fabs((matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2))/Jdet) > tol && fabs(matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)) > tol)
     {
       // cout << " stage 3 \n";

       vector[0]=1;
       
       det=matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2);
       
       vector[1]=1/det*(matrix(1,2)*matrix(2,0)-matrix(2,2)*matrix(1,0));
       vector[2]=1/det*(matrix(2,1)*matrix(1,0)-matrix(1,1)*matrix(2,0));

       /*dummy*/
       vector2[0]=1;
       vector2[1]=0;
       vector2[2]=0;
       
       vector3[0]=1;
       vector3[1]=0;
       vector3[2]=0;

     }
     else
       {
	 // vector not found-additional possible submatrices not yet implemented
	 //cout << "vector seems to exist but is not found \n";
	 /* dummy*/

	 //cout << " stage 4 \n";
	
	 vector[0]=1;
	 vector[1]=0;
	 vector[2]=0;
	 
	 vector2[0]=1;
	 vector2[1]=0;
	 vector2[2]=0;
	 
	 vector3[0]=1;
	 vector3[1]=0;
	 vector3[2]=0;
       }
  }
  

/*normalize vector*/

double norm;

norm=sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);
vector[0]=vector[0]/norm;
vector[1]=vector[1]/norm;
vector[2]=vector[2]/norm;

norm=sqrt(vector2[0]*vector2[0]+vector2[1]*vector2[1]+vector2[2]*vector2[2]);
vector2[0]=vector2[0]/norm;
vector2[1]=vector2[1]/norm;
vector2[2]=vector2[2]/norm;

norm=sqrt(vector3[0]*vector3[0]+vector3[1]*vector3[1]+vector3[2]*vector3[2]);
vector3[0]=vector3[0]/norm;
vector3[1]=vector3[1]/norm;
vector3[2]=vector3[2]/norm;

// cout << "vector = \n";
// cout << vector << endl;

}


/*forms acoustic tensor from rank 4 tangent modulus, normal */
//void dMatrixEXT::formacoustictensor(dMatrixEXT& A, double C [3] [3] [3] [3], dArrayT& normal)
void dMatrixEXT::formacoustictensor(dMatrixEXT& A, dTensor4DT& C, dArrayT& normal)
{
 int k,l,m,n; //counters

	// initialize acoustic tensor A
	A = 0.0;
	
	    
 // A=normal*C*normal
              
 for (m=0;m<3;m++){
   for (n=0;n<3;n++){
     for (k=0;k<3;k++){
       for (l=0;l<3;l++){
	 //A (m,n)+=normal[k]*C [k] [m] [n] [l]*normal[l]; 
	 A (m,n)+=normal[k]*C(k,m,n,l)*normal[l]; 
       }
     }
   }
 }

}
