#include "CubicSplineWindowT.h"
#include "ExceptionT.h"
#include <math.h>

using namespace Tahoe;

const double sqrtPi = sqrt(acos(-1.0));
static double Max(double a, double b) { return (a > b) ? a : b; };

/* constructor */
CubicSplineWindowT::CubicSplineWindowT(double dilation_scaling):
	fDilationScaling(dilation_scaling)
{
	if (fDilationScaling < 0.0) 
		ExceptionT::BadInputValue("CubicSplineWindowT::CubicSplineWindowT");
}

/* "synchronization" of nodal field parameters. */
void CubicSplineWindowT::SynchronizeSupportParameters(dArray2DT& params_1, 
	dArray2DT& params_2) const
{
	/* should be over the same global node set (marked by length) */
	if (params_1.Length() != params_2.Length() ||
	    params_1.MinorDim() != NumberOfSupportParameters() ||
	    params_2.MinorDim() != NumberOfSupportParameters())
		ExceptionT::SizeMismatch("CubicSplineWindowT::SynchronizeSupportParameters");
		
	/* "synchronize" means take max of dmax */
	double* p1 = params_1.Pointer();
	double* p2 = params_2.Pointer();
	int length = params_1.Length();
	for (int i = 0; i < length; i++)
	{
		*p1 = *p2 = Max(*p1, *p2);
		p1++; p2++;
	}
}

void CubicSplineWindowT::WriteParameters(ostream& out) const
{
	/* window function parameters */
	out << " Dilation scaling factor . . . . . . . . . . . . = " << fDilationScaling << '\n';
}

/* Single point evaluations */
bool CubicSplineWindowT::Window(const dArrayT& x_n, const dArrayT& param_n, const dArrayT& x,
		int order, double& w, dArrayT& Dw, dSymMatrixT& DDw, dArrayT& DDDw) 
{
	/* outside the range of influence */
	if (!CubicSplineWindowT::Covers(x_n, x, param_n))
	{
		w = 0.0;
		if (order > 0)
		{
			Dw = 0.0;
			if (order > 1)
			{
				DDw = 0.0;
				if (order > 2) // kyonten
					DDDw = 0.0;	
			}
		}
		
		/* does not cover */
		return false;
	}
  	else
  	{
		/* distance */
		Dw.DiffOf(x_n, x);
		double dist = Dw.Magnitude();
		double a = 0.5*param_n[0]*fDilationScaling; /* window defined over -2 < x < 2 */
		double r = dist/a;
		if (r > 2.0) {
			w = 0.0;
			Dw = 0.0;
			DDw = 0.0;
			DDDw = 0.0; // kyonten
		} else if (r > 1.0) {
			double dr = 2.0 - r;
			w = 1.0/(6.0*a)*dr*dr*dr;
			if (order > 0) {
				double dw = -1.0/(2.0*a)*dr*dr;
				if (order > 1) {
					double ddw = dr/a;
					DDw.Outer(Dw, (ddw/a - dw/dist)/(dist*dist*a));
					DDw.PlusIdentity(dw/(dist*a));
					if (order > 2) // kyonten
	  				{ 
	  					/*  NOTE: DDDw is a [nsd*nsd]x[nsd] matrix (outer product of 3 vectors). 
	  				    	using symmetry we have only 4 components in 2D or 10 components in 3D
	  					*/
	  
						int nsd = x.Length(); 
	  					/* work space */
	  					dSymMatrixT dddw1(nsd);
	  					dMatrixT dddw2(nsd),  dddw3(nsd);
	  					dArrayT dddw1_dia(nsd), I(nsd), colA(nsd), colB(nsd);
	  					dArrayT DDDw2(DDDw.Length());
	  					int offset_length;
	  					
	  					/* initialize */
	  					DDDw2 = 0.0;
	  					
	  					dddw1.Outer(Dw);
	  					/* collect diagonal terms */
	  					for (int i = 0; i < nsd; i++)
	  					{
	  						dddw1_dia[i] = dddw1(i,i); 
	  						I[i] = 1.0;
	  					}
	  					dddw2.Outer(Dw, dddw1_dia);
	  					dddw3.Outer(Dw, I);
	  					/* collect the components of dddw2 and dddw3 in arrays */
	  					for (int i = 0; i < nsd; i++) {
	  						dddw2.CopyColumn(i, colA);
	  						dddw3.CopyColumn(i, colB);
	  						offset_length = i*colA.Length();
	  						if (i == 0) {
	  							DDDw.CopyIn(0, colA);
	  							DDDw2.CopyIn(0, colB);
	  						}
	  						else {
	  							DDDw.CopyIn(offset_length, colA);
	  							DDDw2.CopyIn(offset_length, colB);
	  						}
	  					
	  					}
	  					/* 3D: add the 10th component */
	  					if (nsd == 3)
	  						DDDw[nsd*nsd-1] = dddw1[5]*Dw[2];
	  						
	  					DDDw *= (-dist/(a*a*a*a) - 3.0*(dr)/(a*a*a) - 3.0*dr*dr/(2.0*a*a*dist));
	  					DDDw /= (dist*dist*dist*dist);
	  					DDDw2 *= 3.0*(dr/(a*a*a) + dr*dr/(2.0*a*a*dist));
	  					DDDw /= (dist*dist);
	  					DDDw += DDDw2;
	  				} // if(order > 2)
				}
				Dw *= dw/(dist*a);
			}
		} else {
			w = (2.0/3.0 + r*r*(r/2.0 - 1.0))/a;
			if (order > 0) {
				//double dw = r*(3.0*r/2.0 - 2.0)/a; re-write to avoid division by zero when r = 0.
				double dw_by_r = (3.0*r/2.0 - 2.0)/a;
				if (order > 1) {
					double ddw = (3.0*r - 2.0)/a;
					DDw.Outer(Dw, (ddw/a - dw_by_r*r/dist)/(dist*dist*a));
					DDw.PlusIdentity(dw_by_r*r/(dist*a));
					if (order > 2) // kyonten
	  				{
	  					int nsd = x.Length(); 
	  					/* work space */
	  					dSymMatrixT dddw1(nsd);
	  					dMatrixT dddw2(nsd),  dddw3(nsd);
	  					dArrayT dddw1_dia(nsd), I(nsd), colA(nsd), colB(nsd);
	  					dArrayT DDDw2(DDDw.Length()), DDDw3(DDDw.Length());
	  					int offset_length;
	  					
	  					/* initialize */
	  					DDDw2 = 0.0;
	  					DDDw3 = 0.0;
	  					
	  					dddw1.Outer(Dw);
	  					/* collect diagonal terms */
	  					for (int i = 0; i < nsd; i++)
	  					{
	  						dddw1_dia[i] = dddw1(i,i); 
	  						I[i] = 1.0;
	  					}
	  					dddw2.Outer(Dw, dddw1_dia);
	  					dddw3.Outer(Dw, I);
	  					/* collect the components of dddw2 and dddw3 in arrays */
	  					for (int i = 0; i < nsd; i++) {
	  						dddw2.CopyColumn(i, colA);
	  						dddw3.CopyColumn(i, colB);
	  						offset_length = i*colA.Length();
	  						if (i == 0) {
	  							DDDw.CopyIn(0, colA);
	  							DDDw2.CopyIn(0, colB);
	  							DDDw3.CopyIn(0, colB);
	  						}
	  						else {
	  							DDDw.CopyIn(offset_length, colA);
	  							DDDw2.CopyIn(offset_length, colB);
	  							DDDw3.CopyIn(offset_length, colB);
	  						}
	  					
	  					}
	  					/* 3D: add the 10th component */
	  					if (nsd == 3)
	  						DDDw[nsd*nsd-1] = dddw1[5]*Dw[2];
	  						
	  					DDDw  *= (3.0*dist/(a*a*a) - 3.0*r/(2*a*a) - dw_by_r*r/dist);
	  					DDDw2 *= 2.0*(ddw/a - dw_by_r*r/dist);
	  					DDDw += DDDw2;
	  					DDDw /= (a*dist*dist*dist*dist);
	  					DDDw3 *= (ddw/a - dw_by_r*r/dist);
	  					DDDw3 /= (a*dist*dist);   	
	  					DDDw += DDDw3;
	  				} // (order > 2)
				}
				Dw *= dw_by_r/(a*a);
			}
		}
	
		/* does cover */
		return true;	
	}
}

/* multiple point calculations */
int CubicSplineWindowT::Window(const dArray2DT& x_n, const dArray2DT& param_n, const dArrayT& x,
		int order, dArrayT& w, dArray2DT& Dw, dArray2DT& DDw, dArray2DT& DDDw)
{
	/* allocate */
	int nsd = x.Length();
	fNSD.Dimension(nsd);
	fNSDsym.Dimension(nsd);
	
	if (nsd == 3)
		fNSDArray.Dimension(nsd*nsd+1);
	else
		fNSDArray.Dimension(nsd*nsd);

	/* work space */
	dArrayT x_node, param_node;
	int count = 0;    
	int numwindowpoints = x_n.MajorDim(); 

	for (int i = 0; i < numwindowpoints; i++)
	{
		/* collect nodal values */
		x_n.RowAlias(i, x_node);
		param_n.RowAlias(i, param_node);
      
		if (CubicSplineWindowT::Window(x_node, param_node, x, order, w[i], fNSD, fNSDsym, fNSDArray))
			count++;

		/* store derivatives */
		if (order > 0)
		{
			Dw.SetColumn(i, fNSD);
			if (order > 1)
			{
				DDw.SetColumn(i, fNSDsym);
				if (order > 2) // kyonten
					DDDw.SetColumn(i, fNSDArray); 
			}
		}
	}
	
	return count;
}

bool CubicSplineWindowT::Covers(const dArrayT& x_n, const dArrayT& x, const dArrayT& param_n) const
{
	double dist = dArrayT::Distance(x, x_n);
	return (dist < param_n[0]*fDilationScaling);
}

int CubicSplineWindowT::Covers(const dArray2DT& x_n, const dArrayT& x, 
	const dArray2DT& param_n, ArrayT<bool>& covers) const
{
	int count = 0;
	int numwindowpoints = x_n.MajorDim();
	for (int i = 0; i < numwindowpoints; i++)
	{
		double dist = dArrayT::Distance(x, x_n);
		if (dist < param_n[0]*fDilationScaling) {
			count++;
			covers[i] = true;
		} 
		else
			covers[i] = false;
	}
	
	return count;
}

/* spherical upport size */
double CubicSplineWindowT::SphericalSupportSize(const dArrayT& param_n) const
{
#if __option(extended_errorcheck)
	if (param_n.Length() != 1) ExceptionT::GeneralFail("CubicSplineWindowT::SphericalSupportSize");
#endif
	return 1.0*fDilationScaling*param_n[0];
}

/* rectangular support size */
void CubicSplineWindowT::RectangularSupportSize(const dArrayT& param_n, dArrayT& support_size) const
{
	/* same in all dimensions */
	support_size = SphericalSupportSize(param_n);
}

/* spherical support sizes in batch */
void CubicSplineWindowT::SphericalSupportSize(const dArray2DT& param_n, ArrayT<double>& support_size) const
{
#if __option(extended_errorcheck)
	if (param_n.MinorDim() != 1 ||
	    param_n.MajorDim() != support_size.Length()) 
		ExceptionT::GeneralFail("CubicSplineWindowT::SphericalSupportSize");
#endif

	dArrayT tmp;
	tmp.Alias(support_size);
	tmp.SetToScaled(1.0*fDilationScaling, param_n);
}
