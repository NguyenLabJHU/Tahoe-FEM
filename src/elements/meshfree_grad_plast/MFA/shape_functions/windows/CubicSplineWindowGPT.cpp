#include "CubicSplineWindowGPT.h"
#include "ExceptionT.h"
#include <math.h>


using namespace Tahoe;

const double sqrtPi = sqrt(acos(-1.0));
static double Max(double a, double b) { return (a > b) ? a : b; };

/* constructor */
CubicSplineWindowGPT::CubicSplineWindowGPT(double dilation_scaling):
	fDilationScaling(dilation_scaling)
{
	if (fDilationScaling < 0.0) throw ExceptionT::kBadInputValue;
}

/* "synchronization" of nodal field parameters. */
void CubicSplineWindowGPT::SynchronizeSupportParameters(dArray2DT& params_1, 
	dArray2DT& params_2) const
{
	/* should be over the same global node set (marked by length) */
	if (params_1.Length() != params_2.Length() ||
	    params_1.MinorDim() != NumberOfSupportParameters() ||
	    params_2.MinorDim() != NumberOfSupportParameters())
	{
		cout << "\n CubicSplineWindowGPT::SynchronizeSupportParameters: nodal\n"
		     << " parameters dimension mismatch" << endl;
		throw ExceptionT::kSizeMismatch;
	}
		
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

/* modify nodal shape function parameters */
void CubicSplineWindowGPT::ModifySupportParameters(dArray2DT& nodal_params) const
{
	/* scale supports */
	nodal_params *= fDilationScaling;	
}

void CubicSplineWindowGPT::WriteParameters(ostream& out) const
{
	/* window function parameters */
	out << " Dilation scaling factor . . . . . . . . . . . . = " << fDilationScaling << '\n';
}

/* Single point evaluations */
bool CubicSplineWindowGPT::Window(const dArrayT& x_n, const dArrayT& param_n, const dArrayT& x,
		int order, double& w, dArrayT& Dw, dSymMatrixT& DDw, dMatrixT& DDDw) // added KY (3rd derivative)
{
	/* allocate */
	int nsd = x.Length();
	/* outside the range of influence */
	if (!CubicSplineWindowGPT::Covers(x_n, x, param_n))
	{
		w = 0.0;
		if (order > 0)
		{
			Dw = 0.0;
			if (order > 1)
			{
				DDw = 0.0;
				if (order > 2) // kyonten
				{
					DDDw = 0.0;
				}	
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
		double a = param_n[0];
		double r = dist/a;
		if (r > 2.0) 
		{
			w = 0.0;
			Dw = 0.0;
			DDw = 0.0;
			DDDw = 0.0; // kyonten
		} 
		else if (r > 1.0) 
		{
			double dr = 2.0 - r;
			w = 1.0/(6.0*a)*dr*dr*dr;
			if (order > 0) 
			{
				double dw = -1.0/(2.0*a)*dr*dr;
				if (order > 1) 
				{
					double ddw = dr/a;
					DDw.Outer(Dw, (ddw/a - dw/dist)/(dist*dist*a));
					DDw.PlusIdentity(dw/(dist*a));
					if (order > 2) // kyonten (DDDw)
	  				{ 
	  					dSymMatrixT DDDw1(nsd);
	  					dMatrixT DDDw2(nsd,dSymMatrixT::NumValues(nsd));
	  					dArrayT DDDw1_vec(nsd), I(nsd);
	  					double const1, const2;
	  					// In 3D case: DDDw is a 3x9 matrix (non-symmetric) or a 3x6 (symmetric)
	  					// out of 27 (non-symmetric) or 18 (symmetric) components only 9
	  					// of them are needed for forming B3
	  					// DDDw, thus, becomes a 3x3 unsymmetric matrix
	  					DDDw1.Outer(Dw);
	  					if (nsd == 2)
	  					{
	  						DDDw1_vec[0] = DDDw1(0,0);
	  						DDDw1_vec[1] = DDDw1(1,1);
	  						I[0] = 1.0; I[1] = 1.0;
	  					}
	  					else if (nsd == 3)
	  					{
	  						DDDw1_vec[0]=DDDw1(0,0);
	  						DDDw1_vec[1]=DDDw1(1,1);
	  						DDDw1_vec[2]=DDDw1(2,2);
	  						I[0] = 1.0; I[1] = 1.0; I[2] = 1.0;
	  					}
	  					DDDw.Outer(Dw,DDDw1_vec);  
	  					const1 = (ddw/a - dw/dist);
	  					const1 *= -2./(dist*dist);
	  					const1 -= 1./(a*a*a*dist);
	  					const1 -= dr/(a*a*dist*dist);
	  					const1 += dw/(dist*dist*dist);
	  					DDDw *= const1/(dist*dist*a);
	  					const2 = (ddw/a - dw/dist)/(dist*dist*a);
	  					DDDw2.Outer(Dw,I); 
	  					DDDw2 += DDDw2;
	  					DDDw2 += DDDw2;
	  					DDDw2 *= const2;
	  					DDDw += DDDw2;
	  				}
				}
				Dw *= dw/(dist*a);
			}
		} 
		else 
		{
			w = (2.0/3.0 + r*r*(r/2.0 - 1.0))/a;
			if (order > 0) 
			{
				//double dw = r*(3.0*r/2.0 - 2.0)/a; re-write to avoid division by zero when r = 0.
				double dw_by_r = (3.0*r/2.0 - 2.0)/a;
				if (order > 1) 
				{
					double ddw = (3.0*r - 2.0)/a;
					DDw.Outer(Dw, (ddw/a - dw_by_r*r/dist)/(dist*dist*a));
					DDw.PlusIdentity(dw_by_r*r/(dist*a));
					if (order > 2) // kyonten (DDDw)
	  				{
	  					dSymMatrixT DDDw1(nsd);
	  					dMatrixT DDDw2(nsd,dSymMatrixT::NumValues(nsd));
	  					dArrayT DDDw1_vec(2), I(2);
	  					double const1, const2;
	  					
	  					DDDw1.Outer(Dw);
	  					if (nsd == 2)
	  					{
	  						DDDw1_vec[0] = DDDw1(0,0);
	  						DDDw1_vec[1] = DDDw1(1,1);
	  						I[0] = 1.0; I[1] = 1.0;
	  					}
	  					else if (nsd == 3)
	  					{
	  						DDDw1_vec[0]=DDDw1(0,0);
	  						DDDw1_vec[1]=DDDw1(1,1);
	  						DDDw1_vec[2]=DDDw1(2,2);
	  						I[0] = 1.0; I[1] = 1.0; I[2] = 1.0;
	  					}
	  					DDDw.Outer(Dw,DDDw1_vec); 
	  					const1 = (ddw/a - dw_by_r*r/dist);
	  					const1 *= -2./(dist*dist);
	  					const1 += 3./(a*a*a*dist);
	  					const1 -= (3.*r-1.)/(a*a*dist*dist);
	  					const1 += dw_by_r * r/(dist*dist*dist);
	  					DDDw *= const1/(dist*dist*a);
	  					const2 = (ddw/a - dw_by_r*r/dist)/(dist*dist*a);
	  					DDDw2.Outer(Dw,I); 
	  					DDDw2 += DDDw2;
	  					DDDw2 += DDDw2;
	  					DDDw2 *= const2;
	  					DDDw += DDDw2;
	  				}
				}
				Dw *= dw_by_r/(a*a);
			}
		}
	
		/* does cover */
		return true;	
	}
}

/* multiple point calculations */
int CubicSplineWindowGPT::Window(const dArray2DT& x_n, const dArray2DT& param_n, const dArrayT& x,
		int order, dArrayT& w, dArray2DT& Dw, dArray2DT& DDw, dArray2DT& DDDw) // kyonten (DDDw)
{
	/* allocate */
	int nsd = x.Length();
	fNSD.Dimension(nsd);
	fNSDsym.Dimension(nsd);
	fNSDunsym(nsd,nsd); // for DDDw??

	/* work space */
	dArrayT x_node, param_node;
	int count = 0;    
	int numwindowpoints = x_n.MajorDim(); 

	for (int i = 0; i < numwindowpoints; i++)
	{
		/* collect nodal values */
		x_n.RowAlias(i, x_node);
		param_n.RowAlias(i, param_node);
      
		if (CubicSplineWindowGPT::Window(x_node, param_node, x, order, w[i], fNSD, fNSDsym, fNSDunsym))
			count++;

		/* store derivatives */
		if (order > 0)
		{
			Dw.SetColumn(i, fNSD);
			if (order > 1)
			{
				DDw.SetColumn(i, fNSDsym);
				if (order > 2)
				{
					DDDw.SetColumn(i, fNSDunsym); // for DDDw??
				}
			}
		}
	}
	
	return count;
}

bool CubicSplineWindowGPT::Covers(const dArrayT& x_n, const dArrayT& x, const dArrayT& param_n) const
{
	double dist = dArrayT::Distance(x, x_n);
	return (dist/param_n[0] < 2.0);
}

int CubicSplineWindowGPT::Covers(const dArray2DT& x_n, const dArrayT& x, 
	const dArray2DT& param_n, ArrayT<bool>& covers) const
{
	int count = 0;
	int numwindowpoints = x_n.MajorDim();
	for (int i = 0; i < numwindowpoints; i++)
	{
		double dist = dArrayT::Distance(x, x_n);
		if (dist/param_n[0] < 2.0) {
			count++;
			covers[i] = true;
		} 
		else
			covers[i] = false;
	}
	
	return count;
}
