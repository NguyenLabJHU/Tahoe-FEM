/* $Id: GaussianWindowGPT.cpp,v 1.3 2004-07-14 19:50:31 kyonten Exp $ */

#include "GaussianWindowGPT.h"
#include "ExceptionT.h"
#include <math.h>


using namespace Tahoe;

const double sqrtPi = sqrt(acos(-1.0));
static double Max(double a, double b) { return (a > b) ? a : b; };

/* constructor */
GaussianWindowGPT::GaussianWindowGPT(double dilation_scaling, double sharpening_factor,
	double cut_off_factor):
	fDilationScaling(dilation_scaling),
	fSharpeningFactor(sharpening_factor),
	fCutOffFactor(cut_off_factor)
{
	if (fDilationScaling < 0.0 || fSharpeningFactor < 0.0 || fCutOffFactor < 1.0)
		throw ExceptionT::kBadInputValue;
}

/* "synchronization" of nodal field parameters. */
void GaussianWindowGPT::SynchronizeSupportParameters(dArray2DT& params_1, 
	dArray2DT& params_2) const
{
	/* should be over the same global node set (marked by length) */
	if (params_1.Length() != params_2.Length() ||
	    params_1.MinorDim() != NumberOfSupportParameters() ||
	    params_2.MinorDim() != NumberOfSupportParameters())
	{
		cout << "\n GaussianWindowGPT::SynchronizeSupportParameters: nodal\n"
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
void GaussianWindowGPT::ModifySupportParameters(dArray2DT& nodal_params) const
{
	/* scale supports */
	nodal_params *= fDilationScaling;
}

void GaussianWindowGPT::WriteParameters(ostream& out) const
{
	/* window function parameters */
	out << " Dilation scaling factor . . . . . . . . . . . . = " << fDilationScaling << '\n';
	out << " Window function sharpening factor . . . . . . . = " << fSharpeningFactor << '\n';
	out << " Neighbor cutoff factor. . . . . . . . . . . . . = " << fCutOffFactor << '\n';
}

/* Single point evaluations */
bool GaussianWindowGPT::Window(const dArrayT& x_n, const dArrayT& param_n, const dArrayT& x,
		int order, double& w, dArrayT& Dw, dSymMatrixT& DDw, dMatrixT& DDDw)
{
	/* allocate */  //kyonten
	int nsd = x.Length();
	/* check out of influence range */
	if (!Covers(x_n, x, param_n))
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
    	
    	/* no cover */
    	return false;
  	}
  	else
  	{
		/* distance */
		Dw.DiffOf(x_n, x);
		double dist = Dw.Magnitude();
		
  		/* scalar factors */
    	double adm = param_n[0] * fSharpeningFactor;
    	double adm2 = adm * adm;
    	double q = dist / adm;
    	w = exp(-q * q) / (sqrtPi * adm);
    	if (order > 0)
    	{
    		/* compute second derivative */
      		if (order > 1)
      		{
	  			DDw.Outer(Dw);
	  			DDw *= 4.0 * w / (adm2 * adm2);
	  			DDw.PlusIdentity(-2.0 * w / adm2);
	  			if (order > 2) // kyonten (DDDw)
	  			{
	  				dSymMatrixT DDDw1(nsd);
	  				dMatrixT DDDw2(nsd,dSymMatrixT::NumValues(nsd));
	  				dMatrixT DDDw3(nsd,dSymMatrixT::NumValues(nsd));
	  				dArrayT DDDw1_vec(nsd), I(nsd);
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
	  				DDDw *= -8.0*w/(adm2*adm2*adm2);
	  				DDDw2.Outer(Dw,I);
	  				DDDw2 += DDDw2;
	  				DDDw2 *= 4.0*w/(adm2*adm2);
	  				DDDw += DDDw2;
	  				DDDw3.Outer(Dw,I); 
	  				DDDw3 *= 4.0*w/(adm2*adm2);
	  				DDDw += DDDw3;
	  			}
      		}
      		
      		/* set first derivative */
			Dw *= -2.0 * w / adm2;
    	}
    	
    	/* covers */
    	return true;
  	}
}

/* multiple point calculations */
int GaussianWindowGPT::Window(const dArray2DT& x_n, const dArray2DT& param_n, 
	const dArrayT& x, int order, dArrayT& w, dArray2DT& Dw, dArray2DT& DDw, dArray2DT& DDDw)
{
	/* allocate */
	int nsd = x.Length();
	fNSD.Dimension(nsd);
	fNSDsym.Dimension(nsd);
	fNSDunsym(nsd,nsd); 
	
	/* work space */
	dArrayT x_node, param_node;

	int count = 0;
	int num_points = x_n.MajorDim();
	for (int i = 0; i < num_points; i++)
	{
		/* collect nodal values */
		x_n.RowAlias(i, x_node);
		param_n.RowAlias(i, param_node);
	
		/* single point evaluation (override virtual) */
		if (GaussianWindowGPT::Window(x_node, param_node, x, order, w[i], fNSD, fNSDsym, fNSDunsym))
			count++;
			
		/* store derivatives */
		if (order > 0)
		{
			Dw.SetColumn(i, fNSD);
			if (order > 1)
			{
				DDw.SetColumn(i, fNSDsym);
				if (order > 2) //kyonten
				{
					DDDw.SetColumn(i, fNSDunsym); 
				}	
			}
		}
	}
	return count;
}

bool GaussianWindowGPT::Covers(const dArrayT& x_n, const dArrayT& x, const dArrayT& param_n) const
{
	double dist = dArrayT::Distance(x_n, x);
	if (dist > fCutOffFactor*param_n[0])
		return false;
	else
		return true;
}

int GaussianWindowGPT::Covers(const dArray2DT& x_n, const dArrayT& x, 
	const dArray2DT& param_n, ArrayT<bool>& covers) const
{
	int count = 0;    // # of point covered...
	int numwindows = x_n.MajorDim();
	for (int i = 0; i < numwindows; i++)
  	{
		double dist = dArrayT::Distance(x, x_n);

		/* check out of influence range */
		if (dist > fCutOffFactor*param_n[i])
			covers[i] = false;
		else
		{
			count++;
			covers[i] = true;
		}
  	}
  	return count;
}
