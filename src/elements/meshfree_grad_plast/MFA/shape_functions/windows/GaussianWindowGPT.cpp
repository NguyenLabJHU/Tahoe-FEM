/* $Id: GaussianWindowGPT.cpp,v 1.2 2004-06-30 16:56:21 kyonten Exp $ */

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
	  			    if (nsd == 2)
	  			    {
	  			     len = 3.; // upper triangle only (symmetry)
	  			    }
	  			    else if (nsd == 3)
	  			    {
	  			     len = 6.; // upper triangle only (symmetry)
	  			    }
	  				dMatrixT AA(1, len), II(1, len); 
	  				dMatrixT DDDw1(nsd, nsd), I(nsd, nsd);
	  				I = 0.; II = 0.;
	  				DDDw1.Outer(Dw);
	  				if (nsd == 2)
	  				{
	  				 AA(1,1) = DDDw1(1,1); AA(1,2) = DDDw1(2,2); AA(1,3) = DDDw1(1,2);
	  				 I(1,1) = I(2,2) = 1.;
	  				 II(1,1) = I(1,1); II(1,2) = I(2,2); II(1,3) = I(1,2);
	  				}
	  				else if (nsd == 3)
	  				{
	  				 AA(1,1) = DDDw1(1,1); AA(1,2) = DDDw1(2,2); AA(1,3) = DDDw1(3,3);
	  				 AA(1,4) = DDDw1(2,3); AA(1,5) = DDDw1(1,3): AA(1,6) = DDDw1(1,2);
	  				 I(1,1) = I(2,2) = I(3,3) = 1.;
	  				 II(1,1) = I(1,1); II(1,2) = I(2,2); II(1,3) = I(3,3);
	  				 II(1,4) = I(2,3); II(1,5) = I(1,3): II(1,6) = II(1,2);
	  				}
	  				DDDw.MultAB(Dw,DDDw1); // 3x6
	  				DDDw *= -8.0*w/(adm2*adm2*adm2);
	  				DDDw2.MultAB(Dw,II); // 3x6
	  				DDDw2 += DDDw2;
	  				DDDw2 *= 4.0*w/(adm2*adm2);
	  				DDDw += DDDw2;
	  				DDDw3.MultAB(Dw,II); // 3x6
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
		if (GaussianWindowGPT::Window(x_node, param_node, x, order, w[i], fNSD, fNSDsym))
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
					DDDw.SetColumn(i, fNSDsym); 
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
