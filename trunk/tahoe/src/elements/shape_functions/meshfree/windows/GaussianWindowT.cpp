/* $Id: GaussianWindowT.cpp,v 1.14 2004-11-03 01:21:07 raregue Exp $ */
#include "GaussianWindowT.h"
#include "ExceptionT.h"
#include <math.h>

using namespace Tahoe;

const double sqrtPi = sqrt(acos(-1.0));
static double Max(double a, double b) { return (a > b) ? a : b; };

/* constructor */
GaussianWindowT::GaussianWindowT(double dilation_scaling, double sharpening_factor,
	double cut_off_factor):
	fDilationScaling(dilation_scaling),
	fSharpeningFactor(sharpening_factor),
	fCutOffFactor(cut_off_factor)
{
	if (fDilationScaling < 0.0 || fSharpeningFactor < 0.0 || fCutOffFactor < 1.0)
		ExceptionT::BadInputValue("GaussianWindowT::GaussianWindowT");
}

/* "synchronization" of nodal field parameters. */
void GaussianWindowT::SynchronizeSupportParameters(dArray2DT& params_1, 
	dArray2DT& params_2) const
{
	/* should be over the same global node set (marked by length) */
	if (params_1.Length() != params_2.Length() ||
	    params_1.MinorDim() != NumberOfSupportParameters() ||
	    params_2.MinorDim() != NumberOfSupportParameters())
	ExceptionT::SizeMismatch("GaussianWindowT::SynchronizeSupportParameters"); 
		
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

void GaussianWindowT::WriteParameters(ostream& out) const
{
	/* window function parameters */
	out << " Dilation scaling factor . . . . . . . . . . . . = " << fDilationScaling << '\n';
	out << " Window function sharpening factor . . . . . . . = " << fSharpeningFactor << '\n';
	out << " Neighbor cutoff factor. . . . . . . . . . . . . = " << fCutOffFactor << '\n';
}

/* Single point evaluations */
bool GaussianWindowT::Window(const dArrayT& x_n, const dArrayT& param_n, const dArrayT& x,
		int order, double& w, dArrayT& Dw, dSymMatrixT& DDw, dMatrixT& DDDw) //kyonten
{
	/* check out of influence range */
	if (!GaussianWindowT::Covers(x_n, x, param_n))
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
    	
    	/* no cover */
    	return false;
  	}
  	else
  	{
		/* distance */
		Dw.DiffOf(x_n, x);
		double dist = Dw.Magnitude();
		
  		/* scalar factors */
    	double adm = param_n[0]*fDilationScaling*fSharpeningFactor;
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
	  			if (order > 2) // kyonten
	  			{
					int nsd = x.Length();
					/* work space */
	  				dSymMatrixT DDDw1(nsd);
	  				dMatrixT DDDw2(nsd,nsd);
	  				dMatrixT DDDw3(nsd,nsd);
	  				dArrayT DDDw1_vec(nsd), I(nsd);
	  				/*  DDDw is a [nsd]x[nsd]x[nsd] or [nsd]x[nsd*nsd] matrix. 
	  				    using symmetry it reduces to [nsd]x[nstr]
	  					only the first three (3D) or two (2D) columns (contribution from 
	  					diagonal terms) of the [nsd]x[nstr] matrix are needed for the
	  					calculation of the Laplacian of strain tensor
	  					DDDw, thus, becomes a [nsd]x[nsd] unsymmetric matrix 
	  				*/
	  				DDDw1.Outer(Dw);
	  				for (int i = 0; i < nsd; i++)
	  				{
	  					for (int j = 0; i < nsd; j++)
	  					{
	  						if (i == j)
	  						{
	  							DDDw1_vec[j] = DDDw1(i,j); //collect diagonal terms
	  							I[j] = 1.0;
	  						}
	  					}
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
	  			} // (order > 2)
      		}
      		
      		/* set first derivative */
			Dw *= -2.0 * w / adm2;
    	}
    	
    	/* covers */
    	return true;
  	}
}

/* multiple point calculations */
int GaussianWindowT::Window(const dArray2DT& x_n, const dArray2DT& param_n, 
	const dArrayT& x, int order, dArrayT& w, dArray2DT& Dw, dArray2DT& DDw, dArray2DT& DDDw) //kyonten
{
	/* allocate */
	int nsd = x.Length();
	fNSD.Dimension(nsd);
	fNSDsym.Dimension(nsd);
	fNSDunsym.Dimension(nsd,nsd);
	
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
		if (GaussianWindowT::Window(x_node, param_node, x, order, w[i], fNSD, fNSDsym, fNSDunsym))
			count++;
			
		/* store derivatives */
		if (order > 0)
		{
			Dw.SetColumn(i, fNSD);
			if (order > 1)
			{
				DDw.SetColumn(i, fNSDsym);
				if (order > 2) //kyonten
					DDDw.SetColumn(i, fNSDunsym);
			}
		}
	}
	return count;
}

bool GaussianWindowT::Covers(const dArrayT& x_n, const dArrayT& x, const dArrayT& param_n) const
{
	double dist = dArrayT::Distance(x_n, x);
	if (dist > fCutOffFactor*fDilationScaling*param_n[0])
		return false;
	else
		return true;
}

int GaussianWindowT::Covers(const dArray2DT& x_n, const dArrayT& x, 
	const dArray2DT& param_n, ArrayT<bool>& covers) const
{
	int count = 0; /* # of point covered... */
	int numwindows = x_n.MajorDim();
	for (int i = 0; i < numwindows; i++)
  	{
		double dist = dArrayT::Distance(x, x_n);

		/* check out of influence range */
		if (dist > fCutOffFactor*fDilationScaling*param_n[i])
			covers[i] = false;
		else
		{
			count++;
			covers[i] = true;
		}
  	}
  	return count;
}

/* spherical upport size */
double GaussianWindowT::SphericalSupportSize(const dArrayT& param_n) const
{
#if __option(extended_errorcheck)
	if (param_n.Length() != 1) ExceptionT::GeneralFail("GaussianWindowT::SphericalSupportSize");
#endif
	return fCutOffFactor*fDilationScaling*param_n[0];
}

/* rectangular support size */
void GaussianWindowT::RectangularSupportSize(const dArrayT& param_n, dArrayT& support_size) const
{
	/* same in all dimensions */
	support_size = SphericalSupportSize(param_n);
}

/* spherical support sizes in batch */
void GaussianWindowT::SphericalSupportSize(const dArray2DT& param_n, ArrayT<double>& support_size) const
{
#if __option(extended_errorcheck)
	if (param_n.MinorDim() != 1 ||
	    param_n.MajorDim() != support_size.Length()) 
		ExceptionT::SizeMismatch("GaussianWindowT::SphericalSupportSize");
#endif

	dArrayT tmp;
	tmp.Alias(support_size);
	tmp.SetToScaled(fCutOffFactor*fDilationScaling, param_n);
}
