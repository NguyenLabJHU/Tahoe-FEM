/* $Id: RectGaussianWindowGPT.cpp,v 1.2 2004-07-14 19:50:31 kyonten Exp $ */

#include "RectGaussianWindowGPT.h"
#include "ExceptionT.h"
#include <math.h>


using namespace Tahoe;

const double sqrtPi = sqrt(acos(-1.0));
static double Max(double a, double b) { return (a > b) ? a : b; };

/* constructor */
RectGaussianWindowGPT::RectGaussianWindowGPT(const dArrayT& dilation_scaling, double sharpening_factor,
						       double cut_off_factor):
	fDilationScaling(dilation_scaling),
	fSharpeningFactor(sharpening_factor),
	fCutOffFactor(cut_off_factor)
{
        int count = 0;
        for (int i = 0; i < fDilationScaling.Length(); i++)
	{
	  if (fDilationScaling[i] < 0.0)
	    count++;
	}
	if (count > 0 || fSharpeningFactor < 0.0)
		throw ExceptionT::kBadInputValue;
}

/* "synchronization" of nodal field parameters. */
void RectGaussianWindowGPT::SynchronizeSupportParameters(dArray2DT& params_1, 
	dArray2DT& params_2) const
{
	/* should be over the same global node set (marked by length) */
	if (params_1.Length() != params_2.Length() ||
	    params_1.MinorDim() != NumberOfSupportParameters() ||
	    params_2.MinorDim() != NumberOfSupportParameters())
	{
		cout << "\n RectGaussianWindowGPT::SynchronizeSupportParameters: nodal\n"
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
void RectGaussianWindowGPT::ModifySupportParameters(dArray2DT& nodal_params) const
{
	/* scale supports */
	for (int i = 0; i < fDilationScaling.Length(); i++)
		nodal_params.ScaleColumn(i, fDilationScaling[i]);	
}

void RectGaussianWindowGPT::WriteParameters(ostream& out) const
{
	/* window function parameters */
	out << " Dilation scaling factor . . . . . . . . . . . . :\n" << fDilationScaling << '\n';
	out << " Window function sharpening factor . . . . . . . = " << fSharpeningFactor << '\n';
	out << " Neighbor cutoff factor. . . . . . . . . . . . . = " << fCutOffFactor << '\n';
}

/* Single point evaluations */
bool RectGaussianWindowGPT::Window(const dArrayT& x_n, const dArrayT& param_n, const dArrayT& x,
		int order, double& w, dArrayT& Dw, dSymMatrixT& DDw, dMatrixT& DDDw) // kyonten (DDDw)
{
  /* Compute window function and its derivatives - accomplish by scalar product of individual
   * window functions in x/y/z directions */

  if (!RectGaussianWindowGPT::Covers(x_n, x, param_n))
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
    return false;
  }
  else      
  {
    Dw.DiffOf(x, x_n);
    if (x.Length() == 2)      // 2D calculation
    {
      double admx = param_n[0] * fSharpeningFactor;
      double admy = param_n[1] * fSharpeningFactor;
      double admx2 = admx * admx;
      double admy2 = admy * admy;
      double qx = Dw[0] / admx;
      double qy = Dw[1] / admy;
      double wx = exp(-qx * qx) / (sqrtPi * admx);
      double wy = exp(-qy * qy) / (sqrtPi * admy);
      w = wx * wy;
      if (order > 0)
      {
	if (order > 1)
	{
	  if (order > 2) // kyonten (DDDw)
	  {
	  	DDDw[0] = (-2.0 * w / admx2) * (2.0 * Dw[0] * Dw[0] / admx2 - 1)
	  			*(2.0 * DDw[0] * Dw[0] / admx2)*(2.0 * Dw[0] * DDw[0] / admx2); //double check!!
	  	DDDw[1] = (-2.0 * w / admy2) * (2.0 * Dw[1] * Dw[1] / admy2 - 1)
	  			*(2.0 * DDw[1] * Dw[1] / admy2)*(2.0 * Dw[1] * DDw[1] / admy2); 
	  }
	  DDw[0] = (2.0 * w / admx2) * (2.0 * Dw[0] * Dw[0] / admx2 - 1);
	  DDw[1] = (2.0 * w / admy2) * (2.0 * Dw[1] * Dw[1] / admy2 - 1);
	}
	Dw[0] *= -2.0 * w / admx2;
	Dw[1] *= -2.0 * w / admy2;
      }

    }
    else if (x.Length() == 3)     // 3D calculation
    {
      double admx = param_n[0] * fSharpeningFactor;
      double admy = param_n[1] * fSharpeningFactor;
      double admz = param_n[2] * fSharpeningFactor;
      double admx2 = admx * admx;
      double admy2 = admy * admy;
      double admz2 = admz * admz;
      double qx = Dw[0] / admx;
      double qy = Dw[1] / admy;
      double qz = Dw[2] / admz;
      double wx = exp(-qx * qx) / (sqrtPi * admx);
      double wy = exp(-qy * qy) / (sqrtPi * admy);
      double wz = exp(-qz * qz) / (sqrtPi * admz);
      w = wx * wy * wz;
      if (order > 0)
      {
	if (order > 1)
	{
	  if (order > 2) // kyonten (DDDw)
	  {
	  	DDDw[0] = (2.0 * w / admx2) * (2.0 * Dw[0] * Dw[0] / admx2 - 1)
	  			*(2.0 * DDw[0] * Dw[0] / admx2)*(2.0 * Dw[0] * DDw[0] / admx2); //double check!!
	  	DDDw[1] = (2.0 * w / admy2) * (2.0 * Dw[1] * Dw[1] / admy2 - 1)
	  			*(2.0 * DDw[1] * Dw[1] / admy2)*(2.0 * Dw[1] * DDw[1] / admy2);
	  	DDDw[2] = (2.0 * w / admz2) * (2.0 * Dw[2] * Dw[2] / admz2 - 1)
	  			*(2.0 * DDw[2] * Dw[2] / admz2)*(2.0 * Dw[2] * DDw[2] / admz2); 
	  }
	  DDw[0] = (2.0 * w / admx2) * (2.0 * Dw[0] * Dw[0] / admx2 - 1);
	  DDw[1] = (2.0 * w / admy2) * (2.0 * Dw[1] * Dw[1] / admy2 - 1);
	  DDw[2] = (2.0 * w / admz2) * (2.0 * Dw[2] * Dw[2] / admz2 - 1);
	}
	Dw[0] *= -2.0 * w / admx2;
	Dw[1] *= -2.0 * w / admy2;
	Dw[2] *= -2.0 * w / admz2;
      }
    }
    else
      throw ExceptionT::kGeneralFail;
      
		return true;
  }
  return false;
}

/* multiple point calculations */
int RectGaussianWindowGPT::Window(const dArray2DT& x_n, const dArray2DT& param_n, const dArrayT& x,
		int order, dArrayT& w, dArray2DT& Dw, dArray2DT& DDw, dArray2DT& DDDw)// kyonten (DDDw)
{
  /* compute window function and derivatives for multiple field points */

  /* allocate */
  int nsd = x.Length();
  fNSD.Dimension(nsd);
  fNSDsym.Dimension(nsd);
  fNSDunsym(nsd,nsd); // allocation for DDDw??

  /* work space */
  dArrayT x_node, param_node;
  int count = 0;    
  int numwindowpoints = x_n.MajorDim(); 

  for (int i = 0; i < numwindowpoints; i++)
  {
    /* collect nodal values */
    x_n.RowAlias(i, x_node);
    param_n.RowAlias(i, param_node);
      
    if (RectGaussianWindowGPT::Window(x_node, param_node, x, order, w[i], fNSD, fNSDsym, fNSDunsym)) //last fNSDsym??
      count ++;

    /* store derivatives */
    if (order > 0)
    {
      Dw.SetColumn(i, fNSD);
      if (order > 1)
      {
      	DDw.SetColumn(i, fNSDsym);
      	if (order > 2)
	  	{
	  		DDDw.SetColumn(i, fNSDunsym);	// dimension correct??
	  	}	
      }
    }
  }
  return count;
}

bool RectGaussianWindowGPT::Covers(const dArrayT& x_n, const dArrayT& x, 
	const dArrayT& param_n) const
{
  dArrayT dx(x.Length());
  dx.DiffOf(x, x_n);

  /* check individual directions to see if outside the "box" */
  for (int i = 0; i < x.Length(); i++)
  {
    if (fabs(dx[i]) > fCutOffFactor * param_n[i])
      return false;
  }
  return true;
}

int RectGaussianWindowGPT::Covers(const dArray2DT& x_n, const dArrayT& x, 
	const dArray2DT& param_n, ArrayT<bool>& covers) const
{
  int count = 0;
  int numwindowpoints = x_n.MajorDim();        // Could be MajorDim!
  dArrayT dx(x.Length()), temprow(x.Length());
  for (int i = 0; i < numwindowpoints; i++)
  {
    x_n.RowCopy(i, temprow);          // Could be COLUMN copy!!!
    dx.DiffOf(x, temprow);
    for (int j = 0; j < x.Length(); j++)
    {
      if (fabs(dx[j]) > fCutOffFactor * param_n(i,j))
		count++;
    }
    if (count == 0)
      covers[i] = true;
    else
      covers[i] = false;
  }
  return count;
}
