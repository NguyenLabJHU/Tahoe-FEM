/* $Id: RectCubicSplineWindowT.cpp,v 1.12 2004-10-31 20:48:41 paklein Exp $ */
#include "RectCubicSplineWindowT.h"
#include "ExceptionT.h"
#include <math.h>

using namespace Tahoe;

const double sqrtPi = sqrt(acos(-1.0));
static double Max(double a, double b) { return (a > b) ? a : b; };

/* constructor */
RectCubicSplineWindowT::RectCubicSplineWindowT(const dArrayT& dilation_scaling):
	fDilationScaling(dilation_scaling)
{
	double min_scaling = dilation_scaling.Min();
	if (min_scaling < kSmall)
		ExceptionT::GeneralFail("RectCubicSplineWindowT::RectCubicSplineWindowT");
}

/* "synchronization" of nodal field parameters. */
void RectCubicSplineWindowT::SynchronizeSupportParameters(dArray2DT& params_1, 
	dArray2DT& params_2) const
{
	/* should be over the same global node set (marked by length) */
	if (params_1.Length() != params_2.Length() ||
	    params_1.MinorDim() != NumberOfSupportParameters() ||
	    params_2.MinorDim() != NumberOfSupportParameters())
		ExceptionT::SizeMismatch("RectCubicSplineWindowT::SynchronizeSupportParameters");

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

void RectCubicSplineWindowT::WriteParameters(ostream& out) const
{
	/* window function parameters */
	out << " Dilation scaling factor . . . . . . . . . . . . :\n" << fDilationScaling << '\n';
}

/* Single point evaluations */
bool RectCubicSplineWindowT::Window(const dArrayT& x_n, const dArrayT& param_n, const dArrayT& x,
		int order, double& w, dArrayT& Dw, dSymMatrixT& DDw)
{
  /* Compute window function and its derivatives - accomplish by scalar product of individual
   * window functions in x/y/z directions */

  if (!RectCubicSplineWindowT::Covers(x_n, x, param_n))
  {
    w = 0.0;
    if (order > 0)
    {
      Dw = 0.0;
      if (order > 1)
	DDw = 0.0;
    }
    return false;
  }
  else      
  {
    Dw.DiffOf(x, x_n);
    if (x.Length() == 2)      // 2D calculation
    {
      double wx, wy;
      
      const double ax = 0.5*fDilationScaling[0]*param_n[0];
      const double ay = 0.5*fDilationScaling[1]*param_n[1];
      
      const double rx = Dw[0]/ax;
      const double ry = Dw[1]/ay;

      /* check x-direction */
      if ((rx>=-2)&&(rx<-1)) 
	wx = 1.0 / (6.0 * ax) * (2.0 + rx) * (2.0 + rx) * (2.0 + rx);
      else if ((rx>=-1)&&(rx<0)) 
	wx = (2.0 / 3.0 - rx * rx - .5 * rx * rx * rx) / ax;
      else if ((rx>=0)&&(rx<1)) 
	wx = (2.0 / 3.0 - rx * rx + .5 * rx * rx * rx) / ax;
      else if ((rx>=1)&&(rx<=2)) 
	wx = 1.0 / (6.0 * ax) * (2.0 - rx) * (2.0 - rx) * (2.0 - rx);
      else
	wx = 0.0;

      /* check y-direction */
      if ((ry>=-2)&&(ry<-1))
	wy = 1.0 / (6.0 * ay) * (ry + 2.0) * (ry + 2.0) * (ry + 2.0);
      else if ((ry>=-1)&&(ry<0))
	wy = (2.0 / 3.0 - ry * ry - .5 * ry * ry * ry) / ay;
      else if ((ry>=0)&&(ry<1))
	wy = (2.0 / 3.0 - ry * ry + .5 * ry * ry * ry) / ay;
      else if ((ry>=1)&&(ry<=2))
	wy = 1.0 / (6.0 * ay) * (2.0 - ry) * (2.0 - ry) * (2.0 - ry);
      else
	wy = 0.0;

      w = wx * wy;

      if (order > 0)
      {
	if (order > 1)
	{
	  /* check x-direction */
	  if ((rx >= -2) && (rx < -1))
	    DDw[0] = wy / (ax * ax * ax) * (2.0 + rx);
	  else if ((rx >= -1) && (rx < 0))
	    DDw[0] = - wy / (ax * ax * ax) * (2.0 + 3.0 * rx);
	  else if ((rx >= 0) && (rx < 1))
	    DDw[0] = - wy / (ax * ax * ax) * (2.0 - 3.0 * rx);
	  else if ((rx >= 1) && (rx <= 2))
	    DDw[0] = wy / (ax * ax * ax) * (2.0 - rx);
	  else
	    DDw[0] = 0.0;
	  
	  /* check y-direction */
	  if ((ry >= -2) && (ry < -1))
	    DDw[1] = wx / (ay * ay * ay) * (2.0 + ry);
	  else if ((ry >= -1) && (ry < 0))
	    DDw[1] = - wx / (ay * ay * ay) * (2.0 + 3.0 * ry);
	  else if ((ry >= 0) && (ry < 1))
	    DDw[1] = - wx / (ay * ay * ay) * (2.0 - 3.0 * ry);
	  else if ((ry >= 1) && (ry <= 2))
	    DDw[1] = wx / (ay * ay * ay) * (2.0 - ry);
	  else
	    DDw[1] = 0.0;
	}

	/* check x-direction */
	if ((rx >= -2) && (rx < -1))
	  Dw[0] = - .5 * wy / (ax * ax) * (2.0 + rx) * (2.0 + rx);
	else if ((rx >= -1) && (rx < 0))
	  Dw[0] = wy / (ax * ax) * (2.0 * rx + 1.5 * rx * rx);
	else if ((rx >= 0) && (rx < 1))
	  Dw[0] = wy / (ax * ax) * (2.0 * rx - 1.5 * rx * rx);
	else if ((rx >= 1) && (rx <= 2))
	  Dw[0] = .5 * wy / (ax * ax) * (2.0 - rx) * (2.0 - rx);
	else
	  Dw[0] = 0.0;
	  
	  /* check y-direction */
	if ((ry >= -2) && (ry < -1))
	  Dw[1] = - .5 * wx / (ay * ay) * (2.0 + ry) * (2.0 + ry);
	else if ((ry >= -1) && (ry < 0))
	  Dw[1] = wx / (ay * ay) * (2.0 * ry + 1.5 * ry * ry);
	else if ((ry >= 0) && (ry < 1))
	  Dw[1] = wx / (ay * ay) * (2.0 * ry - 1.5 * ry * ry);
	else if ((ry >= 1) && (ry <= 2))
	  Dw[1] = .5 * wx / (ay * ay) * (2.0 - ry) * (2.0 - ry);
	else
	  Dw[1] = 0.0;

      }
    }
    else if (x.Length() == 3)     // 3D calculation
    {
      double wx, wy, wz;

      const double ax = 0.5*fDilationScaling[0]*param_n[0];
      const double ay = 0.5*fDilationScaling[1]*param_n[1];
      const double az = 0.5*fDilationScaling[2]*param_n[2];
          
      double rx = Dw[0]/ax;
      double ry = Dw[1]/ay;
      double rz = Dw[2]/az;
      
      /* check x-direction */
      if ((rx >= -2) && (rx < -1))
	wx = 1.0 / (6.0 * ax) * (2.0 + rx) * (2.0 + rx) * (2.0 + rx);
      else if ((rx >= -1) && (rx < 0))
	wx = (2.0 / 3.0 - rx * rx - .5 * rx * rx * rx) / ax;
      else if ((rx >= 0) && (rx < 1))
	wx = (2.0 / 3.0 - rx * rx + .5 * rx * rx * rx) / ax;
      else if ((rx >= 1) && (rx <= 2))
	wx = 1.0 / (6.0 * ax) * (2 - rx) * (2 - rx) * (2 - rx);
      else
	wx = 0.0;

      /* check y-direction */
      if ((ry >= -2) && (ry < -1))
	wy = 1.0 / (6.0 * ay) * (ry + 2.0) * (ry + 2.0) * (ry + 2.0);
      else if ((ry >= -1) && (ry < 0))
	wy = 1.0 / ay * (2.0 / 3.0 - ry * ry * (1.0 + .5 * ry));
      else if ((ry >= 0) && (ry < 1))
	wy = 1.0 / ay * (2.0 / 3.0 - ry * ry * (1.0 - .5 * ry));
      else if ((ry >= 1) && (ry <= 2))
	wy = 1.0 / (6.0 * ay) * (2.0 - ry) * (2.0 - ry) * (2.0 - ry);
      else
	wy = 0.0;

      /* check z-direction */
      if ((rz >= -2) && (rz < -1))
	wz = 1.0 / (6.0 * az) * (rz + 2.0) * (rz + 2.0) * (rz + 2.0);
      else if ((rz >= -1) && (rz < 0))
	wz = 1.0 / az * (2.0 / 3.0 - rz * rz * (1.0 + .5 * rz));
      else if ((rz >= 0) && (rz < 1))
	wz = 1.0 / az * (2.0 / 3.0 - rz * rz * (1.0 - .5 * rz));
      else if ((rz >= 1) && (rz <= 2))
	wz = 1.0 / (6.0 * az) * (2.0 - rz) * (2.0 - rz) * (2.0 - rz);
      else
	wz = 0.0;

      w = wx * wy * wz;
      if (order > 0)
      {
	if (order > 1)
	{
	  /* check x-direction */
	  if ((rx >= -2) && (rx < -1))
	    DDw[0] = wy * wz / (ax * ax * ax) * (2.0 + rx);
	  else if ((rx >= -1) && (rx < 0))
	    DDw[0] = - wy * wz / (ax * ax * ax) * (2.0 + 3.0 * rx);
	  else if ((rx >= 0) && (rx < 1))
	    DDw[0] = - wy * wz / (ax * ax * ax) * (2.0 - 3.0 * rx);
	  else if ((rx >= 1) && (rx <= 2))
	    DDw[0] = wy * wz / (ax * ax * ax) * (2.0 - rx);
	  else
	    DDw[0] = 0.0;
	  
	  /* check y-direction */
	  if ((ry >= -2) && (ry < -1))
	    DDw[1] = wx * wz / (ay * ay * ay) * (2.0 + ry);
	  else if ((ry >= -1) && (ry < 0))
	    DDw[1] = - wx * wz / (ay * ay * ay) * (2.0 + 3.0 * ry);
	  else if ((ry >= 0) && (ry < 1))
	    DDw[1] = - wx * wz / (ay * ay * ay) * (2.0 - 3.0 * ry);
	  else if ((ry >= 1) && (ry <= 2))
	    DDw[1] = wx * wz / (ay * ay * ay) * (2.0 - ry);
	  else
	    DDw[1] = 0.0;

	  /* check z-direction */
	  if ((rz >= -2) && (rz < -1))
	    DDw[2] = wx * wy / (az * az * az) * (2.0 + rz);
	  else if ((rz >= -1) && (rz < 0))
	    DDw[2] = - wx * wy / (az * az * az) * (2.0 + 3.0 * rz);
	  else if ((rz >= 0) && (rz < 1))
	    DDw[2] = - wx * wy / (az * az * az) * (2.0 - 3.0 * rz);
	  else if ((rz >= 1) && (rz <= 2))
	    DDw[2] = wx * wy / (az * az * az) * (2.0 - rz);
	  else
	    DDw[2] = 0.0;
	}

	/* check x-direction */
	if ((rx >= -2) && (rx < -1))
	  Dw[0] = - .5 * wy * wz / (ax * ax) * (2.0 + rx) * (2.0 + rx);
	else if ((rx >= -1) && (rx < 0))
	  Dw[0] = wy * wz / (ax * ax) * (2.0 * rx + 1.5 * rx * rx);
	else if ((rx >= 0) && (rx < 1))
	  Dw[0] = wy * wz / (ax * ax) * (2.0 * rx - 1.5 * rx * rx);
	else if ((rx >= 1) && (rx <= 2))
	  Dw[0] = .5 * wy * wz / (ax * ax) * (2.0 - rx) * (2.0 - rx);
	else
	  Dw[0] = 0.0;
	
	/* check y-direction */
	if ((ry >= -2) && (ry < -1))
	  Dw[1] = - .5 * wx * wz / (ay * ay) * (2.0 + ry) * (2.0 + ry);
	else if ((ry >= -1) && (ry < 0))
	  Dw[1] = wx * wz / (ay * ay) * (2.0 * ry + 1.5 * ry * ry);
	else if ((ry >= 0) && (ry < 1))
	  Dw[1] = wx * wz / (ay * ay) * (2.0 * ry - 1.5 * ry * ry);
	else if ((ry >= 1) && (ry <= 2))
	  Dw[1] = .5 * wx * wz / (ay * ay) * (2.0 - ry) * (2.0 - ry);
	else
	  Dw[1] = 0.0;

	/* check z-direction */
	if ((rz >= -2) && (rz < -1))
	  Dw[2] = - .5 * wx * wy / (az * az) * (2.0 + rz) * (2.0 + rz);
	else if ((rz >= -1) && (rz < 0))
	  Dw[2] = wx * wy / (az * az) * (2.0 * rz + 1.5 * rz * rz);
	else if ((rz >= 0) && (rz < 1))
	  Dw[2] = wx * wy / (az * az) * (2.0 * rz - 1.5 * rz * rz);
	else if ((rz >= 1) && (rz <= 2))
	  Dw[2] = .5 * wx * wy / (az * az) * (2.0 - rz) * (2.0 - rz);
	else
	  Dw[2] = 0.0;
      }
    }
    else
		ExceptionT::GeneralFail("RectCubicSplineWindowT::Window");
      
    return true;
  }
  return false;
}

/* multiple point calculations */
int RectCubicSplineWindowT::Window(const dArray2DT& x_n, const dArray2DT& param_n, const dArrayT& x,
		int order, dArrayT& w, dArray2DT& Dw, dArray2DT& DDw)
{
  /* compute window function and derivatives for multiple field points */

  /* allocate */
  int nsd = x.Length();
  fNSD.Dimension(nsd);
  fNSDsym.Dimension(nsd);

  /* work space */
  dArrayT x_node, param_node;
  int count = 0;    
  int numwindowpoints = x_n.MajorDim(); 

  for (int i = 0; i < numwindowpoints; i++)
  {
    /* collect nodal values */
    x_n.RowAlias(i, x_node);
    param_n.RowAlias(i, param_node);
      
    if (RectCubicSplineWindowT::Window(x_node, param_node, x, order, w[i], fNSD, fNSDsym))
      count ++;

    /* store derivatives */
    if (order > 0)
    {
      Dw.SetColumn(i, fNSD);
      if (order > 1)
	DDw.SetColumn(i, fNSDsym);
    }
  }
  return count;
}

bool RectCubicSplineWindowT::Covers(const dArrayT& x_n, const dArrayT& x, 
	const dArrayT& param_n) const
{
	for (int i = 0; i < x.Length(); i++)
		if (fabs(x_n[i] - x[i]) > fDilationScaling[i]*param_n[i])
			return false;

	return true;
}

int RectCubicSplineWindowT::Covers(const dArray2DT& x_n, const dArrayT& x, 
	const dArray2DT& param_n, ArrayT<bool>& covers) const
{
	/* initialize */
	covers = true;
	int count = x_n.MajorDim();
	
	/* search */
	int nsd = x_n.MinorDim();
	int npr = param_n.MinorDim();
	for (int i = 0; i < nsd; i++) /* spatial dimensions */
	{
		double x_i = x[i];
		double scale = fDilationScaling[i];
		const double* px_n = x_n.Pointer(i);
		const double* pparam_n = param_n.Pointer(i);
	
		for (int j = 0; j < x_n.MajorDim(); j++) /* test points */
		{
			/* check */
			if (fabs(px_n[i] - x_i) > scale*pparam_n[j] && covers[j])
			{
				covers[j] = false;
				count--;
			}
			
			/* next */
			px_n += nsd;
			pparam_n += npr;
		}
	}		

	return count;
}

/* spherical upport size */
double RectCubicSplineWindowT::SphericalSupportSize(const dArrayT& param_n) const
{
	dArrayT support_size(param_n.Length());
	RectangularSupportSize(param_n, support_size);
	return support_size.Magnitude();
}

/* rectangular support size */
void RectCubicSplineWindowT::RectangularSupportSize(const dArrayT& param_n, dArrayT& support_size) const
{
	for (int i = 0; i < support_size.Length(); i++)
		support_size[i] = param_n[i]*fDilationScaling[i];
}
