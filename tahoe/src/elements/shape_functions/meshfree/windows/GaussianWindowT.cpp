/* $Id: GaussianWindowT.cpp,v 1.2 2001-06-14 17:18:35 hspark Exp $ */

#include "GaussianWindowT.h"
#include "ExceptionCodes.h"
#include <math.h>

const double sqrtPi = sqrt(acos(-1.0));

/* constructor */
GaussianWindowT::GaussianWindowT(double dilation_scaling, double sharpening_factor):
	fDilationScaling(dilation_scaling),
	fSharpeningFudgeFactor(sharpening_factor)
{
	if (fDilationScaling < 0.0 || fSharpeningFudgeFactor < 0.0)
		throw eBadInputValue;
}

void GaussianWindowT::WriteParameters(ostream& out) const
{
  /* Not sure what to do here */


}

/* Single point evaluations */
void GaussianWindowT::window(const dArrayT& x_0, const dArrayT& param, const dArrayT& x,
		double& w)
{
  /* Compute the window function */
  dArrayT dx(x_0.Length());
  dx.DiffOf(x_0, x);
  double dist = dx.Magnitude();

  /* check out of influence range */
  if (dist > 4.0 * fDilationScaling)
    w = 0.0;
  else
  {
    double adm = fDilationScaling * fSharpeningFudgeFactor;
    double q = dist / adm;
    w = exp(-q*q)/(sqrtPi*adm);
  }
}

void GaussianWindowT::Dwindow(const dArrayT& x_0, const dArrayT& param, const dArrayT& x,
		dArrayT& Dw)
{
  /* Compute derivative of the window function */
  dArrayT dx(x_0.Length());
  dx.DiffOf(x_0, x);
  double dist = dx.Magnitude();

  /* check out of influence range */
  if (dist > 4.0 * fDilationScaling)
    Dw = 0.0;
  else
  {
    double adm = fDilationScaling * fSharpeningFudgeFactor;
    double adm2 = adm * adm;
    double q = dist / adm;
    double w = exp(-q * q) / (sqrtPi * adm);
    Dw.SetToScaled(-2.0 * w / adm2, dx);
  }
}

void GaussianWindowT::DDwindow(const dArrayT& x_0, const dArrayT& param, const dArrayT& x,
		dArrayT& DDw)
{
  /* Compute second derivative of the window function */
  dArrayT dx(x_0.Length());
  dx.DiffOf(x_0, x);
  double dist = dx.Magnitude();

  /* check out of influence range */
  if (dist > 4.0 * fDilationScaling)
    DDw = 0.0;
  else
  {
    dSymMatrixT NSDsym(x_0.Length());
    double adm = fDilationScaling * fSharpeningFudgeFactor;
    double adm2 = adm * adm;
    double q = dist / adm;
    double w = exp(-q * q) / (sqrtPi * adm);
    NSDsym.Outer(dx);
    NSDsym *= 4.0 * w / (adm2 * adm2);
    NSDsym.PlusIdentity(-2.0 * w / adm2);
    //DDw.SetColumn(0, NSDsym);
  }
}

/* multiple point calculations */
void GaussianWindowT::window(const dArrayT& x_0, const dArray2DT& param, const dArray2DT& x,
		dArrayT& w)
{
  /* compute window function for multiple field points */
  int numneighbors = x.MinorDim();        // Could be major dim!
  dArrayT dx(x_0.Length()), row(x_0.Length());
  for (int i = 0; i < numneighbors; i++)
  {
    x.RowCopy(i, row);
    dx.DiffOf(x, x_0);
    double dist = dx.Magnitude();

    /* check out of influence range */
    if (dist > 4.0 * fDilationScaling)
      w = 0.0;
    else
    {
      double adm = fDilationScaling * fSharpeningFudgeFactor;
      double q = dist / adm;
      w[i] = exp(-q*q)/(sqrtPi*adm);
    }
  }
}

void GaussianWindowT::Dwindow(const dArrayT& x_0, const dArray2DT& param, const dArray2DT& x,
		dArray2DT& Dw)
{
  /* compute derivative of window function for multiple field points */
  int numneighbors = x.MinorDim();        // Could be major dim!
  dArrayT dx(x_0.Length()), row(x_0.Length()), tempDw(x_0.Length());  
  for (int i = 0; i < numneighbors; i++)
  {
    x.RowCopy(i, row);
    dx.DiffOf(x, x_0);
    double dist = dx.Magnitude();

    /* check out of influence range */
    if (dist > 4.0 * fDilationScaling)
      Dw.SetColumn(i, 0.0);
    else
    {
      double adm = fDilationScaling * fSharpeningFudgeFactor;
      double q = dist / adm;
      double w = exp(-q * q) / (sqrtPi * adm);
      double adm2 = adm * adm;
      tempDw.SetToScaled(-2.0 * w / adm2, dx);
      Dw.SetColumn(i, tempDw);
    }
  }
}

void GaussianWindowT::DDwindow(const dArrayT& x_0, const dArray2DT& param, const dArray2DT& x,
		dArray2DT& DDw)
{
  /* compute second derivative of window function for multiple field points */
   int numneighbors = x.MinorDim();        // Could be major dim!
  dArrayT dx(x_0.Length()), row(x_0.Length());  
  dSymMatrixT NSDsym(x_0.Length());
  for (int i = 0; i < numneighbors; i++)
  {
    x.RowCopy(i, row);
    dx.DiffOf(x, x_0);
    double dist = dx.Magnitude();

    /* check out of influence range */
    if (dist > 4.0 * fDilationScaling)
      DDw.SetColumn(i, 0.0);
    else
    {
      double adm = fDilationScaling * fSharpeningFudgeFactor;
      double q = dist / adm;
      double adm2 = adm * adm;
      double w = (-q * q) / (sqrtPi * adm);
      NSDsym.Outer(dx);
      NSDsym *= 4.0 * w / (adm2 * adm2);
      NSDsym.PlusIdentity(-2.0 * w / adm2);
      DDw.SetColumn(i, NSDsym);
    }
  }
}
