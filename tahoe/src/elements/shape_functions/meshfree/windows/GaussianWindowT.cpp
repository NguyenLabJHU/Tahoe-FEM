/* $Id: GaussianWindowT.cpp,v 1.5 2001-06-18 18:02:53 hspark Exp $ */

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
  cout << "Dilation scaling factor = " << fDilationScaling << '\n';
  cout << "Fudge factor used = " << fSharpeningFudgeFactor << '\n';
}

/* Single point evaluations */
void GaussianWindowT::window(const dArrayT& x_n, const dArrayT& param_n, const dArrayT& x,
		int order, double& w, dArrayT& Dw, dArrayT& DDw)
{
  /* Compute second derivative of the window function */
  dArrayT dx(x.Length());
  dSymMatrixT NSDsym(x.Length());
  dx.DiffOf(x, x_n);
  double dist = dx.Magnitude();

  /* check out of influence range */
  if (!Covers(x_n, x, param_n))    // WARNING:  param_n[0] may not be correct with only 1 value!
  {
    w = 0.0;
    if (order > 0)
    {
      Dw = 0.0;
      if (order > 1)
	DDw = 0.0;
    }
  }
  else
  {
    double adm = param_n[0] * fSharpeningFudgeFactor;
    double adm2 = adm * adm;
    double q = dist / adm;
    w = exp(-q * q) / (sqrtPi * adm);
    if (order > 0)
    {
      Dw.SetToScaled(-2.0 * w / adm2, dx);
      if (order > 1)
      {
	  NSDsym.Outer(dx);
	  NSDsym *= 4.0 * w / (adm2 * adm2);
	  NSDsym.PlusIdentity(-2.0 * w / adm2);
	  DDw = NSDsym;    // THIS COULD BE WRONG!!
      }
    }
  }
}

/* multiple point calculations */
int GaussianWindowT::window(const dArray2DT& x_n, const dArray2DT& param_n, const dArrayT& x,
		int order, dArrayT& w, dArray2DT& Dw, dArray2DT& DDw)
{
  /* compute window function and derivatives for multiple field points */
  int count = 0;    // # of point covered...
  int numwindows = x_n.MinorDim();        // Could be MajorDim!
  dArrayT dx(x.Length()), row(x.Length()), tempDw(x.Length());
  dSymMatrixT NSDsym(x.Length());
  ArrayT<bool> covers(x.Length());
  Covers(x_n, x, param_n, covers);
  for (int i = 0; i < numwindows; i++)
  {
    x_n.RowCopy(i, row);
    dx.DiffOf(x, x_n);
    double dist = dx.Magnitude();

    /* check out of influence range */
    if (!covers[i])
    {
      w[i] = 0.0;
      if (order > 0)
      {
	Dw.SetColumn(i, 0.0);
	if (order > 1)
	  DDw.SetColumn(i, 0.0);
      }
    }
    else
    {
      double adm = param_n[i] * fSharpeningFudgeFactor;
      double q = dist / adm;
      double adm2 = adm * adm;
      w[i] = (-q * q) / (sqrtPi * adm);
      if (order > 0)
      {
	tempDw.SetToScaled(-2.0 * w[i] / adm2, dx);
	Dw.SetColumn(i, tempDw);
	if (order > 1)
	{
	  NSDsym.Outer(dx);
	  NSDsym *= 4.0 * w[i] / (adm2 * adm2);
	  NSDsym.PlusIdentity(-2.0 * w[i] / adm2);
	  DDw.SetColumn(i, NSDsym);
	}
      }
      count++;
    }
  }
  return count;
}

bool GaussianWindowT::Covers(const dArrayT& x_n, const dArrayT& x, const dArrayT& param_n)
{
  /* is this necessary? */
  dArrayT dx(x.Length());
  dx.DiffOf(x, x_n);
  double dist = dx.Magnitude();
  if (dist > 4.0 * param_n[0])
    return false;
  else
    return true;
}

void GaussianWindowT::Covers(const dArray2DT& x_n, const dArrayT& x, 
			     const dArray2DT& param_n, const ArrayT<bool>& covers)
{
  int count = 0;    // # of point covered...
  int numwindows = x_n.MinorDim();        // Could be MajorDim!
  dArrayT dx(x.Length()), row(x.Length()), tempDw(x.Length());
  dSymMatrixT NSDsym(x.Length());
  for (int i = 0; i < numwindows; i++)
  {
    x_n.RowCopy(i, row);
    dx.DiffOf(x, x_n);
    double dist = dx.Magnitude();

    /* check out of influence range */
    if (dist > 4.0 * param_n[i])
      covers[i] = false;
    else
      covers[i] = true;
  }
}
