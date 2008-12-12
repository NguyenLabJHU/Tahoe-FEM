/* $Id: QuadT.cpp,v 1.22 2008-12-12 17:45:22 lxmota Exp $ */
/* created: paklein (07/03/1996) */
#include "QuadT.h"
#include <math.h>
#include "ExceptionT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"
#include "LocalArrayT.h"

using namespace Tahoe;

/* parameters */
const int kQuadnsd          = 2;
const int kNumVertexNodes	= 4;

const double sqrt3 = sqrt(3.0);

/* cannonical nodal coordinates */
const double ra[9] = {-1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 0.0,-1.0, 0.0};
const double sa[9] = {-1.0,-1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 0.0, 0.0};

/* constructor */
QuadT::QuadT(int numnodes): GeometryBaseT(numnodes, kNumVertexNodes)
{
	const char caller[] = "HexahedronT::HexahedronT";
	fCoords.Dimension(2,numnodes);
	double* x = fCoords(0);
	double* y = fCoords(1);
	const double ra9[9] = {-1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 0.0, -1.0, 0.0};
	const double sa9[9] = {-1.0,-1.0, 1.0, 1.0, -1.0, 0.0, 1.0, 0.0,0.0};
	for (int i = 0; i< numnodes; i++)
	{
		x[i] = ra9[i];
		y[i] = sa9[i];
	}
}
/* evaluate the shape functions and gradients. */
void QuadT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na) const
{
  const char caller[] = "QuadT::EvaluateShapeFunctions";

#if __option(extended_errorcheck)
  if (coords.Length() != 2 || Na.Length() != fNumNodes) ExceptionT::SizeMismatch(
      caller);

  if (fNumNodes != 1 &&
      fNumNodes != kNumVertexNodes &&
      fNumNodes != 8 &&
      fNumNodes != 9) {
    ExceptionT::GeneralFail(caller);
  }

#endif

  /* coordinates */
  double r = coords[0];
  double s = coords[1];
#if 0
  if (r < -1.0 || r> 1.0) ExceptionT::OutOfRange(caller);
  if (s < -1.0 || s> 1.0) ExceptionT::OutOfRange(caller);
#endif

  /* destinations */
  double* na = Na.Pointer();

  if (fNumNodes == 1) {

    na[0] = 1.0;

  } else if (fNumNodes == kNumVertexNodes) {

    // corner nodes
    for (int lnd = 0; lnd < kNumVertexNodes; lnd++) {
      double tempr1 = 1.0 + ra[lnd] * r;
      double temps1 = 1.0 + sa[lnd] * s;
      *na++ = 0.25 * tempr1 * temps1;
    }

  } else if (fNumNodes == 8) {

    // corner nodes
    for (int lnd = 0; lnd < kNumVertexNodes; lnd++) {
      double tempr1 = 1.0 + ra[lnd] * r;
      double temps1 = 1.0 + sa[lnd] * s;
      *na++ = 0.25 * tempr1 * temps1;
    }

    // mid-side nodes
    for (int lnd = kNumVertexNodes; lnd < fNumNodes; lnd++) {
      int off1 = -kNumVertexNodes;
      int off2 = off1 + 1;
      if (lnd == 7) /* the 8th node is the 0th */
      off2 = -7;

      double fac = 1.0;
      if (lnd > 5) /* the 7th and 8th nodes */
      fac = -1.0;

      double tempr1 = 1.0 + r * sa[lnd];
      double temps1 = 1.0 + s * ra[lnd];
      double tempr2 = 1.0 + r * fac;
      double temps2 = 1.0 - s * fac;
      double shape = 0.5 * tempr1 * tempr2 * temps1 * temps2;

      /* set mid-side node functions */
      *na = shape;

      /* correct the corner node functions */
      shape *= 0.5;
      *(na + off1) -= shape;
      *(na + off2) -= shape;

      /* next */
      na++;
    }

  } else if (fNumNodes == 9) {

    double lr1 = 0.5 * r * (r - 1);
    double lr2 = 1.0 - r * r;
    double lr3 = 0.5 * r * (r + 1);

    double ls1 = 0.5 * s * (s - 1);
    double ls2 = 1.0 - s * s;
    double ls3 = 0.5 * s * (s + 1);

    na[0] = lr1 * ls1;
    na[1] = lr3 * ls1;
    na[2] = lr3 * ls3;
    na[3] = lr1 * ls3;
    na[4] = lr2 * ls1;
    na[5] = lr3 * ls2;
    na[6] = lr2 * ls3;
    na[7] = lr1 * ls2;
    na[8] = lr2 * ls2;

  } else {

    ExceptionT::GeneralFail(caller);

  }

}

/* evaluate the shape functions and gradients. */
void QuadT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na,
    dArray2DT& DNa) const
{
  const char caller[] = "HexahedronT::EvaluateShapeFunctions";

#if __option(extended_errorcheck)
  if (coords.Length() != 2 || Na.Length() != fNumNodes || DNa.MajorDim() != 2
      || DNa.MinorDim() != fNumNodes) ExceptionT::SizeMismatch(caller);

  if (fNumNodes != 1 && fNumNodes != kNumVertexNodes && fNumNodes != 8
      && fNumNodes != 9) {
    ExceptionT::GeneralFail(caller);
  }

#endif

  /* coordinates */
  double r = coords[0];
  double s = coords[1];

#if 0
  if (r < -1.0 || r> 1.0) ExceptionT::OutOfRange(caller, "r = %15.12e", r);
  if (s < -1.0 || s> 1.0) ExceptionT::OutOfRange(caller, "s = %15.12e", s);
#endif

  /* destinations */
  double* na = Na.Pointer();
  double* nax = DNa(0);
  double* nay = DNa(1);

  if (fNumNodes == 1) {

    na[0] = 1.0;
    nax[0] = 0.0;
    nay[0] = 0.0;

  } else if (fNumNodes == kNumVertexNodes) {

    for (int lnd = 0; lnd < kNumVertexNodes; lnd++) {
      double tempr1 = 1.0 + ra[lnd] * r;
      double temps1 = 1.0 + sa[lnd] * s;

      *na++ = 0.25 * tempr1 * temps1;
      *nax++ = 0.25 * ra[lnd] * temps1;
      *nay++ = 0.25 * tempr1 * sa[lnd];
    }

  } else if (fNumNodes == 8) {

    for (int lnd = 0; lnd < kNumVertexNodes; lnd++) {
      double tempr1 = 1.0 + ra[lnd] * r;
      double temps1 = 1.0 + sa[lnd] * s;

      *na++ = 0.25 * tempr1 * temps1;
      *nax++ = 0.25 * ra[lnd] * temps1;
      *nay++ = 0.25 * tempr1 * sa[lnd];
    }

    // mid-side nodes
    for (int lnd = kNumVertexNodes; lnd < fNumNodes; lnd++) {
      int off1 = -kNumVertexNodes;
      int off2 = off1 + 1;
      if (lnd == 7) /* the 8th node is the 0th */
      off2 = -7;

      double fac = 1.0;
      if (lnd > 5) /* the 7th and 8th nodes */
      fac = -1.0;

      double tempr1 = 1.0 + r * sa[lnd];
      double temps1 = 1.0 + s * ra[lnd];
      double tempr2 = 1.0 + r * fac;
      double temps2 = 1.0 - s * fac;

      double shape = 0.5 * tempr1 * tempr2 * temps1 * temps2;
      double shapex, shapey;

      if ((lnd + 1) % 2 == 0) {
        /* even numbered mid-side nodes */
        shapex = 0.5 * ra[lnd] * temps1 * temps2;
        shapey = -s * tempr2;
      } else {
        /* odd numbered mid-side nodes */
        shapex = -r * temps2;
        shapey = 0.5 * sa[lnd] * tempr1 * tempr2;
      }

      /* set mid-side node functions */
      *na = shape;
      *nax = shapex;
      *nay = shapey;

      /* correct the corner node functions */
      shapex *= 0.5;
      shapey *= 0.5;
      shape *= 0.5;

      *(na + off1) -= shape;
      *(na + off2) -= shape;
      *(nax + off1) -= shapex;
      *(nax + off2) -= shapex;
      *(nay + off1) -= shapey;
      *(nay + off2) -= shapey;

      na++;
      nax++;
      nay++;
    }

  } else if (fNumNodes == 9) {

    double lr1 = 0.5 * r * (r - 1.0);
    double lr2 = 1.0 - r * r;
    double lr3 = 0.5 * r * (r + 1.0);

    double ls1 = 0.5 * s * (s - 1.0);
    double ls2 = 1.0 - s * s;
    double ls3 = 0.5 * s * (s + 1.0);

    na[0] = lr1 * ls1;
    nax[0] = (r - 0.5) * ls1;
    nay[0] = lr1 * (s - 0.5);
    na[1] = lr3 * ls1;
    nax[1] = (r + 0.5) * ls1;
    nay[1] = lr3 * (s - 0.5);
    na[2] = lr3 * ls3;
    nax[2] = (r + 0.5) * ls3;
    nay[2] = lr3 * (s + 0.5);
    na[3] = lr1 * ls3;
    nax[3] = (r - 0.5) * ls3;
    nay[3] = lr1 * (s + 0.5);
    na[4] = lr2 * ls1;
    nax[4] = (-2.0 * r) * ls1;
    nay[4] = lr2 * (s - 0.5);
    na[5] = lr3 * ls2;
    nax[5] = (r + 0.5) * ls2;
    nay[5] = lr3 * (-2.0 * s);
    na[6] = lr2 * ls3;
    nax[6] = (-2.0 * r) * ls3;
    nay[6] = lr2 * (s + 0.5);
    na[7] = lr1 * ls2;
    nax[7] = (r - 0.5) * ls2;
    nay[7] = lr1 * (-2.0 * s);
    na[8] = lr2 * ls2;
    nax[8] = (-2.0 * r) * ls2;
    nay[8] = lr2 * (-2.0 * s);

  } else {

    ExceptionT::GeneralFail(caller);
  }

}

/* evaluate the shape functions and first and second gradients. */
void QuadT::EvaluateShapeFunctions(const dArrayT& coords, dArrayT& Na, dArray2DT& DNa, dArray2DT& DDNa) const
{
#pragma unused(coords)
#pragma unused(Na)
#pragma unused(DNa)
#pragma unused(DDNa)

	cout << "\n QuadT::EvaluateShapeFunctions: not implemented" << endl;
	throw ExceptionT::kGeneralFail;
}

/* sets first and second derivative of shape functions */
void QuadT::SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x, ArrayT<dArray2DT>& Na_xx, dArrayT& weights) const
{
#pragma unused(Na)
#pragma unused(Na_x)
#pragma unused(Na_xx)
#pragma unused(weights)

	cout << "\n QuadT::SetLocalShape: not implemented" << endl;
	throw ExceptionT::kGeneralFail;
}


/* compute local shape functions and derivatives */
void QuadT::SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x,
	dArrayT& weights) const
{
	const char caller[] = "QuadT::SetLocalShape";

	/* dimensions */
	int numnodes  = Na.MinorDim();
	int numint    = weights.Length();
	int nsd       = Na_x[0].MajorDim();

	/* dimension checks */
	if (numnodes < 1 || numnodes >9)
		ExceptionT::GeneralFail(caller, "unsupported number of element nodes: %d", numnodes);

	if (numint != 1 &&
	    numint != 4 &&
	    numint != 5 &&
	    numint != 9 &&
	    numint != 16)
		ExceptionT::GeneralFail(caller, "unsupported number of integration points: %d", numint);

	if (numnodes == 9 && (numint !=4 && numint != 9))
		ExceptionT::GeneralFail(caller, "number of integration points is not supported by quadratic Lagrange element: %d", numint);


	if (nsd != kQuadnsd) ExceptionT::GeneralFail(caller);

	/* initialize */
	Na = 0.0;
	for (int i = 0; i < Na_x.Length(); i++)
		Na_x[i] = 0.0;

	double	ra_5[5] = {-1.0, 1.0, 1.0,-1.0, 0.0};
	double  sa_5[5] = {-1.0,-1.0, 1.0, 1.0, 0.0};

	double xa_16[16], ya_16[16];
	const double *xa, *ya;
	double 	g;

	/* integration weights */
	switch (numint)
	{
		case 1:
		{
			xa = ra;
			ya = sa;
			g = 0.0;
			weights[0] = 4.0;
			break;
		}
		case 4:
		{
			xa = ra;
			ya = sa;
			g = 1.0/sqrt3;

			weights[0] = 1.0;
			weights[1] = 1.0;
			weights[2] = 1.0;
			weights[3] = 1.0;
			break;
		}
		case 5:
		{
			xa = ra_5;
			ya = sa_5;
			g = sqrt(3.0/5.0);

			double w0 = 16.0/9.0;
			double w1 = 5.0/9.0;

			weights[0] = w1;
			weights[1] = w1;
			weights[2] = w1;
			weights[3] = w1;
			weights[4] = w0;
			break;
		}
		case 9:
		{
			xa = ra;
			ya = sa;
			g = sqrt(3.0/5.0);

			double a = 25.0/81.0;
			double b = 40.0/81.0;
			double c = 64.0/81.0;

			weights[0] = a;
			weights[1] = a;
			weights[2] = a;
			weights[3] = a;
			weights[4] = b;
			weights[5] = b;
			weights[6] = b;
			weights[7] = b;
			weights[8] = c;
			break;
		}

		case 16:
		{
			/* coordinates */
			double b1 = sqrt((3.0 - 2.0*sqrt(6.0/5.0))/7.0);
			double b2 = sqrt((3.0 + 2.0*sqrt(6.0/5.0))/7.0);
			double b_1D[4] = {-b2,-b1, b1, b2};

			/* weights */
			double w1 = (18.0 + sqrt(30.0))/36.0;
			double w2 = (18.0 - sqrt(30.0))/36.0;
			double w_1D[4] = {w2, w1, w1, w2};

			int x_i = 0;
			int y_i = 0;
			for (int i = 0; i < 16; i++)
			{
				xa_16[i]   = b_1D[x_i];
				ya_16[i]   = b_1D[y_i];
				weights[i] = w_1D[x_i]*w_1D[y_i];

				if (++y_i == 4)
				{
					x_i++;
					y_i = 0;
				}
			}

			xa = xa_16;
			ya = ya_16;
			g  = 1.0;
			break;
		}
		default:

			ExceptionT::GeneralFail(caller);
	}

	/* evaluate shape functions at integration points */
	dArrayT coords(2), na;
	for (int i = 0; i < numint; i++)
	{
		/* ip coordinates */
		coords[0] = g*xa[i];
		coords[1] = g*ya[i];

		/* shape function values */
		Na.RowAlias(i, na);

		/* evaluate (static binding) */
		QuadT::EvaluateShapeFunctions(coords, na, Na_x[i]);
	}
}

/* compute gradients of the "bubble" modes */
void QuadT::BubbleModeGradients(ArrayT<dArray2DT>& Na_x) const
{
	const char caller[] = "QuadT::BubbleModeGradients";

	/* limit integration rules */
	int nip = Na_x.Length();
	if (nip != 4 && nip != 5)
		ExceptionT::GeneralFail(caller, "only 4 or 5 point rules defined: %d", nip);

	/* integration rules */
	double ra[5] = {-1.0, 1.0, 1.0,-1.0, 0.0};
	double sa[5] = {-1.0,-1.0, 1.0, 1.0, 0.0};
	double g = (nip == 4) ? 1.0/sqrt3 : sqrt(3.0/5.0);

	/* set gradients */
	for (int i = 0; i < nip; i++)
	{
		dArray2DT& na_x = Na_x[i];
		if (na_x.MajorDim() != 2 || na_x.MinorDim() != 2)
			ExceptionT::GeneralFail(caller, "gradients array must be 2x2: %dx%d",
				na_x.MajorDim(), na_x.MinorDim());

		/* integration point coordinates */
		double r = g*ra[i];
		double s = g*sa[i];

		/* set gradients */
		double* nax = na_x(0);
		double* nay = na_x(1);
		nax[0] = r;
		nax[1] = 0.0;
		nay[0] = 0.0;
		nay[1] = s;
	}
}

/* set the values of the nodal extrapolation matrix */
void QuadT::SetExtrapolation(dMatrixT& extrap) const
{
	const char caller[] = "QuadT::SetExtrapolation";

	/* dimensions */
	int numnodes = extrap.Rows();
	int numint   = extrap.Cols();

	/* dimension checks */
	if (numnodes < 1 || numnodes > 9) ExceptionT::GeneralFail(caller);

	if(numnodes == 9 && (numint !=4 && numint != 9))
		ExceptionT::GeneralFail(caller, "requested number of integration points is not supported by quadratic Lagrange element: %d", numint);

	/* initialize */
	extrap = 0.0;

	switch (numint)
	{
		case 1:

		extrap = 1.0;
		break;

		case 4:
		{
			double a = 1.0 + sqrt3/2.0;
			double b =-1.0/2.0;
			double c = 1.0 - sqrt3/2.0;

			/* vertex nodes */
			extrap(0,0) = a;
			extrap(0,1) = b;
			extrap(0,2) = c;
			extrap(0,3) = b;

			extrap(1,0) = b;
			extrap(1,1) = a;
			extrap(1,2) = b;
			extrap(1,3) = c;

			extrap(2,0) = c;
			extrap(2,1) = b;
			extrap(2,2) = a;
			extrap(2,3) = b;

			extrap(3,0) = b;
			extrap(3,1) = c;
			extrap(3,2) = b;
			extrap(3,3) = a;

			/* midside nodes */
			if (numnodes > kNumVertexNodes)
			{
				double d = (1.0 + sqrt3)/4.0;
				double e = (1.0 - sqrt3)/4.0;
                double f = 0.25;

				double smooth[5][4] = {{d,d,e,e},
				                       {e,d,d,e},
				                       {e,e,d,d},
				                       {d,e,e,d},
				                       {f,f,f,f}};

				for (int i = kNumVertexNodes; i < numnodes; i++)
					for (int j = 0; j < numint; j++)
						extrap(i,j) = smooth[i-kNumVertexNodes][j];
			}

			break;
		}
		case 5:
		{
			/* data for the smoothing  matrix */
			double smooth_5[8*5] = {
 3.8660254037844     ,1.5                 ,2.1339745962156     ,1.5                 ,-8.
,1.5                 ,3.8660254037844     ,1.5                 ,2.1339745962156     ,-8.
,2.1339745962156     ,1.5                 ,3.8660254037844     ,1.5                 ,-8.
,1.5                 ,2.1339745962156     ,1.5                 ,3.8660254037844     ,-8.
,0.43301270189222    ,0.43301270189222    ,-0.43301270189222   ,-0.43301270189222   ,1.
,-0.43301270189222   ,0.43301270189222    ,0.43301270189222    ,-0.43301270189222   ,1.
,-0.43301270189222   ,-0.43301270189222   ,0.43301270189222    ,0.43301270189222    ,1.
,0.43301270189222    ,-0.43301270189222   ,-0.43301270189222   ,0.43301270189222    ,1.};

			/* copy in */
			dMatrixT smooth(8, 5, smooth_5);
			smooth.CopyBlock(0, 0, extrap);
			break;
		}
		case 9:
		{
			/* constants */
			double a = 5.0/18.0;
			double b = (20.0 + 5.0*sqrt(15.0))/18.0;
			double c = (20.0 - 5.0*sqrt(15.0))/18.0;

			double d = (-10.0 + 2.0*sqrt(15.0))/18.0;
			double e = (-10.0 - 2.0*sqrt(15.0))/18.0;

			/* vertex nodes */
			extrap(0,0) = b;
			extrap(0,1) = a;
			extrap(0,2) = c;
			extrap(0,3) = a;
			extrap(0,4) = e;
			extrap(0,5) = d;
			extrap(0,6) = d;
			extrap(0,7) = e;
			extrap(0,8) = 4.0/9.0;

			extrap(1,0) = a;
			extrap(1,1) = b;
			extrap(1,2) = a;
			extrap(1,3) = c;
			extrap(1,4) = e;
			extrap(1,5) = e;
			extrap(1,6) = d;
			extrap(1,7) = d;
			extrap(1,8) = 4.0/9.0;

			extrap(2,0) = c;
			extrap(2,1) = a;
			extrap(2,2) = b;
			extrap(2,3) = a;
			extrap(2,4) = d;
			extrap(2,5) = e;
			extrap(2,6) = e;
			extrap(2,7) = d;
			extrap(2,8) = 4.0/9.0;

			extrap(3,0) = a;
			extrap(3,1) = c;
			extrap(3,2) = a;
			extrap(3,3) = b;
			extrap(3,4) = d;
			extrap(3,5) = d;
			extrap(3,6) = e;
			extrap(3,7) = e;
			extrap(3,8) = 4.0/9.0;

			/* midside nodes */
			if (numnodes > kNumVertexNodes)
			{
				double f = (5.0 + sqrt(15.0))/6.0;
				double g = (5.0 - sqrt(15.0))/6.0;
				double h =-2.0/3.0;

				double smooth[5][9] = {{0,0,0,0,f,0,g,0,h},
				                       {0,0,0,0,0,f,0,g,h},
				                       {0,0,0,0,g,0,f,0,h},
				                       {0,0,0,0,0,g,0,f,h},
				                       {0,0,0,0,0,0,0,0,1}};

				for (int i = kNumVertexNodes; i < numnodes; i++)
					for (int j = 0; j < numint; j++)
						extrap(i,j) = smooth[i-kNumVertexNodes][j];
			}
			break;
		}

		case 16:
		{
			/* data for the smoothing  matrix */
			double smooth_16[8*16] = {
2.3310819800373e+00, -1.7392742256873e-01,  1.2977127608750e-02, -1.7392742256873e-01,
-1.4096315416317e-01,  1.0517587236621e-02,  1.0517587236621e-02, -1.4096315416317e-01,
-1.2422443623633e+00,  9.2686727449598e-02, -4.5653628771612e-02,  6.1187793035203e-01,
7.5119916442126e-02, -6.7476185377616e-02, -3.7000947956309e-02,  9.0435721689180e-01,
6.1187793035203e-01, -4.5653628771610e-02,  9.2686727449599e-02, -1.2422443623634e+00,
-3.7000947956310e-02, -6.7476185377616e-02,  7.5119916442126e-02,  9.0435721689180e-01,
-1.7392742256873e-01,  1.2977127608750e-02, -1.7392742256873e-01,  2.3310819800373e+00,
1.0517587236621e-02,  1.0517587236621e-02, -1.4096315416317e-01, -1.4096315416317e-01,
-1.2422443623634e+00,  6.1187793035203e-01, -4.5653628771615e-02,  9.2686727449597e-02,
9.0435721689180e-01, -3.7000947956310e-02, -6.7476185377616e-02,  7.5119916442126e-02,
6.6199776285810e-01, -3.2607257743128e-01,  1.6060979616250e-01, -3.2607257743128e-01,
-4.8193614118559e-01,  2.3738170811213e-01,  2.3738170811213e-01, -4.8193614118559e-01,
-3.2607257743127e-01,  1.6060979616250e-01, -3.2607257743127e-01,  6.6199776285809e-01,
2.3738170811213e-01,  2.3738170811213e-01, -4.8193614118559e-01, -4.8193614118559e-01,
9.2686727449598e-02, -4.5653628771611e-02,  6.1187793035203e-01, -1.2422443623633e+00,
-6.7476185377616e-02, -3.7000947956309e-02,  9.0435721689180e-01,  7.5119916442126e-02,
6.1187793035203e-01, -1.2422443623634e+00,  9.2686727449600e-02, -4.5653628771612e-02,
9.0435721689180e-01,  7.5119916442126e-02, -6.7476185377617e-02, -3.7000947956310e-02,
-3.2607257743128e-01,  6.6199776285809e-01, -3.2607257743127e-01,  1.6060979616251e-01,
-4.8193614118559e-01, -4.8193614118559e-01,  2.3738170811214e-01,  2.3738170811214e-01,
1.6060979616250e-01, -3.2607257743127e-01,  6.6199776285810e-01, -3.2607257743127e-01,
2.3738170811214e-01, -4.8193614118559e-01, -4.8193614118559e-01,  2.3738170811214e-01,
-4.5653628771611e-02,  9.2686727449597e-02, -1.2422443623634e+00,  6.1187793035203e-01,
-6.7476185377616e-02,  7.5119916442125e-02,  9.0435721689180e-01, -3.7000947956310e-02,
-1.7392742256873e-01,  2.3310819800373e+00, -1.7392742256873e-01,  1.2977127608749e-02,
-1.4096315416317e-01, -1.4096315416317e-01,  1.0517587236621e-02,  1.0517587236621e-02,
9.2686727449598e-02, -1.2422443623633e+00,  6.1187793035203e-01, -4.5653628771610e-02,
7.5119916442126e-02,  9.0435721689180e-01, -3.7000947956308e-02, -6.7476185377616e-02,
-4.5653628771610e-02,  6.1187793035203e-01, -1.2422443623634e+00,  9.2686727449597e-02,
-3.7000947956310e-02,  9.0435721689180e-01,  7.5119916442125e-02, -6.7476185377616e-02,
1.2977127608750e-02, -1.7392742256873e-01,  2.3310819800373e+00, -1.7392742256873e-01,
1.0517587236621e-02, -1.4096315416317e-01, -1.4096315416317e-01,  1.0517587236621e-02};

				/* copy in */
				dMatrixT smooth(8,16, smooth_16);
				smooth.CopyBlock(0, 0, extrap);
			break;
		}
		default:
			ExceptionT::GeneralFail(caller, "no nodal extrapolation with Gauss rule: %d", numint);
	}
}

/* integration point gradient matrix */
void QuadT::IPGradientTransform(int ip, dMatrixT& transform) const
{
	const char caller[] = "QuadT::IPGradientTransform";

	/* dimensions */
	int nsd = transform.Rows();
	int nip = transform.Cols();
	if (nsd != 2) ExceptionT::SizeMismatch(caller);

	if (nip == 1)
		transform = 0.0;
	else if (nip == 4) {
		double a = sqrt(3.0)/2.0;
		double m0[2*4] = {-a, -a, a, 0.0, 0.0, 0.0, 0.0, a};
		double m1[2*4] = {-a, 0.0, a, -a, 0.0, a, 0.0, 0.0};
		double m2[2*4] = {0.0, 0.0, 0.0, -a, a, a, -a, 0.0};
		double m3[2*4] = {0.0, -a, 0.0, 0.0, a, 0.0, -a, a};
		double* m[4] = {m0, m1, m2, m3};
		ArrayT<double*> m_array(4, m);
		dMatrixT trans(2, 4, m_array[ip]);
		transform = trans;
	}
	else if (nip == 9) {

		double a = 3.0*sqrt(3.0/5.0)/2.0;
		double b = 2.0*sqrt(3.0/5.0);
		double c = sqrt(3.0/5.0)/2.0;

		double m0[2*9] = {-a, -a, -c, 0, 0, 0, 0, -c, b, 0, 0, 0, 0, 0, 0, b, 0, 0};
		double m1[2*9] = {c, 0, a, -a, 0, -c, 0, 0, -b, 0, 0, b, 0, 0, 0, 0, 0, 0};
		double m2[2*9] = {0, 0, 0, c, a, a, c, 0, 0, 0, 0, -b, -b, 0, 0, 0, 0, 0};
		double m3[2*9] = {0, c, 0, 0, -c, 0, -a, a, 0, 0, 0, 0, b, 0, 0, -b, 0, 0};
		double m4[2*9] = {-c, 0, c, 0, 0, 0, 0, 0, 0, -a, 0, 0, 0, -c, 0, 0, 0, b};
		double m5[2*9] = {0, 0, 0, -c, 0, c, 0, 0, 0, 0, a, 0, 0, 0, c, 0, -b, 0};
		double m6[2*9] = {0, 0, 0, 0, c, 0, -c, 0, 0, c, 0, 0, 0, a, 0, 0, 0, -b};
		double m7[2*9] = {0, -c, 0, 0, 0, 0, 0, c, 0, 0, -c, 0, 0, 0, -a, 0, b, 0};
		double m8[2*9] = {0, 0, 0, 0, 0, 0, 0, 0, 0, -c, c, 0, 0, c, -c, 0, 0, 0};
		double* m[9] = {m0, m1, m2, m3, m4, m5, m6, m7, m8};
		ArrayT<double*> m_array(9, m);
		dMatrixT trans(2, 9, m_array[ip]);
		transform = trans;
	}
	else
		ExceptionT::GeneralFail(caller, "unsupported number of integration points %d", nip);
}

/* return the local node numbers for each facet of the element
* numbered to produce at outward normal in the order: vertex
* nodes, mid-edge nodes, mid-face nodes */
void QuadT::NodesOnFacet(int facet, iArrayT& facetnodes) const
{
// TEMP: not implemented with midside nodes
	const char caller[] = "QuadT::NodesOnFacet";
	if (fNumNodes != 4 && fNumNodes != 8 && fNumNodes != 9)
		ExceptionT::GeneralFail(caller, "only implemented for 4, 8, and 9 element nodes: %d", fNumNodes);

#if __option(extended_errorcheck)
	if (facet < 0 || facet > 3) ExceptionT::OutOfRange(caller);
#endif

	/* nodes-facet data */
	int dat4[] = {0,1,1,2,2,3,3,0};
	int dat8[] = {0,1,4,1,2,5,2,3,6,3,0,7};

	/* collect facet data */
	iArrayT tmp;
	if (fNumNodes == 4)
		tmp.Set(2, dat4 + facet*2);
	else
		tmp.Set(3, dat8 + facet*3);

	/* (allocate and) copy in */
	facetnodes = tmp;
}

void QuadT::NumNodesOnFacets(iArrayT& num_nodes) const
{
// TEMP: not implemented with midside nodes
	if (fNumNodes != 4 && fNumNodes != 8 && fNumNodes != 9)
		ExceptionT::GeneralFail("QuadT::NumNodesOnFacet", "only implemented for 4, 8, and 9 element nodes: %d", fNumNodes);

	num_nodes.Dimension(4);
	if (fNumNodes == 4)
		num_nodes = 2;
	else
		num_nodes = 3;
}

/* return the local node numbers for each edge of element */
void QuadT::NodesOnEdges(iArray2DT& nodes_on_edges) const {
	nodes_on_edges.Dimension(0,0);
}

/* returns the nodes on each facet needed to determine neighbors
* across facets */
void QuadT::NeighborNodeMap(iArray2DT& facetnodes) const
{
	int dat4[] = {0,1,1,2,2,3,3,0};
	iArray2DT temp(4, 2, dat4);

	facetnodes = temp;
}

/* return geometry and number of nodes on each facet */
void QuadT::FacetGeometry(ArrayT<CodeT>& facet_geom, iArrayT& facet_nodes) const
{
	facet_geom.Dimension(fNumFacets);
	facet_geom = kLine;

	facet_nodes.Dimension(fNumFacets);
	facet_nodes = 2;
/*	for (int i = 0; i < (fNumNodes - kNumVertexNodes); i++)
		facet_nodes[i] = 3;*/
	if (fNumNodes==8||fNumNodes==9)
	    facet_nodes=3;
}

/* return true if the given point is within the domain */
bool QuadT::PointInDomain(const LocalArrayT& coords, const dArrayT& point) const
{
	/* method: run around the perimeter of the element and see if
	 *         the point always lies to the left of segment a-b */
	int nen = coords.NumberOfNodes();
	int a = nen - 1;
	int b = 0;
	bool in_domain = true;
	for (int i = 0; in_domain && i < nen; i++)
	{
		double ab_0 = coords(b,0) - coords(a,0);
		double ab_1 = coords(b,1) - coords(a,1);
		double L2 = ab_0*ab_0 + ab_1*ab_1;

		double ap_0 = point[0] - coords(a,0);
		double ap_1 = point[1] - coords(a,1);

		double cross = ab_0*ap_1 - ab_1*ap_0;
		in_domain = (cross/L2) > -kSmall;
		a++;
		b++;
		if (a == nen) a = 0;
	}
	return in_domain;
}

/* return the integration point whose domain contains the given point in the
 * parent domain coordinates */
int QuadT::IPDomain(int nip, const dArrayT& coords) const
{
	const char caller[] = "QuadT::IPDomain";

	/* domain check */
	if (coords[0] < -1.0 || coords[0] > 1.0 || coords[1] < -1.0 || coords[1] > 1.0)
		ExceptionT::OutOfRange(caller, "{%g,%g} outside domain", coords[0], coords[1]);

	if (nip == 1)
		return 0;
	else if (nip == 4) {
		if (coords[1] > 0.0) /* upper */ {
			if (coords[0] > 0.0) /* right */
				return 2;
			else /* left */
				return 3;
		}
		else /* lower-half */ {
			if (coords[0] > 0.0) /* right */
				return 1;
			else /* left */
				return 0;
		}
	}
	else
		ExceptionT::GeneralFail(caller, "%d integration points not supported", nip);

	/* dummy */
	return -1;
}

/* subdomain geometry */
GeometryT::CodeT QuadT::NodalSubDomainGeometry(void) const
{
	/* limited support */
	if (fNumNodes != 4)
		ExceptionT::GeneralFail("QuadT::NodalSubDomainGeometry",
			"unsupported number of nodes %d", fNumNodes);

	return GeometryT::kQuadrilateral;
}

/* number of nodes defining the nodal subdomain */
int QuadT::NodalSubDomainNumPoints(void) const
{
	/* limited support */
	if (fNumNodes != 4)
		ExceptionT::GeneralFail("QuadT::NodalSubDomainGeometry",
			"unsupported number of nodes %d", fNumNodes);

	return 4;
}

/* compute the coordinates of the points defining the nodal subdomain */
void QuadT::NodalSubDomainCoordinates(const LocalArrayT& coords, int node,
	LocalArrayT& subdomain_coords) const
{
	const char caller[] = "QuadT::NodalSubDomainCoordinates";

#if __option(extended_errorcheck)
	/* limited support */
	if (fNumNodes != 4)
		ExceptionT::GeneralFail(caller, "unsupported number of nodes %d", fNumNodes);

	/* checks */
	if (coords.NumberOfNodes() != fNumNodes ||
		coords.MinorDim() != 2 ||
		node < 0 || node >= fNumNodes ||
		subdomain_coords.MinorDim() != 2 ||
		subdomain_coords.NumberOfNodes() != QuadT::NodalSubDomainNumPoints())
		ExceptionT::SizeMismatch(caller);
#endif

	int next_node[4] = {1,2,3,0};
	int back_node[4] = {3,0,1,2};
	int diag_node[4] = {2,3,0,1};

	const double* px = coords(0);
	const double* py = coords(1);

	int back = back_node[node];
	subdomain_coords(back,0) = 0.5*(px[node] + px[back]);
	subdomain_coords(back,1) = 0.5*(py[node] + py[back]);

	subdomain_coords(node,0) = px[node];
	subdomain_coords(node,1) = py[node];

	int next = next_node[node];
	subdomain_coords(next,0) = 0.5*(px[node] + px[next]);
	subdomain_coords(next,1) = 0.5*(py[node] + py[next]);

	int diag = diag_node[node];
	subdomain_coords(diag,0) = 0.25*(px[0] + px[1] + px[2] + px[3]);
	subdomain_coords(diag,1) = 0.25*(py[0] + py[1] + py[2] + py[3]);
}
