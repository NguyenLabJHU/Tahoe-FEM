/* $Id: least_square.cpp,v 1.1 2004-11-12 21:09:57 paklein Exp $ */
#include "PiecewiseLinearT.h"
#include "dArray2DT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

using namespace Tahoe;

/* read list of ordered pairs */
void read_points(ifstreamT& in, dArray2DT& points, double& abs_max);

/* integrated the squared difference between the two functions between the limits given */
double squared_difference(const dArray2DT& pts1, const dArray2DT& pts2);

int main(int argc, char** argv)
{
	ofstreamT::format_stream(cout);
	cout.precision(12);

	/* echo command line arguments */
	cout << "arguments:\n";
	for (int i = 0; i < argc; i++)
		cout << setw(5) << i << ": " << argv[i] << '\n';
	cout.flush();

	/* caller */
	const char* caller = argv[0];

	/* check command line arguments */
	if (argc < 4) {
		cout << "\n usage: " << caller << " [results 1 (ref)] [results 2] [output file]\n" << endl;
		return 1;
	}

	ifstreamT in;

	/* read function 1 */
	in.open(argv[1]);	
	if (!in.is_open())
		ExceptionT::GeneralFail(caller, "could not open file \"%s\"", in.filename());
		
	dArray2DT pts1;
	double max1;
	read_points(in, pts1, max1);
	in.close();
	
//	cout << "function 1: " << max1 << '\n' << pts1 << endl;

	/* read function 2 */
	in.open(argv[2]);	
	if (!in.is_open())
		ExceptionT::GeneralFail(caller, "could not open file \"%s\"", in.filename());
		
	dArray2DT pts2;
	double max2;
	read_points(in, pts2, max2);
	in.close();
	
//	cout << "function 2: " << max2 << '\n' << pts2 << endl;

	/* construct piecewise linear functions */
	PiecewiseLinearT func1(pts1);
	PiecewiseLinearT func2(pts2);
	
	/* integration bounds */
	double x1 = (pts1(0,0) > pts2(0,0)) ? pts1(0,0) : pts2(0,0);
	int n1 = pts1.MajorDim();
	int n2 = pts2.MajorDim();
	double x2 = (pts1(n1-1,0) < pts2(n2-1,0)) ? pts1(n1-1,0) : pts2(n2-1,0);
	
	/* analytical scheme */
	double ref = (x2 - x1)*max1;
	double e_sqr = squared_difference(pts1, pts2);

	/* results */
	cout << "   x_start: " << x1 << '\n';
	cout << "     x_end: " << x2 << '\n';
	cout << "       ref: " << ref << '\n';
	cout << "   error^2: " << e_sqr << '\n';
	cout << "norm error: " << sqrt(e_sqr)/ref << endl;
	
	/* write results */
	ofstreamT out(argv[3]);
	out.precision(12);
	int d_width = OutputWidth(out,&ref);
	out << setw(d_width) << sqrt(e_sqr)/ref << ' ' << "error_norm" << endl;
	
	return 0;
}

/* read list of ordered pairs */
void read_points(ifstreamT& in, dArray2DT& points, double& abs_max)
{
	/* Results file format:
	 * x0 y0
	 * x1 y1
	 * ...
	 * x(np-1) y(np-1) */

	/* initial size */
	points.Dimension(100, 2);
	int count = 0;
	abs_max = 0.0;
	double x, y;
	in >> x >> y;
	while (in.good())
	{
		/* keep track of biggest value */
		abs_max = (fabs(y) > abs_max) ? fabs(y) : abs_max;

		/* need more memory */
		if (count+1 == points.MajorDim())
			points.Resize(2*points.MajorDim(),2);
	
		points(count,0) = x;
		points(count,1) = y;
		count++;
		
		in >> x >> y;
	}
	
	/* restore stream */
	in.clear();
	
	/* trim */
	points.Resize(count,2);
}

double squared_difference(const dArray2DT& pts1, const dArray2DT& pts2)
{
	const char caller[] = "squared_difference";
	double sum = 0.0;

	/* dimensions */
	int n1 = pts1.MajorDim();
	int n2 = pts2.MajorDim();

	/* expecting at least 2 points in each set */
	if (n1 < 2 || n2 < 2)
		ExceptionT::GeneralFail(caller, "at least 2 points needed in each set");
	
	/* smallest intervals */
	double dx1_min = pts1(n1-1,0) - pts1(0,0);
	for (int i = 1; i < pts1.MajorDim(); i++)
	{
		double dx = pts1(i,0) - pts1(i-1,0);
		dx1_min = (dx < dx1_min) ? dx : dx1_min;
	}
	double dx2_min = pts2(n2-1,0) - pts2(0,0);
	for (int i = 1; i < pts2.MajorDim(); i++)
	{
		double dx = pts2(i,0) - pts2(i-1,0);
		dx2_min = (dx < dx2_min) ? dx : dx2_min;
	}
	double small = ((dx2_min < dx1_min) ? dx2_min : dx1_min)/10.0;
	if (small < 1.0e-12)
		ExceptionT::GeneralFail(caller, "intervals too small");
		
	/* bounds of the integral - the overlap */
	double x1 = (pts1(0,0) > pts2(0,0)) ? pts1(0,0) : pts2(0,0);
	double x2 = (pts1(n1-1,0) < pts2(n2-1,0)) ? pts1(n1-1,0) : pts2(n2-1,0);
	
	int dex1 = 0;
	int dex2 = 0;
	double m1, b1;
	double m2, b2;
	
	/* advance dex1 to first point beyond x1 */
	if (pts1(dex1,0) < x1 + small)
	{
		/* next interval */
		while (dex1 < n1 && pts1(dex1,0) < x1 + small)
			dex1++;
		
		/* compute slope and intercept */
		if (dex1 < n1) {
			m1 = (pts1(dex1,1) - pts1(dex1-1,1))/(pts1(dex1,0) - pts1(dex1-1,0));
			b1 = pts1(dex1,1) - m1*pts1(dex1,0);
		}
	}
	
	/* advance dex2 to first point beyond x1 */
	if (pts2(dex2,0) < x1 + small)
	{
		/* next interval */
		while (dex2 < n2 && pts2(dex2,0) < x1 + small)
			dex2++;
		
		/* compute slope and intercept */
		if (dex2 < n2) {
			m2 = (pts2(dex2,1) - pts2(dex2-1,1))/(pts2(dex2,0) - pts2(dex2-1,0));
			b2 = pts2(dex2,1) - m2*pts2(dex2,0);
		}
	}
		
	while (dex1 < n1 && dex2 < n2)
	{
		/* interval width */
		double dx1 = pts1(dex1,0) - x1;
		double dx2 = pts2(dex2,0) - x1;
		double dx = (dx1 < dx2) ? dx1 : dx2;
	
		/* coefficience of the difference function */
		double a0 = b1*b1 - 2.0*b1*b2 + b2*b2;
		double a1 = 2.0*(b1*m1 + b2*m2 - b1*m2 - b2*m1);
		double a2 = m1*m1 - 2.0*m1*m2 + m2*m2;
		
		/* integrated */
		double c1 = a0 + a1*x1 + a2*x1*x1;
		double c2 = 0.5*a1 + a2*x1;
		double c3 = a2/3.0;
		sum += dx*(c1 + dx*(c2 + dx*c3));
		
		/* advance */
		x1 += dx;

		/* advance dex1 to first point beyond x1 */
		if (pts1(dex1,0) < x1 + small)
		{
			/* next interval */
			while (dex1 < n1 && pts1(dex1,0) < x1 + small)
				dex1++;
		
			/* compute slope and intercept */
			if (dex1 < n1) {
				m1 = (pts1(dex1,1) - pts1(dex1-1,1))/(pts1(dex1,0) - pts1(dex1-1,0));
				b1 = pts1(dex1,1) - m1*pts1(dex1,0);
			}
		}
	
		/* advance dex2 to first point beyond x1 */
		if (pts2(dex2,0) < x1 + small)
		{
			/* next interval */
			while (dex2 < n2 && pts2(dex2,0) < x1 + small)
				dex2++;
		
			/* compute slope and intercept */
			if (dex2 < n2) {
				m2 = (pts2(dex2,1) - pts2(dex2-1,1))/(pts2(dex2,0) - pts2(dex2-1,0));
				b2 = pts2(dex2,1) - m2*pts2(dex2,0);
			}
		}
		
	}
	
	return sum;
}
