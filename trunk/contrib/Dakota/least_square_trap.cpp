/* $Id: least_square_trap.cpp,v 1.1 2004-11-12 21:09:57 paklein Exp $ */
#include "PiecewiseLinearT.h"
#include "dArray2DT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

using namespace Tahoe;

/* read list of ordered pairs */
void read_points(ifstreamT& in, dArray2DT& points, double& abs_max);

/* integrated the squared difference between the two functions between the limits given */
double squared_difference(const C1FunctionT& f1, C1FunctionT& f2, double x1, double x2, int nx);

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
	if (argc < 5) {
		cout << "\n usage: " << caller << " [results 1 (ref)] [results 2] [n intervals] [output file]\n" << endl;
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
	
	/* trapezoidal scheme */
	double ref = (x2 - x1)*max1;
	int n_int = atoi(argv[3]);
	double e_sqr = squared_difference(func1, func2, x1, x2, n_int);

	/* results */
	cout << "   x_start: " << x1 << '\n';
	cout << "     x_end: " << x2 << '\n';
	cout << "     n_int: " << n_int << '\n';
	cout << "       ref: " << ref << '\n';
	cout << "   error^2: " << e_sqr << '\n';
	cout << "norm error: " << sqrt(e_sqr)/ref << endl;

	/* write results */
	ofstreamT out(argv[4]);
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

/* integrated the squared difference between the two functions between the limits given */
double squared_difference(const C1FunctionT& f1, C1FunctionT& f2, double x1, double x2, int nx)
{
//	ofstreamT debug("int.out");

	double sum = 0.0;

	double dx = (x2 - x1)/nx;
	double x = x1;
	double y1 = f1.Function(x);
	double y2 = f2.Function(x);
	double error = y2 - y1;
	double e1 = error*error;
	for (int i = 0; i < nx; i++)
	{
		x += dx;
		y1 = f1.Function(x);
		y2 = f2.Function(x);
		error = y2 - y1;
		double e2 = error*error;	
		sum += (e1 + e2)*0.5*dx; /* trapezoidal rule */

		/* output */
//		debug << x << ' ' << y1 << ' '  << y2 << ' ' << error << ' ' << e2 << '\n';

		e1 = e2;
	}
//	debug.close();
	return sum;
}
