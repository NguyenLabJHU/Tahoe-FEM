/* $Id: area_diff.cpp,v 1.1 2005-02-17 17:13:34 paklein Exp $ */
#include "PiecewiseLinearT.h"
#include "dArray2DT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

using namespace Tahoe;

/* read list of ordered pairs */
void read_points(ifstreamT& in, dArray2DT& points, double& abs_max);

/* add zero function extensions to read given bounds */
void set_lower_bound(double x_lower, dArray2DT& points);
void set_upper_bound(double x_upper, dArray2DT& points);

/* integrated the squared difference between the two functions between the limits given */
double difference(const dArray2DT& pts1, const dArray2DT& pts2);
double area(const dArray2DT& pts);

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
		cout << "\n usage: " << caller << " [results 1 (ref)] [results 2] [output file] (x_min [min]) (x_max [max])\n" << endl;
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

	/* read function 2 */
	in.open(argv[2]);	
	if (!in.is_open())
		ExceptionT::GeneralFail(caller, "could not open file \"%s\"", in.filename());
		
	dArray2DT pts2;
	double max2;
	read_points(in, pts2, max2);
	in.close();

	/* integration bounds */
	double x1 = (pts1(0,0) > pts2(0,0)) ? pts1(0,0) : pts2(0,0);
	int n1 = pts1.MajorDim();
	int n2 = pts2.MajorDim();
	double x2 = (pts1(n1-1,0) < pts2(n2-1,0)) ? pts1(n1-1,0) : pts2(n2-1,0);

	/* user-defined bounds */
	if (argc > 4) {
		for (int i = 4; i < argc; i++)
		{
			if (strcmp(argv[i], "x_min") == 0) {
				x1 = atof(argv[i+1]);
				cout << "\n re-setting lower bound to " << x1 << endl;
				set_lower_bound(x1, pts1);
				set_lower_bound(x1, pts2);
			}
			else if (strcmp(argv[i], "x_max") == 0) {
				x2 = atof(argv[i+1]);
				cout << "\n re-setting upper bound to " << x2 << endl;
				set_upper_bound(x2, pts1);
				set_upper_bound(x2, pts2);
			}
		}
	}

	/* areas */
	double a1 = area(pts1);
	double a2 = area(pts2);

	/* results */
	cout << "   x_start: " << x1 << '\n';
	cout << "     x_end: " << x2 << '\n';
	cout << "        a1: " << a1 << '\n';
	cout << "        a2: " << a2 << '\n';
	cout << "norm error: " << (a2 - a1)/a1 << endl;
	
	/* write results */
	ofstreamT out(argv[3]);
	out.precision(12);
	int d_width = OutputWidth(out, &a1);
	out << setw(d_width) << (a2 - a1)/a1 << ' ' << "error_norm" << endl;
	
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

double difference(const dArray2DT& pts1, const dArray2DT& pts2)
{
	const char caller[] = "difference";
	double sum = 0.0;

	/* dimensions */
	int n1 = pts1.MajorDim();
	int n2 = pts2.MajorDim();

	/* expecting at least 2 points in each set */
	if (n1 < 2 || n2 < 2)
		ExceptionT::GeneralFail(caller, "at least 2 points needed in each set");
	
	/* smallest intervals */
	double dx1_min = pts1(n1-1,0) - pts1(0,0);
	double dx1_avg = dx1_min/(pts1.MajorDim() - 1);
	for (int i = 1; i < pts1.MajorDim(); i++)
	{
		double dx = pts1(i,0) - pts1(i-1,0);
		dx1_min = (dx < dx1_min) ? dx : dx1_min;
	}
	double dx2_min = pts2(n2-1,0) - pts2(0,0);
	double dx2_avg = dx2_min/(pts2.MajorDim() - 1);
	for (int i = 1; i < pts2.MajorDim(); i++)
	{
		double dx = pts2(i,0) - pts2(i-1,0);
		dx2_min = (dx < dx2_min) ? dx : dx2_min;
	}
	double small = (dx2_min < dx1_min) ? dx2_min : dx1_min;
	dx1_min /= dx1_avg;
	dx2_min /= dx2_avg;	
	double test = (dx2_min < dx1_min) ? dx2_min : dx1_min;
	if (test < 1.0e-12)
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
	
		/* integrate */
		sum += 0.5*(m2 - m1)*(2.0*x1*dx + dx*dx) + (b2 - b1)*dx;
		
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

double area(const dArray2DT& pts)
{
	const char caller[] = "difference";
	double sum = 0.0;
	for (int i = 1; i < pts.MajorDim(); i++) {
		double dx = pts(i,0) - pts(i-1,0);
		sum += 0.5*dx*(pts(i,1) + pts(i-1,1));
	}
	return sum;
}

/* add zero function extensions to read given bounds */
void set_lower_bound(double x_lower, dArray2DT& points)
{
	/* average interval */
	double dx_avg = (points(points.MajorDim()-1,0) - points(0,0))/(points.MajorDim()-1);

	/* less than an interval */
	if (fabs(points(0,0) - x_lower) < dx_avg) 
		return;
	else if (x_lower < points(0,0))
	{
		/* prepend extension */
		dArray2DT pts_tmp(points.MajorDim()+2, 2);
		pts_tmp.BlockRowCopyAt(points,2);

		pts_tmp(0,0) = x_lower;
		pts_tmp(0,1) = 0.0;

		pts_tmp(1,0) = points(0,0) - 1.0e-06*(points(0,0) - x_lower);
		pts_tmp(1,1) = 0.0;

		/* return */
		points.Swap(pts_tmp);
	}
	else /* chop beginning */
	{
		/* check */
		if (x_lower > points(points.MajorDim()-1,0))
			ExceptionT::GeneralFail("set_lower_bound",
				"lower bound %g exceeds upper range of data %g",
				x_lower, points(points.MajorDim()-1,0));

		/* find points to chop */
		int keep = 0;
		for (int i = 1; keep == i-1 && i < points.MajorDim(); i++)
			if (points(i,0) < x_lower)
				keep = i;
		
		/* resize */
		dArray2DT pts_tmp(points.MajorDim() - keep, 2);
		dArray2DT tmp(pts_tmp.MajorDim(), 2, points(keep));
		pts_tmp = tmp;
	
		/* linearly extrapolate to the bound */
		pts_tmp(0,0) = x_lower;
		if (fabs(points(keep,0) - x_lower) < kSmall)
			pts_tmp(0,1) = points(keep,1);
		else
		{
			double m = (points(keep+1,1) - points(keep,1))/(points(keep+1,0) - points(keep,0));
			double dx = x_lower - points(keep,0);
			pts_tmp(0,1) = points(keep,1) + dx*m;
		}
		
		/* return */
		pts_tmp.Swap(points);
	}
}

void set_upper_bound(double x_upper, dArray2DT& points)
{
	/* average interval */
	double dx_avg = (points(points.MajorDim()-1,0) - points(0,0))/(points.MajorDim()-1);

	/* nothing to do */
	if (fabs(points(points.MajorDim()-1,0) - x_upper) < dx_avg)
		return;
	else if (x_upper > points(points.MajorDim()-1,0))
	{
		/* append extension */
		dArray2DT pts_tmp(points.MajorDim()+2, 2);
		pts_tmp.BlockRowCopyAt(points,0);

		int index = points.MajorDim();
		pts_tmp(index,0) = points(index-1,0) + 1.0e-06*(x_upper - points(index-1,0));
		pts_tmp(index,1) = 0.0;

		index++;
		pts_tmp(index,0) = x_upper;
		pts_tmp(index,1) = 0.0;

		/* return */
		points.Swap(pts_tmp);
	}
	else /* chop tail */
	{
		/* check */
		if (x_upper < points(0,0))
			ExceptionT::GeneralFail("set_upper_bound",
				"upper bound %g exceeds lower range of data %g",
				x_upper, points(0,0));

		/* find points to chop */
		int keep = points.MajorDim()-1;
		for (int i = keep-1; keep == i+1 && i > -1; i--)
			if (points(i,0) > x_upper)
				keep = i;
		
		/* resize */
		dArray2DT pts_tmp(keep+1, 2);
		pts_tmp.BlockRowCopyAt(points, 0, pts_tmp.MajorDim());
	
		/* linearly extrapolate to the bound */
		pts_tmp(keep,0) = x_upper;
		if (fabs(points(keep,0) - x_upper) < kSmall)
			pts_tmp(keep,1) = points(keep,1);
		else
		{
			double m = (points(keep,1) - points(keep-1,1))/(points(keep,0) - points(keep-1,0));
			double dx = x_upper - points(keep,0);
			pts_tmp(keep,1) = points(keep,1) + dx*m;
		}
		
		/* return */
		pts_tmp.Swap(points);
	}
}
