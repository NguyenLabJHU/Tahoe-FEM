/* $Id: extract_CT.cpp,v 1.1 2004-11-12 21:09:57 paklein Exp $ */
#include "ModelManagerT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

using namespace Tahoe;

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
		cout << "\n usage: " << caller << " [pin KBC file] [pin fixed node file] [displacement file] [output file]\n" << endl;
		return 1;
	}

	/* initialize pin file */
	const char* pin_file = argv[1];
	IOBaseT::FileTypeT pin_file_type = IOBaseT::name_to_FileTypeT(pin_file);
	InputBaseT* pin_input = IOBaseT::NewInput(pin_file_type, cout);
	pin_input->Open(pin_file);
	
	/* locate "F_Y" */
	ArrayT<StringT> pin_nlabels;
	pin_input->ReadNodeLabels(pin_nlabels);
	int F_Y_index = -1;
	for (int i = 0; F_Y_index == -1 && i < pin_nlabels.Length(); i++)
		if (pin_nlabels[i] == "F_Y")
			F_Y_index = i;

	//TEMP	
	cout << "F_Y_index = " << F_Y_index << endl;

	/* initialize fixed node file */
	const char* pin_fixed_file = argv[2];
	IOBaseT::FileTypeT pin_fixed_file_type = IOBaseT::name_to_FileTypeT(pin_fixed_file);
	InputBaseT* pin_fixed_input = IOBaseT::NewInput(pin_fixed_file_type, cout);
	pin_fixed_input->Open(pin_fixed_file);

	/* locate "F_D_Y" */
	ArrayT<StringT> pin_fixed_nlabels;
	pin_fixed_input->ReadNodeLabels(pin_fixed_nlabels);
	int F_D_Y_index = -1;
	for (int i = 0; F_D_Y_index == -1 && i < pin_fixed_nlabels.Length(); i++)
		if (pin_fixed_nlabels[i] == "F_D_Y")
			F_D_Y_index = i;

	//TEMP	
	cout << "F_D_Y_index = " << F_D_Y_index << endl;

	/* initialize displacement file */
	const char* disp_file = argv[3];
	IOBaseT::FileTypeT disp_file_type = IOBaseT::name_to_FileTypeT(disp_file);
	InputBaseT* disp_input = IOBaseT::NewInput(disp_file_type, cout);
	disp_input->Open(disp_file);

	/* locate "D_Y" */
	ArrayT<StringT> disp_nlabels;
	disp_input->ReadNodeLabels(disp_nlabels);
	int D_Y_index = -1;
	for (int i = 0; D_Y_index == -1 && i < disp_nlabels.Length(); i++)
		if (disp_nlabels[i] == "D_Y")
			D_Y_index = i;

	//TEMP	
	cout << "D_Y_index = " << D_Y_index << endl;

	/* workspace */
	dArray2DT pin_values(pin_input->NumNodes(), pin_nlabels.Length());
	dArray2DT pin_fixed_values(pin_fixed_input->NumNodes(), pin_fixed_nlabels.Length());
	dArray2DT disp_values(disp_input->NumNodes(), disp_nlabels.Length());

	dArrayT f_y(pin_values.MajorDim());
	dArrayT f_d_y(pin_fixed_values.MajorDim());
	dArrayT d_y(disp_values.MajorDim());

	/* open output file */
	ofstreamT out(argv[4]);
	int num_steps = pin_input->NumTimeSteps();
	
	//TEMP
	cout << "num_steps = " << num_steps << endl;
	
	for (int i = 0; i < num_steps; i++)
	{
		/* read values */
		pin_input->ReadAllNodeVariables(i, pin_values);
		pin_fixed_input->ReadAllNodeVariables(i, pin_fixed_values);
		disp_input->ReadAllNodeVariables(i, disp_values);
	
		/* extract column */
		pin_values.ColumnCopy(F_Y_index, f_y);
		pin_fixed_values.ColumnCopy(F_D_Y_index, f_d_y);
		disp_values.ColumnCopy(D_Y_index, d_y);

		/* output */
		out << d_y.Max() << ' ' << f_y.Sum() + f_d_y.Sum() << '\n';
	}
	out.flush();

#if 0
	ifstreamT in;

	/* read function 1 */
	in.open(argv[1]);	
	if (!in.is_open())
		ExceptionT::GeneralFail(caller, "could not open file \"%s\"", in.filename());
		
	dArray2DT pts1;
	double max1;
	read_points(in, pts1, max1);
	in.close();
	
	cout << "function 1: " << max1 << '\n' << pts1 << endl;

	/* read function 2 */
	in.open(argv[2]);	
	if (!in.is_open())
		ExceptionT::GeneralFail(caller, "could not open file \"%s\"", in.filename());
		
	dArray2DT pts2;
	double max2;
	read_points(in, pts2, max2);
	in.close();
	
	cout << "function 2: " << max2 << '\n' << pts2 << endl;

	/* construct piecewise linear functions */
	PiecewiseLinearT func1(pts1);
	PiecewiseLinearT func2(pts2);
	
	/* integration bounds */
	double x1 = (pts1(0,0) > pts2(0,0)) ? pts1(0,0) : pts2(0,0);
	int n1 = pts1.MajorDim();
	int n2 = pts2.MajorDim();
	double x2 = (pts1(n1-1,0) < pts2(n2-1,0)) ? pts1(n1-1,0) : pts2(n2-1,0);
	
	/* trapezoidal scheme */
	double ref = (x2 - x1)*((max1 > max2) ? max1 : max2);
	int n_int = atoi(argv[3]);
	double e_sqr = squared_difference(func1, func2, x1, x2, n_int);

	/* results */
	cout << "   x_start: " << x1 << '\n';
	cout << "     x_end: " << x2 << '\n';
	cout << "     n_int: " << n_int << '\n';
	cout << "       ref: " << ref << '\n';
	cout << "   error^2: " << e_sqr << '\n';
	cout << "norm error: " << sqrt(e_sqr)/ref << endl;

	/* analytical scheme */
	e_sqr = squared_difference(pts1, pts2);

	/* results */
	cout << "   error^2: " << e_sqr << '\n';
	cout << "norm error: " << sqrt(e_sqr)/ref << endl;
	
#endif

	return 0;
}
