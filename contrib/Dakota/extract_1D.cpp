/* $Id: extract_1D.cpp,v 1.1 2004-11-12 21:09:57 paklein Exp $ */
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
	if (argc < 4) {
		cout << "\n usage: " << caller << " [force file] [displacement file] [output file]\n" << endl;
		return 1;
	}

	/* initialize fixed node file */
	const char* pin_fixed_file = argv[1];
	IOBaseT::FileTypeT pin_fixed_file_type = IOBaseT::name_to_FileTypeT(pin_fixed_file);
	InputBaseT* pin_fixed_input = IOBaseT::NewInput(pin_fixed_file_type, cout);
	pin_fixed_input->Open(pin_fixed_file);

	/* locate "F_D_Z" */
	ArrayT<StringT> pin_fixed_nlabels;
	pin_fixed_input->ReadNodeLabels(pin_fixed_nlabels);
	int F_D_Z_index = -1;
	for (int i = 0; F_D_Z_index == -1 && i < pin_fixed_nlabels.Length(); i++)
		if (pin_fixed_nlabels[i] == "F_D_Z")
			F_D_Z_index = i;

	//TEMP	
	cout << "F_D_Z_index = " << F_D_Z_index << endl;

	/* initialize displacement file */
	const char* disp_file = argv[2];
	IOBaseT::FileTypeT disp_file_type = IOBaseT::name_to_FileTypeT(disp_file);
	InputBaseT* disp_input = IOBaseT::NewInput(disp_file_type, cout);
	disp_input->Open(disp_file);

	/* locate "D_Z" */
	ArrayT<StringT> disp_nlabels;
	disp_input->ReadNodeLabels(disp_nlabels);
	int D_Z_index = -1;
	for (int i = 0; D_Z_index == -1 && i < disp_nlabels.Length(); i++)
		if (disp_nlabels[i] == "D_Z")
			D_Z_index = i;

	//TEMP	
	cout << "D_Z_index = " << D_Z_index << endl;

	/* workspace */
	dArray2DT pin_fixed_values(pin_fixed_input->NumNodes(), pin_fixed_nlabels.Length());
	dArray2DT disp_values(disp_input->NumNodes(), disp_nlabels.Length());

	dArrayT f_d_y(pin_fixed_values.MajorDim());
	dArrayT d_y(disp_values.MajorDim());

	/* open output file */
	ofstreamT out(argv[3]);
	int num_steps = pin_fixed_input->NumTimeSteps();
	
	//TEMP
	cout << "num_steps = " << num_steps << endl;
	
	for (int i = 0; i < num_steps; i++)
	{
		/* read values */
		pin_fixed_input->ReadAllNodeVariables(i, pin_fixed_values);
		disp_input->ReadAllNodeVariables(i, disp_values);
	
		/* extract column */
		pin_fixed_values.ColumnCopy(F_D_Z_index, f_d_y);
		disp_values.ColumnCopy(D_Z_index, d_y);

		/* output */
		out << d_y.Max() << ' ' << f_d_y.Sum() << '\n';
	}
	out.flush();

	return 0;
}
