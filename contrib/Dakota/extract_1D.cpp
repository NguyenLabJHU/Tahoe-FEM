/* $Id: extract_1D.cpp,v 1.2 2004-12-06 16:30:50 rjones Exp $ */
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

	/* pick component */
	char dir[3];
        strcpy(dir,"_Z");
	if     (strstr(argv[3],"_X"))  strcpy(dir,"_X");
	else if(strstr(argv[3],"_Y"))  strcpy(dir,"_X"); 
	else if(strstr(argv[3],"_Z"))  strcpy(dir,"_Z");
	cout << " dir "  << dir << endl;

	/* locate "F_D_X" */
	ArrayT<StringT> pin_fixed_nlabels;
	pin_fixed_input->ReadNodeLabels(pin_fixed_nlabels);
	char F_D_X_label[6];
	strcpy(F_D_X_label,"F_D"); strcat(F_D_X_label,dir);
	int F_D_X_index = -1;
	for (int i = 0; F_D_X_index == -1 && i < pin_fixed_nlabels.Length(); i++)
		if (pin_fixed_nlabels[i] ==  F_D_X_label)
			F_D_X_index = i;

	//TEMP	
	cout << " force " <<  F_D_X_label << " index = " << F_D_X_index << endl;

	/* initialize displacement file */
	const char* disp_file = argv[2];
	IOBaseT::FileTypeT disp_file_type = IOBaseT::name_to_FileTypeT(disp_file);
	InputBaseT* disp_input = IOBaseT::NewInput(disp_file_type, cout);
	disp_input->Open(disp_file);

	/* locate "D_X" */
	ArrayT<StringT> disp_nlabels;
	disp_input->ReadNodeLabels(disp_nlabels);
	char D_X_label[4];
	strcpy(D_X_label,"D"); strcat(D_X_label,dir);
	int D_X_index = -1;
	for (int i = 0; D_X_index == -1 && i < disp_nlabels.Length(); i++)
		if (disp_nlabels[i] == D_X_label)
			D_X_index = i;

	//TEMP	
	cout << " displacement " << D_X_label << " index = " << D_X_index << endl;

	/* workspace */
	dArray2DT pin_fixed_values(pin_fixed_input->NumNodes(), pin_fixed_nlabels.Length());
	dArray2DT disp_values(disp_input->NumNodes(), disp_nlabels.Length());

	dArrayT f_d_x(pin_fixed_values.MajorDim());
	dArrayT d_x(disp_values.MajorDim());

	/* open output file */
	ofstreamT out(argv[3]);
	int num_steps = pin_fixed_input->NumTimeSteps();
	dArrayT time_steps;
        time_steps.Allocate(num_steps);
	pin_fixed_input->ReadTimeSteps(time_steps);
	
	//TEMP
	cout << "num_steps = " << num_steps << endl;

	bool time = 0;
	if(strstr(argv[3],".fvst")) { time = 1;  cout << "f vs t\n";}
	else if(strstr(argv[3],".fvsd")) cout << "f vs d\n";
	
	for (int i = 0; i < num_steps; i++)
	{
	//	cout << i << ", time : " << time_steps[i] << "\n";
		double time =  time_steps[i];

		/* read values */
		pin_fixed_input->ReadAllNodeVariables(i, pin_fixed_values);
		disp_input->ReadAllNodeVariables(i, disp_values);
	
		/* extract column */
		pin_fixed_values.ColumnCopy(F_D_X_index, f_d_x);
		disp_values.ColumnCopy(D_X_index, d_x);

		/* output */
		if (time) out << time  << ' ' << f_d_x.Sum() << '\n';
	  	else      out << d_x.Max() << ' ' << f_d_x.Sum() << '\n';
	}
	out.flush();

	return 0;
}
