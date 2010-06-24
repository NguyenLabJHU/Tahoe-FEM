/* $Id: extract_Disp.cpp,v 1.3 2010-06-24 14:08:33 tdnguye Exp $ */

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
	if (argc <4) {
		cout << "\n usage: " << caller << " [tahoe_output_file] [start_step] [end_step]\n\n";
		return 1;
	}

	/* initialize fixed node file */
	const char* input_file = argv[1];
	IOBaseT::FileTypeT input_file_type = IOBaseT::name_to_FileTypeT(input_file);
	InputBaseT* input = IOBaseT::NewInput(input_file_type, cout);
	input->Open(input_file);
	
	/*spatial dimensions*/
	int start_step = atoi(argv[2]);
	int end_step = atoi(argv[3]);
	/*read time*/
	int tot_num_steps = input->NumTimeSteps();
	dArrayT time_steps, times;
        time_steps.Allocate(tot_num_steps);
	input->ReadTimeSteps(time_steps);
	cout << "start_step = " << start_step << endl;
	cout << "end_step = " << end_step << endl;	

	/* locate "D_X" */
	ArrayT<StringT> disp_nlabels;
	input->ReadNodeLabels(disp_nlabels);
	iArrayT D_X_index;
	D_X_index.Dimension(3);
	
	ArrayT<StringT> D_X_label;
	D_X_label.Dimension(3);
	D_X_label[0] = "D_X";
	D_X_label[1] = "D_Y";
	D_X_label[2] = "D_Z";

	int nout = 0;
	for (int k = 0; k < 3; k++)
	{
		/* pick component */	
		D_X_index[k] = -1;
		for (int i = 0; D_X_index[k] == -1 && i < disp_nlabels.Length(); i++)
			if (disp_nlabels[i] == D_X_label[k])
			{
				D_X_index[k] = i;
				nout++;
			}
		if (D_X_index[k] > -1)
			cout << " displacement " << D_X_label[k] << " index = " << D_X_index[k] << endl;
	}
	/* workspace */
	const int nsd = input->NumDimensions();
	const int nnd = input->NumNodes();
	dArray2DT nodal_values(nnd,disp_nlabels.Length());
	dArray2DT coords(nnd, nsd);
	dArray2DT d_x(nodal_values.MajorDim(), nout);
	iArrayT node_id(nnd);
	d_x = 0.0;

	StringT file = argv[1];
	StringT path = file.FilePath();
	file = argv[1];
	if (path.StringLength()>0) { file.Drop(path.Length()-1);}
	input->ReadNodeID(node_id);
	for (int nstep = start_step; nstep<= end_step; nstep++)
	{
		/*initialize output file */
		char outfile[1000];
		
		strcpy(outfile,file); 
		char* suff;
			suff = strstr(outfile,".io");
		StringT temp="."; 	
		StringT suffnew = temp.Append(nstep, 4);
		strncpy(suff,suffnew,20);
		cout << "\noutfile: "<<outfile<<endl;
		ofstreamT out(outfile);

		double time =  time_steps[nstep];
		cout << "step: "<<nstep << ", time : " << time_steps[nstep] << "\n";

		/* read values */
		input->ReadAllNodeVariables(nstep, nodal_values);
		input->ReadCoordinates(coords);
	
		/* extract column */
		out << setprecision(12)<<d_x.MajorDim()<<endl;
		out << setprecision(12)<<d_x.MinorDim()<<endl;
		out << setprecision(12)<<time<<endl;
	
		for (int k = 0; k < nsd; k++)
			d_x.ColumnCopy(k, nodal_values, D_X_index[k]);
		for (int i = 0; i < d_x.MajorDim(); i++)
		{
			out << node_id[i];
			for (int j = 0; j < nsd; j++)
				out <<' '<<setprecision(12)<< coords(i,j);

			for (int j = 0; j < nout; j++)
				out <<' '<< setprecision(12)<<d_x(i,j);
			out <<endl;
		}
		out.flush();
	}

	return 0;
}
