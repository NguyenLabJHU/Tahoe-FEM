/* $Id: calc_disp_cost.cpp,v 1.1 2010-01-28 22:21:09 tdnguye Exp $ */

#include "ModelManagerT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"

using namespace Tahoe;

/* reads in 2D array of data points to points and time stamp */
void read_points(ifstreamT& in, double& time, dArray2DT& points);

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
	if (argc <7) {
		cout << "\n usage: " << caller << " [tahoe_output_file] [data_file] [dakota_output_file] [start_step] [end_step] [tol]\n\n";
		cout << "\n data_file format: ";
		cout << "\n \t npoints:";
		cout << "\n \t nsd: ";
		cout << "\n \t time: ";
		cout << "\n \t int x0, y0, (z0), dx0, dy0, (dz0)";
		cout << "\n calculates the differences between displacement [dx0, dy0,dz0] of tahoe output file and data file,";
		cout << "\n for each time step from [start_step] to [end_step].  ";
		cout << "\n The program checks that the npoints, time, and coordinates [x0, y0, z0] match within the [tol]";
		cout << "\n The RMS difference is written to [dakota_output_file]";
		cout << "\n\n";
		return 1;
	}

	/*read tahoe_output_file*/
	/* initialize fixed node file */
	const char* input_file = argv[1];
	IOBaseT::FileTypeT input_file_type = IOBaseT::name_to_FileTypeT(input_file);
	InputBaseT* input = IOBaseT::NewInput(input_file_type, cout);
	input->Open(input_file);
	
	/*spatial dimensions*/
	int start_step = atoi(argv[4]);
	int end_step = atoi(argv[5]);
	/*read time*/
	int tot_num_steps = input->NumTimeSteps();
	dArrayT time_steps, times;
        time_steps.Allocate(tot_num_steps);
	input->ReadTimeSteps(time_steps);
	cout << "start_step = " << start_step << endl;
	cout << "end_step = " << end_step << endl;	
	
	double ctol = atof(argv[6]);
	cout << "tolerance = " << ctol << endl;	
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
	
	/*calculates cost function*/
	double cost = 0.0;
	for (int nstep = start_step; nstep<= end_step; nstep++)
	{
		/*read data file*/
		ifstreamT in;
		StringT datafile;
		strcpy(datafile, argv[2]);
		StringT data_temp = ".";
		StringT data_suff = data_temp.Append(nstep,4);
		strcat(datafile,data_suff);
		cout << "\nopening datafile: "<<datafile << endl;
		
		in.open(datafile);
		if (!in.is_open())
			ExceptionT::GeneralFail(caller, "could not open file \"%s\"", in.filename());
		dArray2DT data_pts;
		double data_time;
		read_points(in, data_time, data_pts);

		in.close();
		
		//////////////////////////////////////////////////
		/*reads input file*/

		double time =  time_steps[nstep];
		cout << "step: "<<nstep << ", time : " << time_steps[nstep] << "\n";
		if (fabs(time - data_time) > ctol)
			ExceptionT::GeneralFail(caller, "time stamp of data file and tahoe output file does not match, %f", data_time);
		if (nnd != data_pts.MajorDim())
			ExceptionT::GeneralFail(caller, "number of nodes  in data file and tahoe output file does not match");
		if (2*nsd != data_pts.MinorDim())
			ExceptionT::GeneralFail(caller, "number of spatial dimensions of data file and tahoe output file does not match");

		dArray2DT calc_pts;
		calc_pts.Dimension(nnd,nsd);
		
		/* read values */
		input->ReadAllNodeVariables(nstep, nodal_values);
		input->ReadCoordinates(coords);
	
		for (int k = 0; k < nsd; k++)
			d_x.ColumnCopy(k, nodal_values, D_X_index[k]);

		double error=0.0;
		for (int i; i<nnd; i++)
		{
			/*check for points that are congruent*/
			double distance=0.0;
			for (int j; j<nsd; j++)
			{
				distance += (data_pts(i,j) - coords(i,j))*(data_pts(i,j) - coords(i,j));
				error += (data_pts(i,j+nsd) - d_x(i,j))*(data_pts(i,j+nsd) - d_x(i,j));
			}
			if (distance > ctol)
				ExceptionT::GeneralFail(caller, "coords of point %d do not match", i);
		}
		cost += sqrt(error/nnd);
	}
		
/*
	char outfile[1000];
	strcpy(outfile,file); 
	char* suff;
	suff = strstr(outfile,".io");
	StringT temp=".out"; 	
	strncpy(suff,temp,20);
*/
	StringT outfile;
	strcpy(outfile, argv[3]);

	cout << "\noutfile: "<<outfile<<endl;
	ofstreamT out(outfile);
	out << cost/tot_num_steps;
	out<<endl;
	out.flush();
	return 0;
}

void read_points(ifstreamT& in,  double& data_time, dArray2DT& points)
{
	int num_points;
	int nsd;

	in >> num_points;
	in >> nsd;
	in >> data_time;
	
	/* file format*/
	/*n x0, y0, (z0), dx0, dy0, (dz0)*/ 
	/*dimension data array npoints x 2*nsd*/
	points.Dimension(num_points, 2.0*nsd);
	points = 0.0;
	int nn;
	double p1,p2,p3,p4,p5,p6;
	
	for (int i = 0; i < num_points; i++)
	{
		if (nsd ==2)
		{
			in >> nn >> p1 >> p2 >>p3>>p4;
			points(i,0) = p1;
			points(i,1) = p2;
			points(i,2) = p3;
			points(i,3) = p4;
		}
		else if (nsd ==3)
		{
			in >> nn >> p1 >> p2 >>p3>>p4;
			points(i,0) = p1;
			points(i,1) = p2;
			points(i,2) = p3;
			points(i,3) = p4;
			points(i,4) = p5;
			points(i,5) = p6;
		}
	}
	/*restore error flags*/
	in.clear();
//	cout << "data_time: "<<data_time<<endl;
//	cout << "data_points: "<<points<<endl;
	
}