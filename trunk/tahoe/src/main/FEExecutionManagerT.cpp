/* $Id: FEExecutionManagerT.cpp,v 1.68 2004-07-22 08:32:54 paklein Exp $ */
/* created: paklein (09/21/1997) */
#include "FEExecutionManagerT.h"

#include <iostream.h>
#include <iomanip.h>
#include <time.h>
#include <ctype.h>
#include <stdlib.h>

#if defined(__MWERKS__) && __option(profile)
#include <Profiler.h>
#endif

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "FEManagerT.h"
#include "FEManagerT_mpi.h"
#include "IOManager_mpi.h"
#include "ModelManagerT.h"
#include "CommManagerT.h"
#include "StringT.h"
#include "GraphT.h"
#include "PartitionT.h"
#include "OutputSetT.h"
#include "ModelFileT.h"
#include "ExodusT.h"
#include "JoinOutputT.h"
#include "dArrayT.h"
#include "OutputBaseT.h"
#include "CommunicatorT.h"

/* parameters */
#include "ParameterListT.h"
#include "ParameterTreeT.h"
#include "XML_Attribute_FormatterT.h"
#include "DotLine_FormatterT.h"

/* needed for bridging calculations FEExecutionManagerT::RunBridging */
#ifdef BRIDGING_ELEMENT
#include "FEManagerT_bridging.h"
#include "MultiManagerT.h"
#ifdef __DEVELOPMENT__
#include "FEManagerT_THK.h"
#endif
#include "TimeManagerT.h"
#include "NodeManagerT.h"
#include "dSPMatrixT.h"
#include "FieldT.h"
#include "IntegratorT.h"
#include "ElementBaseT.h"
#include "EAMFCC3D.h"
#endif

using namespace Tahoe;

/* Constructor */
FEExecutionManagerT::FEExecutionManagerT(int argc, char* argv[], char job_char,
	char batch_char, CommunicatorT& comm):
	ExecutionManagerT(argc, argv, job_char, batch_char, comm, 0)
{
#ifdef __CPLANT__
	/* if not prescribed as joined, write separate files */
	if (!CommandLineOption("-join_io")) AddCommandLineOption("-split_io");
#endif

	/* set communicator log level */
	if (CommandLineOption("-verbose")) comm.SetLogLevel(CommunicatorT::kLow);

	/* check for -dtd */
	if (CommandLineOption("-dtd")) AddCommandLineOption("-dtd");

	/* check for -xsd */
	if (CommandLineOption("-xsd")) AddCommandLineOption("-xsd");
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* MUST be overloaded */
void FEExecutionManagerT::RunJob(ifstreamT& in, ostream& status)
{
	const char caller[] = "FEExecutionManagerT::RunJob";
	
	/* set the path to the root file */
	StringT root;
	root.FilePath(in.filename());
	fstreamT::SetRoot(root);

	/* mode - job by default */
	ModeT mode = kJob;

	/* check command line options */
	for (int i = 0; i < fCommandLineOptions.Length(); i++)
	{
		if (fCommandLineOptions[i] == "-decomp")
			mode = kDecompose;
		else if (fCommandLineOptions[i] == "-join")
			mode = kJoin;
	}

	/* check second char */
	if (mode == kJob) {
		
		/* peek at next char */
		char a = in.next_char();
		if (a == fJobChar) {
			mode = kBridging;
			in >> a; /* clear character */	
		}
	}

	//TEMP - look for command line option indicating THK calculation
	if (CommandLineOption("-thk")) mode = kTHK;

	switch (mode)
	{
		case kJob:
		{
			if (fComm.Size() == 1) {
				cout << "\n RunJob_serial: " << in.filename() << endl;
				RunJob_serial(in.filename(), status);
			} else {
				cout << "\n RunJob_parallel: " << in.filename() << endl;
				RunJob_parallel(in.filename(), status);
			}
			break;
		}
		case kDecompose:
		{
			cout << "\n RunDecomp_serial: " << in.filename() << endl;
			if (fComm.Size() > 1) cout << " RunDecomp_serial: SERIAL ONLY" << endl;

			/* 'serial' communicator */
			int rank = fComm.Rank();
			CommunicatorT comm(fComm, (rank == 0) ? rank : CommunicatorT::kNoColor);
			
			/* decompose on rank 0 */
			if (rank == 0) RunDecomp_serial(in.filename(), status, comm);

			/* synch */
			fComm.Barrier();
			break;
		}
		case kJoin:
		{
			cout << "\n RunJoin_serial: " << in.filename() << endl;
			if (fComm.Size() > 1) cout << " RunJoin_serial: SERIAL ONLY" << endl;
			
			/* 'serial' communicator */
			int rank = fComm.Rank();
			CommunicatorT comm(fComm, (rank == 0) ? rank : CommunicatorT::kNoColor);

			/* join using rank 0 */
			if (rank == 0) RunJoin_serial(in.filename(), status);

			/* synch */
			fComm.Barrier();
			break;
		}
#pragma message("delete me")
#if 0
		case kBridging:
		{
			cout << "\n RunBridging: " << in.filename() << endl;
			if (fComm.Size() > 1) ExceptionT::GeneralFail(caller, "RunBridging for SERIAL ONLY");

			RunBridging(in, status);
			break;
		}
#endif
		case kTHK:
		{
#ifdef __DEVELOPMENT__
			cout << "\n RunTHK: " << in.filename() << endl;
			if (fComm.Size() > 1) ExceptionT::GeneralFail(caller, "RunTHK for SERIAL ONLY");

			RunTHK(in, status);
			break;
#else
			ExceptionT::BadInputValue(caller, "development module needed to run THK");
#endif
		}
		default:
			ExceptionT::GeneralFail("FEExecutionManagerT::RunJob", "unknown mode: %d", mode);
	}

	/* clear the path to the root file */
	fstreamT::SetRoot(NULL);
}

/**********************************************************************
 * Protected
 **********************************************************************/

bool FEExecutionManagerT::AddCommandLineOption(const char* str)
{
	/* remove mutually exclusive command line options */
	StringT option(str);
	if (option == "-run")
	{
		RemoveCommandLineOption("-decomp");		
		RemoveCommandLineOption("-join");
	}
	else if (option == "-decomp")
	{
		RemoveCommandLineOption("-run");
		RemoveCommandLineOption("-join");
	}
	else if (option == "-join")
	{
		RemoveCommandLineOption("-decomp");		
		RemoveCommandLineOption("-run");
	}
	else if (option == "-join_io")
		RemoveCommandLineOption("-split_io");
	else if (option == "-split_io")
		RemoveCommandLineOption("-join_io");
	else if (option == "-decomp_method") {
	
		/* clear existing setting */
		int index;
		if (CommandLineOption("-decomp_method", index))
		{
			/* clear option */
			fCommandLineOptions.DeleteAt(index);

			/* clear number */
			const StringT& opt = fCommandLineOptions[index];
			if (strlen(opt) > 1 && isdigit(opt[1]))
				fCommandLineOptions.DeleteAt(index);
		}
	}
	else if (option == "-dtd") {
		RunWriteDescription(XML_Attribute_FormatterT::DTD);
		return true;
	}
	else if (option == "-xsd") {
		RunWriteDescription(XML_Attribute_FormatterT::XSD);
		return true;
	}
	
	/* inherited */
	return ExecutionManagerT::AddCommandLineOption(option);
}

/* remove the command line option to the list */
bool FEExecutionManagerT::RemoveCommandLineOption(const char* str)
{
	StringT option(str);
	if (option == "-decomp")
	{
		int index;
		if (CommandLineOption("-decomp", index))
		{
			/* check for following number option */
			if (fCommandLineOptions.Length() > index+1) {
				const StringT& opt = fCommandLineOptions[index+1];
				if (strlen(opt) > 1 && isdigit(opt[1]))
					RemoveCommandLineOption(opt);
			}
			
			/* remove decomp */
			return ExecutionManagerT::RemoveCommandLineOption(option);
		}
		return false;
	}
	else if (option == "-join")
	{
		int index;
		if (CommandLineOption("-join", index))
		{
			/* check for following number option */
			if (fCommandLineOptions.Length() > index+1) {
				const StringT& opt = fCommandLineOptions[index+1];
				if (strlen(opt) > 1 && isdigit(opt[1]))
					RemoveCommandLineOption(opt);
			}
			
			/* remove decomp */
			return ExecutionManagerT::RemoveCommandLineOption(option);		
		}
		return false;
	} 
	else /* inherited */
		return ExecutionManagerT::RemoveCommandLineOption(option);
}

/**********************************************************************
 * Private
 **********************************************************************/

#ifdef BRIDGING_ELEMENT
#ifdef __DEVELOPMENT__
void FEExecutionManagerT::RunDynamicBridging(FEManagerT_bridging& continuum, FEManagerT_THK& atoms, ofstream& log_out, ifstreamT& in) const
{
	const char caller[] = "FEExecutionManagerT::RunDynamicBridging";
	/* configure ghost nodes */
	ModelManagerT* model = continuum.ModelManager();
	int nsd = model->NumDimensions();		
	int group = 0;
	int order1 = 0;	// For InterpolateField, 3 calls to obtain displacement/velocity/acceleration
	int order2 = 1;
	int order3 = 2;
	double dissipation = 0.0;
	dArray2DT field_at_ghosts, totalu, fubig, fu, projectedu, boundghostdisp, boundghostvel, boundghostacc;
	dArray2DT thkforce, totaldisp, elecdens, embforce, wavedisp;
	dSPMatrixT ntf;
	iArrayT activefenodes, boundaryghostatoms;
	StringT bridging_field = "displacement";
//	atoms.InitGhostNodes(continuum.ProjectImagePoints());
#pragma message("fix me")

	bool makeinactive = false;	
	/* figure out boundary atoms for use with THK boundary conditions, 
	   ghost atoms for usage with MD force calculations */
	if (nsd == 2)
		boundaryghostatoms = atoms.InterpolationNodes2D();
	else
		boundaryghostatoms = atoms.InterpolationNodes3D();
		
	int numgatoms = (atoms.GhostNodes()).Length();	// total number of ghost atoms
	int numbatoms = boundaryghostatoms.Length() - numgatoms;	// total number of boundary atoms
	dArray2DT gadisp(numgatoms,nsd), gavel(numgatoms,nsd), gaacc(numgatoms,nsd);
	dArray2DT badisp(numbatoms,nsd), bavel(numbatoms,nsd), baacc(numbatoms,nsd);
	iArrayT allatoms(boundaryghostatoms.Length()), gatoms(numgatoms), batoms(numbatoms), boundatoms(numbatoms);
	allatoms.SetValueToPosition();
	batoms.CopyPart(0, allatoms, numgatoms, numbatoms);
	gatoms.CopyPart(0, allatoms, 0, numgatoms);      
	boundatoms.CopyPart(0, boundaryghostatoms, numgatoms, numbatoms);
	//elecdens.Dimension(gatoms.Length(), 1);
	//embforce.Dimension(gatoms.Length(), 1);
	continuum.InitInterpolation(bridging_field, boundaryghostatoms, *atoms.NodeManager());
	//dArrayT mdmass;
	//atoms.LumpedMass(atoms.NonGhostNodes(), mdmass);	// acquire array of MD masses to pass into InitProjection, etc...
	continuum.InitProjection(bridging_field, *atoms.CommManager(), atoms.NonGhostNodes(), *atoms.NodeManager(), makeinactive);		

	/* nodes to include/exclude in calculation of the atomistic displacements. Dimension
	 * of these matricies is the number atom types */
	nMatrixT<int> ghostonmap(2), ghostoffmap(2);  // 3D wave propagation
	//nMatrixT<int> ghostonmap(5), ghostoffmap(5);  // for 2D fracture problem
	//nMatrixT<int> ghostonmap(4), ghostoffmap(4);    // for planar wave propagation problem
	//nMatrixT<int> ghostonmap(7), ghostoffmap(7);	// 3D fracture problem
	ghostonmap = 0;
	ghostoffmap = 0;
	ghostoffmap(0,1) = ghostoffmap(1,0) = 1;	// 3D wave propagation
	//ghostonmap(4,5) = ghostonmap(5,4) = 1;	// 3D fracture problem
	//ghostoffmap(0,6) = ghostoffmap(6,0) = ghostoffmap(1,6) = ghostoffmap(6,1) = 1;
	//ghostoffmap(2,6) = ghostoffmap(6,2) = ghostoffmap(3,6) = ghostoffmap(6,3) = 1;
	//ghostoffmap(4,6) = ghostoffmap(6,4) = ghostoffmap(5,6) = ghostoffmap(6,5) = 1;
	//ghostoffmap(1,0) = ghostoffmap(0,1) = 1;  // for wave propagation problem
	//ghostoffmap(4,0) = ghostoffmap(0,4) = ghostoffmap(4,1) = ghostoffmap(1,4) = 1;  // center MD crack
	//ghostoffmap(4,2) = ghostoffmap(2,4) = ghostoffmap(2,3) = ghostoffmap(3,2) = 1;
	//ghostoffmap(4,3) = ghostoffmap(3,4) = 1;
	//ghostonmap(2,3) = ghostonmap(3,2) = 1;
	//ghostoffmap(1,0) = ghostoffmap(0,1) = ghostoffmap(3,0) = ghostoffmap(0,3) = 1; // left edge MD crack
	//ghostoffmap(1,3) = ghostoffmap(3,1) = ghostoffmap(2,3) = ghostoffmap(3,2) = 1;
	//ghostonmap(1,0) = ghostonmap(0,1) = 1;

	if (nsd == 3)
	{
		/* set pointers to embedding force/electron density in FEManagerT_bridging atoms */
		//atoms.SetExternalElecDensity(elecdens, atoms.GhostNodes());
		//atoms.SetExternalEmbedForce(embforce, atoms.GhostNodes());
		//continuum.ElecDensity(gatoms.Length(), elecdens, embforce);
	}
	
	/* time managers */
	TimeManagerT* atom_time = atoms.TimeManager();
	TimeManagerT* continuum_time = continuum.TimeManager();

//	atom_time->Top();
//	continuum_time->Top();
	int d_width = OutputWidth(log_out, field_at_ghosts.Pointer());
	
//	while (atom_time->NextSequence() && continuum_time->NextSequence())
//	{		
		/* set to initial condition */
		atoms.InitialCondition();
		
		/* calculate fine scale part of MD displacement and total displacement u */
		continuum.InitialProject(bridging_field, *atoms.NodeManager(), projectedu, order1);
	
		/* solve for initial FEM force f(u) as function of fine scale + FEM */
		/* use projected totalu instead of totalu for initial FEM displacements */
		nMatrixT<int>& promap = atoms.PropertiesMap(0);   // element group for particles = 0
		promap = ghostoffmap;  // turn ghost atoms off for f(u) calculation
		fubig = InternalForce(projectedu, atoms);
		promap = ghostonmap;   // turn ghost atoms back on for MD force calculations
		
		/* calculate global interpolation matrix ntf */
		continuum.Ntf(ntf, atoms.NonGhostNodes(), activefenodes);
		
		/* compute FEM RHS force as product of ntf and fu */
		dArrayT fx(ntf.Rows()), fy(ntf.Rows()), fz(ntf.Rows()), tempx(ntf.Cols()), tempy(ntf.Cols()), tempz(ntf.Cols());
		dArray2DT ntfproduct(ntf.Rows(), nsd);
		fu.Dimension((atoms.NonGhostNodes()).Length(), nsd);
		fu.RowCollect(atoms.NonGhostNodes(), fubig);   
		fu.ColumnCopy(0,tempx);
		ntf.Multx(tempx,fx);
		ntfproduct.SetColumn(0,fx);
		fu.ColumnCopy(1,tempy);
		ntf.Multx(tempy,fy);
		ntfproduct.SetColumn(1,fy);
		if (nsd == 3)
		{
			fu.ColumnCopy(2,tempz);
			ntf.Multx(tempz,fz);
			ntfproduct.SetColumn(2,fz);
		}

		/* Add FEM RHS force to RHS using SetExternalForce */
		continuum.SetExternalForce(bridging_field, ntfproduct, activefenodes);
	
		/* now d0, v0 and a0 are known after InitialCondition */
		continuum.InitialCondition();
		
		if (nsd == 3)
		{
			/* Calculate EAM electron density/embedding terms for ghost atoms using continuum information */
			//continuum.ElecDensity(gatoms.Length(), elecdens, embforce);
		}
		
		/* Interpolate FEM values to MD ghost nodes which will act as MD boundary conditions */
		continuum.InterpolateField(bridging_field, order1, boundghostdisp);
		continuum.InterpolateField(bridging_field, order2, boundghostvel);
		continuum.InterpolateField(bridging_field, order3, boundghostacc);
		
		/* sort boundary + ghost atom info into separate arrays */
		gadisp.RowCollect(gatoms, boundghostdisp);
		gavel.RowCollect(gatoms, boundghostvel);
		gaacc.RowCollect(gatoms, boundghostacc);
		badisp.RowCollect(batoms, boundghostdisp);
		bavel.RowCollect(batoms, boundghostvel);
		baacc.RowCollect(batoms, boundghostacc);

		/* Removed atoms.SetFieldValues(bridging_field, atoms.GhostNodes(), order1, gadisp); here */	
			
		if (nsd == 2)
		{
			/* store initial MD boundary displacement histories */
			thkforce = atoms.THKForce(badisp);
			atoms.SetExternalForce(bridging_field, thkforce, boundatoms);  // sets pointer to thkforce 
		}
		else
		{
			/* thkdisp = fine scale part of ghost atom displacement */
			thkforce = atoms.THKDisp(badisp);
			atoms.SetExternalForce(bridging_field, thkforce, boundatoms);
			//totaldisp+=gadisp;	// add FEM coarse scale part of displacement
			//atoms.SetFieldValues(bridging_field, atoms.GhostNodes(), order1, totaldisp);
		}
		
		/* figure out timestep ratio between fem and md simulations */
		int nfesteps = continuum_time->NumberOfSteps();
		double mddt = atom_time->TimeStep();
		double fedt = continuum_time->TimeStep();
		double d_ratio = fedt/mddt;		
		int ratio = int((2.0*d_ratio + 1.0)/2.0);
		
		/* running status flag */
		ExceptionT::CodeT error = ExceptionT::kNoError;	

		for (int i = 0; i < nfesteps; i++)	
		{
			for (int j = 0; j < ratio; j++)	// MD update first
			{
				atom_time->Step();	
								
				/* initialize step */
				if (1 || error == ExceptionT::kNoError) error = atoms.InitStep();
				
				/* update FEM solution interpolated at boundary atoms and ghost atoms assuming 
				constant acceleration - because of constant acceleration assumption, predictor and 
				corrector are combined into one function */
				atoms.BAPredictAndCorrect(mddt, badisp, bavel, baacc);
				atoms.BAPredictAndCorrect(mddt, gadisp, gavel, gaacc);

				if (nsd == 2)
				{
					/* Write interpolated FEM values at MD ghost nodes into MD field - displacements only */
					atoms.SetFieldValues(bridging_field, atoms.GhostNodes(), order1, gadisp);	
				
					/* calculate THK force on boundary atoms, update displacement histories */
					thkforce = atoms.THKForce(badisp);  // SetExternalForce set via pointer
				}
				else
				{
					/* Write interpolated FEM values at MD ghost nodes into MD field - displacements only */
					atoms.SetFieldValues(bridging_field, atoms.GhostNodes(), order1, gadisp);
					
					/* calculate thkforces */
					thkforce = atoms.THKDisp(badisp);
				}
				
				/* solve MD equations of motion */
				if (1 || error == ExceptionT::kNoError) {
						atoms.ResetCumulativeUpdate(group);
						error = atoms.SolveStep();
				}

				/* close  md step */
				if (1 || error == ExceptionT::kNoError) error = atoms.CloseStep();    
			}
			continuum_time->Step();

			/* initialize step */
			if (1 || error == ExceptionT::kNoError) error = continuum.InitStep();
            
			/* calculate total displacement u = FE + fine scale here using updated FEM displacement */
			continuum.BridgingFields(bridging_field, *atoms.NodeManager(), *continuum.NodeManager(), totalu);
		
			/* calculate FE internal force as function of total displacement u here */
			promap = ghostoffmap;  // turn off ghost atoms for f(u) calculations
			fubig = InternalForce(totalu, atoms);
			promap = ghostonmap;   // turn on ghost atoms for MD force calculations
			fu.RowCollect(atoms.NonGhostNodes(), fubig); 
			fu.ColumnCopy(0,tempx);
			ntf.Multx(tempx,fx);
			ntfproduct.SetColumn(0,fx);
			fu.ColumnCopy(1,tempy);
			ntf.Multx(tempy,fy);
			ntfproduct.SetColumn(1,fy);	// SetExternalForce updated via pointer
			if (nsd == 3)
			{
				fu.ColumnCopy(2,tempz);
				ntf.Multx(tempz,fz);
				ntfproduct.SetColumn(2,fz);
			}
			
			/* no need to call SetExternalForce again due to pointers to ntfproduct */
			
			/* solve FE equation of motion using internal force just calculated */
			if (1 || error == ExceptionT::kNoError) {
					continuum.ResetCumulativeUpdate(group);
					error = continuum.SolveStep();
			}

			/* Interpolate FEM values to MD ghost nodes which will act as MD boundary conditions */
			continuum.InterpolateField(bridging_field, order1, boundghostdisp);
			continuum.InterpolateField(bridging_field, order2, boundghostvel);
			continuum.InterpolateField(bridging_field, order3, boundghostacc);
			
			/* sort boundary + ghost atom info into separate arrays */
			gadisp.RowCollect(gatoms, boundghostdisp);
			gavel.RowCollect(gatoms, boundghostvel);
			gaacc.RowCollect(gatoms, boundghostacc);
			badisp.RowCollect(batoms, boundghostdisp);
			bavel.RowCollect(batoms, boundghostvel);
			baacc.RowCollect(batoms, boundghostacc);
			
			if (nsd == 3)
			{
				/* Calculate EAM electron density/embedding terms for ghost atoms using updated continuum information */
				//continuum.ElecDensity(gatoms.Length(), elecdens, embforce);
			}
			
			/* close fe step */
			if (1 || error == ExceptionT::kNoError) error = continuum.CloseStep();
                        
		}

		/* check for error */
		if (0)
			ExceptionT::GeneralFail(caller, "hit error %d", error);
                
//	}
}
#endif

/* calculate MD internal force as function of total bridging scale displacement u */
const dArray2DT& FEExecutionManagerT::InternalForce(dArray2DT& totalu, FEManagerT_bridging& atoms) const
{
	/* first obtain the MD displacement field */
	StringT bridging_field = "displacement";
	FieldT* atomfield = atoms.NodeManager()->Field(bridging_field);
	dArray2DT mddisp = (*atomfield)[0];	// temporarily store permanent MD displacements
	int group = 0;	// assume particle group number = 0
	
	/* obtain atom node list - can calculate once and store... */
	int nnd = totalu.MajorDim();
	iArrayT nodes(nnd);
	nodes.SetValueToPosition();
		
	/* now write total bridging scale displacement u into field */
	int order = 0;	// write displacement only
	atoms.SetFieldValues(bridging_field, nodes, order, totalu);
		
	/* compute RHS - ParticlePairT fForce calculated by this call */
	atoms.FormRHS(group);
	
	/* write actual MD displacements back into field */
	atoms.SetFieldValues(bridging_field, nodes, order, mddisp);

	/* get the internal force contribution associated with the last call to FormRHS */
	return atoms.InternalForce(group);
}

#else /* bridging element not enabled */

void FEExecutionManagerT::RunBridging(ifstreamT& in, ostream& status) const
{
#pragma unused(in)
#pragma unused(status)

	const char caller[] = "FEExecutionManagerT::RunBridging";
	ExceptionT::GeneralFail(caller, "BRIDGING_ELEMENT not enabled");
}
#endif

#ifdef __DEVELOPMENT__
#ifdef BRIDGING_ELEMENT
void FEExecutionManagerT::RunTHK(ifstreamT& in, ostream& status) const
{
	const char caller[] = "FEExecutionManagerT::RunTHK";

	/* output stream */
	StringT outfilename;
	outfilename.Root(in.filename());
	outfilename.Append(".out");
	ofstreamT out;
	out.open(outfilename);

#ifdef __MWERKS__
	if (!out.is_open())
	{
		cout << "\n " << caller << " : could not open file: " << outfilename << endl;
		return;
	}
#endif

	clock_t t0 = 0, t1 = 0, t2 = 0;

	/* start day/date info */
	time_t starttime;
	time(&starttime);

	int phase; // job phase
	try
	{
		t0 = clock();

		/* construction */
		phase = 0;
		in.set_marker('#');
		ifstreamT dummy_bridging_input; // TEMP - this would normally be input about ghost nodes
		FEManagerT_THK thk(in.filename(), out, fComm, fCommandLineOptions, dummy_bridging_input);
		thk.Initialize();

		t1 = clock();

		/* solution */
		phase = 1;
		thk.Solve();

		t2 = clock();
	}

	/* job failure */
	catch (ExceptionT::CodeT code)
	{
		status << "\n \"" << in.filename() << "\" exit on exception during the\n";
		if (phase == 0)
		{
			status << " construction phase. Check the input file for errors." << endl;
		
			/* echo some lines from input */
			if (code == ExceptionT::kBadInputValue) Rewind(in, status);
		}
		else
		{
			status << " solution phase. See \"" << outfilename << "\" for a list";
			status << " of the codes.\n";
		}
		
		/* fix clock values */
		if (t1 == 0) t1 = clock();
		if (t2 == 0) t2 = clock();		

		out << endl;
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);

	/* output timing */
	status << "\n     Filename: " << in.filename() << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   " Construction: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "     Solution: " << double(t2 - t1)/CLOCKS_PER_SEC << " sec.\n";
	
	out << "\n   Start time: " << ctime(&starttime);
	out <<   " Construction: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	out <<   "     Solution: "   << double(t2 - t1)/CLOCKS_PER_SEC << " sec.\n";
	status << "    Stop time: " << ctime(&stoptime);
	out   << "    Stop time: " << ctime(&stoptime);

	status << "\n End Execution\n" << endl;
	out    << "\n End Execution\n" << endl;
}
#else /* bridging element not enabled */
void FEExecutionManagerT::RunTHK(ifstreamT& in, ostream& status) const
{
#pragma unused(in)
#pragma unused(status)

	const char caller[] = "FEExecutionManagerT::RunTHK";
	ExceptionT::GeneralFail(caller, "BRIDGING_ELEMENT not enabled");
}
#endif /* BRIDGING_ELEMENT */
#endif /* __DEVELOPMENT__ */

/* dump current DTD or XSD file */
void FEExecutionManagerT::RunWriteDescription(int doc_type) const
{
	try {
//TEMP - parameters currently needed to construct an FEManagerT
	StringT input;
	ofstreamT output;
	CommunicatorT comm;
//TEMP

	/* collect parameters */
	cout << " collecting parameters..." << endl;
	ParameterTreeT tree;
	FEManagerT fe(input, output, comm, fCommandLineOptions);
	tree.BuildDescription(fe);
#ifdef BRIDGING_ELEMENT
	MultiManagerT multi(input, output, comm, fCommandLineOptions);
	tree.BuildDescription(multi);
#endif

	/* write description */
	cout << " writing description..." << endl;
	StringT out_path("tahoe");
	XML_Attribute_FormatterT::DocTypeT doc_type_ = XML_Attribute_FormatterT::Undefined;
	if (doc_type == XML_Attribute_FormatterT::DTD) {
		doc_type_ = XML_Attribute_FormatterT::DTD;
		out_path.Append(".dtd");
	}
	else if (doc_type == XML_Attribute_FormatterT::XSD) {
		doc_type_ = XML_Attribute_FormatterT::XSD;
		out_path.Append(".xsd");
	}
	XML_Attribute_FormatterT attribute(doc_type_);

	ofstreamT out;
	out.open(out_path);
	attribute.InitDescriptionFile(out);
	
	const ArrayT<ParameterListT*>& branches = tree.Branches();
	for (int i = 0; i < branches.Length(); i++)
		attribute.WriteDescription(out, *(branches[i]));

	/* write statistics */
	cout << " " << attribute.ElementCount() << " elements" << '\n';
	cout << " " <<  attribute.AttributeCount() << " attributes" << '\n';
	cout << " " <<  attribute.LimitCount() << " limits" << '\n';

	attribute.CloseDescriptionFile(out);
	out.close();
	cout << " wrote \"" << out_path << '"' << endl;	
	}
	
	catch (ExceptionT::CodeT exc) {
		cout << "\n FEExecutionManagerT::RunDTD: caught exception: " 
		     << ExceptionT::ToString(exc) << endl;
	}
}

/* recursive dispatch */
void FEExecutionManagerT::JobOrBatch(ifstreamT& in, ostream& status)
{
	StringT ext;
	ext.Suffix(in.filename());
	if (ext == ".xml")
		RunJob(in, status);
	else /* inherited */
		ExecutionManagerT::JobOrBatch(in, status);
}

void FEExecutionManagerT::RunJob_serial(const StringT& input_file, ostream& status) const
{
	const char* caller = "FEExecutionManagerT::RunJob_serial_XML";

	/* output stream */
	StringT outfilename;
	outfilename.Root(input_file);
	outfilename.Append(".out");
	ofstreamT out;
	out.open(outfilename);

#ifdef __MWERKS__
	if (!out.is_open()) {
		cout << "\n " << caller << ": could not open file: " << outfilename << endl;
		return;
	}
#endif

	clock_t t0 = 0, t1 = 0, t2 = 0;

	/* start day/date info */
	time_t starttime;
	time(&starttime);

	int phase; // job phase
	try
	{
		t0 = clock();
		phase = 0;

		/* generate validated parameter list */
		ParameterListT valid_list;
		FEManagerT::ParseInput(input_file, valid_list, true, true, true, fCommandLineOptions);

		/* write the validated list as formatted text */
		if (true) {
			DotLine_FormatterT pp_format;
			pp_format.SetTabWidth(4);
			pp_format.InitParameterFile(out);
			pp_format.WriteParameterList(out, valid_list);
			pp_format.CloseParameterFile(out);
			out << endl;
		}
		
		/* analysis type */
		if (valid_list.Name() == "tahoe")
		{
			/* construction */
			FEManagerT analysis1(input_file, out, fComm, fCommandLineOptions);
			analysis1.TakeParameterList(valid_list);
			t1 = clock();

#if defined(__MWERKS__) && __option(profile)
			/* start recording profiler information */
			ProfilerSetStatus(1);
#endif
		
			/* solution */
			phase = 1;
			analysis1.Solve();

#if defined(__MWERKS__) && __option(profile)
			/* stop recording profiler information */
			ProfilerSetStatus(0);
#endif
			t2 = clock();
		}
#ifdef BRIDGING_ELEMENT
		else if (valid_list.Name() == "tahoe_multi")
		{
			/* construction */
			MultiManagerT analysis1(input_file, out, fComm, fCommandLineOptions);
			analysis1.TakeParameterList(valid_list);
			t1 = clock();

#if defined(__MWERKS__) && __option(profile)
			/* start recording profiler information */
			ProfilerSetStatus(1);
#endif
		
			/* solution */
			phase = 1;
			analysis1.Solve();

#if defined(__MWERKS__) && __option(profile)
			/* stop recording profiler information */
			ProfilerSetStatus(0);
#endif
			t2 = clock();
		}
#endif
		else
			ExceptionT::GeneralFail(caller, "unexpected parameter list \"%s\"", valid_list.Name().Pointer());
	}

	/* job failure */
	catch (ExceptionT::CodeT code)
	{
		status << "\n \"" << input_file << "\" exit on exception during the\n";
		if (phase == 0)
		{
			status << " construction phase. Check the input file for errors." << endl;
		}
		else
		{
			status << " solution phase. See \"" << outfilename << "\" for a list";
			status << " of the codes.\n";
		}
		
		/* fix clock values */
		if (t1 == 0) t1 = clock();
		if (t2 == 0) t2 = clock();		

		out << endl;
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);

	/* output timing */
	status << "\n     Filename: " << input_file << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   " Construction: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "     Solution: " << double(t2 - t1)/CLOCKS_PER_SEC << " sec.\n";
	
	out << "\n   Start time: " << ctime(&starttime);
	out <<   " Construction: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	out <<   "     Solution: "   << double(t2 - t1)/CLOCKS_PER_SEC << " sec.\n";
	status << "    Stop time: " << ctime(&stoptime);
	out   << "    Stop time: " << ctime(&stoptime);

	status << "\n End Execution\n" << endl;
	out    << "\n End Execution\n" << endl;
}	

/* generate decomposition files */
void FEExecutionManagerT::RunDecomp_serial(const StringT& input_file, ostream& status,CommunicatorT& comm, int size) const
{
	/* look for size */
	int index;
	if (CommandLineOption("-decomp", index) && fCommandLineOptions.Length() > index+1) {
		const char* opt = fCommandLineOptions[index+1];
		if (strlen(opt) > 1 && isdigit(opt[1]))
			size = atoi(opt+1); /* opt[0] = '-' */
	}

	/* look for method */
	int method = -1;
	if (CommandLineOption("-decomp_method", index))
	{	
		if (fCommandLineOptions.Length() > index+1)
		{
			const char* opt = fCommandLineOptions[index+1];
			if (strlen(opt) > 1 && isdigit(opt[1]))
				method = atoi(opt+1); /* opt[0] = '-' */
		}
	}
	
	/* prompt for size */
	if (size == -1) 
	{
		/* prompt for decomp size */
		int count = 0;
		while (count == 0 || (count < 5 && size < 2))
		{
			count++;

			/* number of partitions */					
			cout << "\n Enter number of partitions > 1 (0 to quit): ";
#if (defined __SGI__ && defined __TAHOE_MPI__)
			cout << '\n';
#endif
			cin >> size;
		
			/* clear to end of line */
			fstreamT::ClearLine(cin);

			if (size == 0) break;
		}
	}
	if (size < 2) return;

	/* prompt for method */
	if (method == -1)
	{
		cout << "\n Select partitioning method:\n"
		     << '\t' << PartitionT::kGraph   << ": graph\n"
		     << '\t' << PartitionT::kAtom    << ": atom\n"
		     << '\t' << PartitionT::kSpatial << ": spatial\n";
		cout << "\n method: "; 
#if (defined __SGI__ && defined __TAHOE_MPI__)
		cout << '\n';
#endif					
		cin >> method;
	
		/* clear to end of line */
		fstreamT::ClearLine(cin);
	}

	/* time markers */
	clock_t t0 = 0, t1 = 0;

	/* start day/date info */
	time_t starttime;
	time(&starttime);
	try
	{
		t0 = clock();

		/* path to parameters file */
		StringT path;
		path.FilePath(input_file);
		
		/* generate validated parameter list */
		ParameterListT valid_list;
		FEManagerT::ParseInput(input_file, valid_list, true, false, false, fCommandLineOptions);
		
		/* extract model file and file format */
		int i_format = valid_list.GetParameter("geometry_format");
		IOBaseT::FileTypeT format = IOBaseT::int_to_FileTypeT(i_format);
		StringT model_file = valid_list.GetParameter("geometry_file");

		/* name translation */
		model_file.ToNativePathName();      
		model_file.Prepend(path);

		/* set output map and and generate decomposition */
		Decompose(input_file, size, method, comm, model_file, format);
		t1 = clock();
	}

	/* job failure */
	catch (ExceptionT::CodeT code)
	{
		/* fix clock values */
		if (t1 == 0) t1 = clock();
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);

	/* output timing */
	status << "\n     Filename: " << input_file << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   "Decomposition: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "    Stop time: " << ctime(&stoptime);	
}

/* join parallel results files */
void FEExecutionManagerT::RunJoin_serial(const StringT& input_file, ostream& status, int size) const
{
	/* time markers */
	clock_t t0 = 0, t1 = 0;

	/* start day/date info */
	time_t starttime;
	time(&starttime);
	try
	{
		t0 = clock();

		/* path to parameters file */
		StringT path;
		path.FilePath(input_file);

		/* generate validated parameter list */
		ParameterListT valid_list;
		FEManagerT::ParseInput(input_file, valid_list, true, false, false, fCommandLineOptions);
		
		/* model file parameters */
		int i_format = valid_list.GetParameter("geometry_format");
		IOBaseT::FileTypeT model_format = IOBaseT::int_to_FileTypeT(i_format);
		i_format = valid_list.GetParameter("output_format");
		IOBaseT::FileTypeT results_format = IOBaseT::int_to_FileTypeT(i_format);
		StringT model_file = valid_list.GetParameter("geometry_file");
		if (results_format == IOBaseT::kTahoe ||
		    results_format == IOBaseT::kTahoeII)
			results_format = IOBaseT::kTahoeResults;

		/* name translation */
		model_file.ToNativePathName();      
		model_file.Prepend(path);

		int index;
		if (CommandLineOption("-join", index) && fCommandLineOptions.Length() > index+1) {
			const char* opt = fCommandLineOptions[index+1];
			if (strlen(opt) > 1 && isdigit(opt[1]))
				size = atoi(opt+1);
		}

		/* prompt for decomp size */
		if (size == -1)
		{
			cout << "\n Enter number of partitions > 1 (0 to quit): ";
#if (defined __SGI__ && defined __TAHOE_MPI__)
			cout << '\n';
#endif					
			cin >> size;

			/* clear to end of line */
			fstreamT::ClearLine(cin);
		}

		if (size < 2) return;
		
		/* construct output formatter */
		StringT program = "tahoe";
		StringT version = FEManagerT::Version();
		StringT title;
		const ParameterT* title_param = valid_list.Parameter("title");
		if (title_param)
			title = *title_param;
		StringT input = input_file;
		OutputBaseT* output = IOBaseT::NewOutput(program, version, title, input, 
			results_format, cout);

		/* construct joiner */
		JoinOutputT output_joiner(input_file, model_file, model_format, 
			results_format, output, size);

		/* join files */
		output_joiner.Join();
		
		/* free output formatter */
		delete output;

		t1 = clock();
	}

	/* job failure */
	catch (ExceptionT::CodeT code)
	{
		/* fix clock values */
		if (t1 == 0) t1 = clock();
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);

	/* output timing */
	status << "\n     Filename: " << input_file << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   "         Join: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "    Stop time: " << ctime(&stoptime);	
}

/* testing for distributed execution */
void FEExecutionManagerT::RunJob_parallel(const StringT& input_file, ostream& status) const
{
	const char caller[] = "::RunJob_parallel";

	/* get rank and size */
	int rank = Rank();
	int size = Size();

	/* time markers */
	clock_t t0 = 0, t1 = 0, t2 = 0, t3 = 0;

	/* start day/date info */
	time_t starttime;
	time(&starttime);

	/* output stream */
	StringT out_file;
	out_file.Root(input_file);
	out_file.Append(".p", rank);
	out_file.Append(".out");
	ofstreamT out;
	out.open(out_file);

	int token; // for run time check sums
	try {
	t0 = clock();
	
	/* generate validated parameter list */
	ParameterListT valid_list;
	FEManagerT::ParseInput(input_file, valid_list, true, false, false, fCommandLineOptions);
	if (true) /* write the validated list as formatted text */ {
		DotLine_FormatterT pp_format;
		pp_format.SetTabWidth(4);
		pp_format.InitParameterFile(out);
		pp_format.WriteParameterList(out, valid_list);
		pp_format.CloseParameterFile(out);
		out << endl;
	}
		
	/* extract model file and file format */
	int i_format = valid_list.GetParameter("geometry_format");
	IOBaseT::FileTypeT format = IOBaseT::int_to_FileTypeT(i_format);
	StringT model_file = valid_list.GetParameter("geometry_file");

	/* name translation */
	StringT path;
	path.FilePath(input_file);
	model_file.ToNativePathName();      
	model_file.Prepend(path);

	/* generate decomposition if needed */
	token = 1;
	if (NeedDecomposition(model_file, size) || !CommandLineOption("-split_io"))
	{
		/* 'serial' communicator */
		int rank = fComm.Rank();
		CommunicatorT comm(fComm, (rank == 0) ? rank : CommunicatorT::kNoColor);

		/* decompose on rank 0 */
		if (rank == 0) RunDecomp_serial(input_file, status, comm, fComm.Size());
	
		/* synch */
		fComm.Barrier();	
	}

	/* synch and check status */
	if (fComm.Sum(token) != size) ExceptionT::GeneralFail();

	/* read partition information */
	PartitionT partition;
	StringT part_file;
	part_file.Root(model_file);
	part_file.Append(".n", size);
	part_file.Append(".part", rank);
	ifstreamT part_in('#', part_file);
	token = 1;
	if (!part_in.is_open())
	{
		cout << "\n " << caller << ": could not open file: " << part_file << endl;
		token = 0;
	}
	else
	{
		part_in >> partition;
		if (partition.ID() != rank)
		{
			cout << "\n " << caller << ": rank and partition ID mismatch" << endl;
			token = 0;
		}
		
		/* set correct numbering scope */
		partition.SetScope(PartitionT::kLocal);
	}
	
	/* synch and check status */
	if (fComm.Sum(token) != size) ExceptionT::GeneralFail();
		
	/* write partial geometry files (if needed) */
	StringT suffix;
	suffix.Suffix(model_file);
	StringT partial_file;
	partial_file.Root(model_file);
	partial_file.Append(".n", size);
	partial_file.Append(".p", rank);
	partial_file.Append(suffix);
	if (NeedModelFile(partial_file, format))
	{
		/* original model file */
		ModelManagerT model_ALL(cout);
		if (!model_ALL.Initialize(format, model_file, true))
			ExceptionT::GeneralFail(caller, 
				"error opening file: %s", (const char*) model_file);

		cout << "\n " << caller << ": writing partial geometry file: " << partial_file << endl;
		EchoPartialGeometry(partition, model_ALL, partial_file, format);
		cout << " " << caller << ": writing partial geometry file: partial_file: "
			 << partial_file << ": DONE" << endl;
	}
	
	/* construct local problem (Initialize() changes the file name) */
	t1 = clock();
	token = 1;

	/* construct solver */
	FEManagerT_mpi FEman(input_file, out, fComm, fCommandLineOptions, &partition, FEManagerT_mpi::kRun);
	try { 
		FEman.TakeParameterList(valid_list); 
	}
	catch (ExceptionT::CodeT code) {
		status << "\n \"" << input_file << "\" exit on exception " << code << " during the\n";
		status << " construction phase. Check the input file for errors." << endl;
		token = 0;
	}
	
	if (fComm.Sum(token) != size)
	{
		/* gather tokens to rank 0 */
		if (rank == 0)
		{
			iArrayT tokens(size);
			fComm.Gather(token, tokens);
		
			status << "\n The following processes exit on exception during construction:\n";
			for (int i = 0; i < tokens.Length(); i++)
				if (tokens[i] != 1)
					status << setw(kIntWidth) << i << '\n';
			status.flush();
		}
		else fComm.Gather(token, 0);
		
		ExceptionT::GeneralFail(caller);
	}
	
	/* external IO */
	token = 1;
	IOManager_mpi* IOMan = NULL;
	if (partition.DecompType() == PartitionT::kGraph && !CommandLineOption("-split_io"))
	{
		try {
			/* set-up local IO */
			IOMan = new IOManager_mpi(input_file, fComm, *(FEman.OutputManager()), FEman.Partition(), model_file, format);
			if (!IOMan) throw ExceptionT::kOutOfMemory;
		
			/* set external IO */
			FEman.SetExternalIO(IOMan);
		}
	
		catch (ExceptionT::CodeT code)
		{
			token = 0;
			status << "\n \"" << input_file << "\" exit on exception " << code 
			       << " setting the external IO" << endl;
		}
	}

	if (fComm.Sum(token) != size)
	{
		/* gather tokens to rank 0 */
		if (rank == 0)
		{
			iArrayT tokens(size);
			fComm.Gather(token, tokens);			
			status << "\n The following processes exit on exception during construction:\n";
			for (int i = 0; i < tokens.Length(); i++)
				if (tokens[i] != 1)
					status << setw(kIntWidth) << i << '\n';
			status.flush();
		}
		else 
			fComm.Gather(token, 0);
		
		ExceptionT::GeneralFail(caller);
	}

	/* solve */
	t2 = clock();
	token = 1;
	try { FEman.Solve(); }
	catch (ExceptionT::CodeT code)
	{
		status << "\n \"" << input_file << "\" exit on exception " << code << " during the\n";
		status << " solution phase. See \"" << out_file << "\" for a list";
		status << " of the codes.\n";
		token = 0;
	}

	/* free external IO */
	delete IOMan;

	if (fComm.Sum(token) != size)
	{
		/* gather tokens to rank 0 */
		if (rank == 0)
		{
			iArrayT tokens(size);
			fComm.Gather(token, tokens);

			status << "\n The following processes exit on exception during solution:\n";
			for (int i = 0; i < tokens.Length(); i++)
				if (tokens[i] != 1)
					status << setw(kIntWidth) << i << '\n';
			status.flush();
		}
		else fComm.Gather(token, 0);
		
		ExceptionT::GeneralFail(caller);
	}
	t3 = clock(); } // end try

	/* job failure */
	catch (ExceptionT::CodeT code)
	{
		/* fix clock values */
		if (t1 == 0) t1 = clock();
		if (t2 == 0) t2 = clock();		
		if (t3 == 0) t3 = clock();		
	}

	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);

	/* output timing */
	status << "\n     Filename: " << input_file << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   "Decomposition: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   " Construction: " << double(t2 - t1)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "     Solution: " << double(t3 - t2)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "    Stop time: " << ctime(&stoptime);
	
	out << "\n   Start time: " << ctime(&starttime);
	out <<   "Decomposition: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	out <<   " Construction: " << double(t2 - t1)/CLOCKS_PER_SEC << " sec.\n";
	out <<   "     Solution: " << double(t3 - t2)/CLOCKS_PER_SEC << " sec.\n";
	out <<   "    Stop time: " << ctime(&stoptime);

	status << "\n End Execution\n" << endl;
	out    << "\n End Execution\n" << endl;
}	

/* print message on exception */
void FEExecutionManagerT::Rewind(ifstreamT& in, ostream& status) const
{
	/* reset stream */
	ifstream& istr = in;
	if (!istr.good()) istr.clear();
			
	/* rewind a couple of lines */
	in.rewind(4);
			
	status << " Error occurred while reading input near:\n\n";
	int i = 8;
	char line[255];
	istr.getline(line, 254);
	while(istr.good() && i-- > 0)
	{
		status << line << '\n';
		istr.getline(line, 254);
	}
	status << endl;
}

void FEExecutionManagerT::Decompose(const StringT& input_file, int size, int decomp_type, CommunicatorT& comm,
	const StringT& model_file, IOBaseT::FileTypeT format) const
{	
	/* dispatch */
	switch (decomp_type)
	{
		case PartitionT::kGraph:
			Decompose_graph(input_file, size, comm, model_file, format);
			break;

		case PartitionT::kAtom:
			Decompose_atom(input_file, size, model_file, format);
			break;
			
		case PartitionT::kSpatial:
			cout << "\n FEExecutionManagerT::Decompose: spatial decomposition not implemented yet" << endl;
			break;
						
		default:
			cout << "\n FEExecutionManagerT::Decompose: unrecognized method: " << decomp_type << endl;
	}
}

void FEExecutionManagerT::Decompose_atom(const StringT& input_file, int size,
	const StringT& model_file, IOBaseT::FileTypeT format) const
{
#pragma unused(input_file)
	const char caller[] = "FEExecutionManagerT::Decompose_atom";

	/* files exist */
	bool need_decomp = NeedDecomposition(model_file, size);
	if (!need_decomp)
	{
		cout << "\n " << caller <<": decomposition files exist" << endl;
		return;
	}
	
	/* model manager */
	ModelManagerT model(cout);
	if (!model.Initialize(format, model_file, true))
		ExceptionT::BadInputValue(caller, "could not open model file: %s", (const char*) model_file);

	/* dimensions */
	int nnd = model.NumNodes();
	int nsd = model.NumDimensions();

	/* node-to-partition map */
	iArrayT part_map(nnd);
	
	/* number of nodes per processor - p_i = i/proc_size */
	int part_size = nnd/size;
	if (part_size < 1) ExceptionT::GeneralFail();
	
	/* labels nodes */
	int part = 0;
	int count = 0;
	for (int i = 0; i < nnd; i++)
	{
		part_map[i] = part;
		if (++count == part_size) {
			if (part < size - 1) {
				part++;
				count = 0;
			}
		}
	}

	/* get pointers to all blocks */
	const ArrayT<StringT>& IDs = model.ElementGroupIDs();
	ArrayT<const iArray2DT*> connects_1(IDs.Length());
	model.ElementGroupPointers(IDs, connects_1);
	ArrayT<const RaggedArray2DT<int>*> connects_2;

	/* set partition information and write partial geometry files*/
	for (int i = 0; i < size; i++)
	{
		/* partition data */
		PartitionT partition;

		/* mark nodes */
		partition.Set(size, i, part_map, connects_1, connects_2);
		
		/* set elements */
		const ArrayT<StringT>& elem_ID = model.ElementGroupIDs();
		partition.InitElementBlocks(elem_ID);
		for (int j = 0; j < elem_ID.Length(); j++)
		{
			const iArray2DT& elems = model.ElementGroup(elem_ID[j]);
			partition.SetElements(elem_ID[j], elems);
		}

		/* set to local scope */
		partition.SetScope(PartitionT::kLocal);
		
		/* set decomposition type */
		partition.SetDecompType(PartitionT::kAtom);
	
		/* output file name */
		StringT geom_file, suffix;
		suffix.Suffix(model_file);
		geom_file.Root(model_file);
		geom_file.Append(".n", size);
		geom_file.Append(".p", i);
		geom_file.Append(suffix);
				
		cout << "     Writing partial model file: " << geom_file << endl;
		try { EchoPartialGeometry(partition, model, geom_file, format); }
		catch (ExceptionT::CodeT error)
		{
			ExceptionT::Throw(error, caller, "exception writing file: %s", (const char*) geom_file);
		}
		
		/* partition information */
		StringT part_file;
		part_file.Root(model_file);
		part_file.Append(".n", size);
		part_file.Append(".part", i);

		ofstream part_out(part_file);
		part_out << "# data for partition: " << i << '\n';
		part_out << partition << '\n';
		part_out.close();
	}
}

void FEExecutionManagerT::Decompose_spatial(const StringT& input_file, int size,
	const StringT& model_file, IOBaseT::FileTypeT format) const
{
#pragma unused(input_file)
#pragma unused(size)
#pragma unused(model_file)
#pragma unused(format)
	cout << "\n FEExecutionManagerT::Decompose_spatial: not implemented" << endl;
}

/* graph-based decomposition */
void FEExecutionManagerT::Decompose_graph(const StringT& input_file, int size,
	CommunicatorT& comm, const StringT& model_file, IOBaseT::FileTypeT format) const
{
	const char caller[] = "FEExecutionManagerT::Decompose_graph";

	bool need_decomp = NeedDecomposition(model_file, size);
	if (need_decomp)
	{
		/* echo stream */
		StringT decomp_file;
		decomp_file.Root(input_file);
		decomp_file.Append(".out");
		ofstreamT decomp_out;
		decomp_out.open(decomp_file);
		
		/* generate validated parameter list */
		ParameterListT valid_list;
		FEManagerT::ParseInput(input_file, valid_list, true, false, false, fCommandLineOptions);

		/* construct global problem */
		FEManagerT_mpi global_FEman(input_file, decomp_out, comm, fCommandLineOptions, NULL, FEManagerT_mpi::kDecompose);
		try { 
			global_FEman.TakeParameterList(valid_list); 
		}
		catch (ExceptionT::CodeT code) {
			ExceptionT::Throw(code, caller, "exception \"%s\" constructing global problem",
				ExceptionT::ToString(code));
		}

		/* decompose */
		ArrayT<PartitionT> partition(size);
		if (need_decomp)
		{
			int method = 0;

/* use METIS by default */
#ifdef __METIS__
			if (!CommandLineOption("-no_metis"))
			{
				cout << "\n using METIS. Disable with command-line option \"-no_metis\"" << endl;
				method = 1;
			}
#endif
			/* graph object */
			GraphT graph;	
			try {
				cout << "\n Decomposing: " << model_file << endl;
				global_FEman.Decompose(partition, graph, true, method);
				cout << " Decomposing: " << model_file << ": DONE"<< endl;
			}
			catch (ExceptionT::CodeT code) {
				ExceptionT::Throw(code, caller, "exception \"%s\" during decomposition",
					ExceptionT::ToString(code));
			}
			
			/* write partition data out */
			for (int q = 0; q < partition.Length(); q++)
			{
				/* set to local scope */
				partition[q].SetScope(PartitionT::kLocal);

				/* set decomposition type */
				partition[q].SetDecompType(PartitionT::kGraph);

				StringT file_name;
				file_name.Root(model_file);
				file_name.Append(".n", partition.Length());
				file_name.Append(".part", q);
		
				ofstream out_q(file_name);
				out_q << "# data for partition: " << q << '\n';
				out_q << partition[q] << '\n';
				out_q.close();
			}
			
			/* write decomposition map */
			if (format == IOBaseT::kExodusII)
			{
				/* file name */
				StringT map_file;
				map_file.Root(model_file);
				map_file.Append(".n", size);
				map_file.Append(".decomp.exo");
				cout << " Node map file: " << map_file << endl;
			
				/* database */
				ExodusT exo(cout);
				StringT str;
				ArrayT<StringT> str_list;
				const dArray2DT& coords = global_FEman.ModelManager()->Coordinates();
				int nnd = coords.MajorDim();
				int nsd = coords.MinorDim();
				exo.Create(map_file, str, str_list, str_list, nsd, nnd,
					nnd, size, 0, 0);
			
				/* coordinates */
				exo.WriteCoordinates(coords);
				
				/* data in decomp file */
				int n_internal = 0;
				int n_border = 0;
				int n_external = 0;
				dArrayT part(nnd); part = -1;
				dArrayT inex(nnd); inex = 0.0;
				for (int i = 0; i < size; i++)
				{
					/* "owned" nodes */
					const iArrayT& nd_i = partition[i].Nodes_Internal();
					const iArrayT& nd_b = partition[i].Nodes_Border();
					iArray2DT connects(nd_i.Length() + nd_b.Length(), 1);
					connects.CopyPart(0, nd_i, 0, nd_i.Length());
					connects.CopyPart(nd_i.Length(), nd_b, 0, nd_b.Length());
					
					/* increment counts */
					n_internal += nd_i.Length();
					n_border += nd_b.Length();
					n_external += partition[i].Nodes_External().Length();
					
					/* convert to global numbering */
					if (partition[i].NumberScope() != PartitionT::kGlobal)
						partition[i].SetNodeScope(PartitionT::kGlobal, connects);
				
					/* label part */
					for (int j = 0; j < connects.Length(); j++)
						part[connects[j]] = i;
				
					/* label border nodes */
					int* pnd_b = connects.Pointer(nd_i.Length());
					for (int k = 0; k < nd_b.Length(); k++)
						inex[*pnd_b++] = 1;
				
					/* write to file */
					connects++;
					exo.WriteConnectivities(i+1, GeometryT::kPoint, connects);
				}
				
				/* degree of each node */
				iArrayT i_degree(nnd);
				int shift;
				graph.Degrees(i_degree, shift);
				if (shift != 0)				
					ExceptionT::GeneralFail(caller, "unexpected node number shift: %d", shift);

				dArrayT degree(nnd);
				for (int j = 0; j < nnd; j++)
					degree[j] = i_degree[j];

				/* part labels */
				ArrayT<StringT> labels(3);
				labels[0] = "part";
				labels[1] = "inex";
				labels[2] = "degree";
				exo.WriteLabels(labels, ExodusT::kNode);
				exo.WriteTime(1, 0.0);
				exo.WriteNodalVariable(1, 1, part);
				exo.WriteNodalVariable(1, 2, inex);
				exo.WriteNodalVariable(1, 3, degree);
				
				/* write statistics */
				cout << " Statistics:\n" 
				     << "     total number of nodes = " << coords.MajorDim() << '\n'
				     << "            internal nodes = " << n_internal << '\n'
				     << "              border nodes = " << n_border << '\n'
				     << "            external nodes = " << n_external << '\n';

				/* done */
				cout << " Node map file: " << map_file << ": DONE"<< endl;
			}
		}
		else
		{
			/* read partition information from stream */
			for (int i = 0; i < partition.Length(); i++)
			{
				StringT part_file;
				part_file.Root(model_file);
				part_file.Append(".n", size);
				part_file.Append(".part", i);
				ifstreamT part_in('#', part_file);
				
				part_in >> partition[i];
				partition[i].SetScope(PartitionT::kLocal);
			}
		}
		
		/* write partial geometry files */
		if (CommandLineOption("-decomp"))
		{
			/* model manager for the total geometry - can't use the one from global_FEman 
			 * because it may contain runtime-generated connectivities not in the original 
			 * model file */
			ModelManagerT model_ALL(cout);
			if (!model_ALL.Initialize(format, model_file, true))
				ExceptionT::GeneralFail(caller, "error opening file: %s", (const char*) model_file);
	
			/* write partial geometry files */
			for (int i = 0; i < partition.Length(); i++)
			{
				StringT partial_file, suffix;
				suffix.Suffix(model_file);
				partial_file.Root(model_file);
				partial_file.Append(".n", size);
				partial_file.Append(".p", i);
				partial_file.Append(suffix);
				
				if (NeedModelFile(partial_file, format))
				{			
					cout << "     Writing partial model file: " << partial_file << endl;
					try { 
						EchoPartialGeometry(partition[i], model_ALL, partial_file, format); 
					}
					catch (ExceptionT::CodeT error) {
						ExceptionT::Throw(error, caller, "exception writing file \"%s\"",
							partial_file.Pointer());
					}
				}
			}
		}
	}
	else
		cout << "\n " << caller << ": decomposition files exist" << endl;
}

/* returns 1 if a new decomposition is needed */
bool FEExecutionManagerT::NeedDecomposition(const StringT& model_file,
	int size) const
{
	/* model file root */
	StringT root;
	root.Root(model_file);

	for (int i = 0; i < size; i++)
	{
		/* partition data file name */
		StringT part_file = root;
		part_file.Append(".n", size);
		part_file.Append(".part", i);
	
		/* open partition file */
		ifstreamT part_in(PartitionT::CommentMarker(), part_file);
		if (part_in.is_open())
		{
			StringT version;
			part_in >> version;
		
			int num_parts, part_ID;
			part_in >> num_parts;
			part_in >> part_ID;

			if (!PartitionT::CheckVersion(version) ||
			     num_parts != size ||
			     part_ID != i) return true;
		}
		else
			return true;
	}
	return false;
}

/* returns true if the global output model file is not found */
bool FEExecutionManagerT::NeedModelFile(const StringT& model_file,
	IOBaseT::FileTypeT format) const
{
	switch (format)
	{
		case IOBaseT::kTahoeII:
		{
			ModelFileT file;
			if (file.OpenRead(model_file) == ModelFileT::kOK)
				return false;
			else
				return true;
		}	
		case IOBaseT::kExodusII:
		{
			ExodusT file(cout);
			if (file.OpenRead(model_file))
				return false;			
			else
				return true;
		}
		default:
		{
			ifstreamT test(model_file);
			if (test.is_open())
				return false;
			else
				return true;
		}
	}
}

/* write partial geometry files */
void FEExecutionManagerT::EchoPartialGeometry(const PartitionT& partition,
	ModelManagerT& model_ALL, const StringT& partial_file,
	IOBaseT::FileTypeT format) const
{
	switch (format)
	{
		case IOBaseT::kExodusII:
			EchoPartialGeometry_ExodusII(partition, model_ALL, partial_file);
			break;
	
		case IOBaseT::kTahoeII:
			EchoPartialGeometry_TahoeII(partition, model_ALL, partial_file);
			break;

		default:	
			ExceptionT::GeneralFail("FEExecutionManagerT::EchoPartialGeometry", 
				"unsupported file format: %d", format);
	}
}

void FEExecutionManagerT::EchoPartialGeometry_ExodusII(const PartitionT& partition,
	ModelManagerT& model_ALL, const StringT& partial_file) const
{
	/* partition */
	int part = partition.ID();

	/* collect file creation information */
	StringT title = "partition file: ";
	title.Append(part);
	ArrayT<StringT> nothing;
	const iArrayT& node_map = partition.NodeMap();
	int num_node = node_map.Length();
	int num_dim  = model_ALL.NumDimensions();
	int num_blks = model_ALL.NumElementGroups();
	int num_ns   = model_ALL.NumNodeSets();
	int num_ss   = model_ALL.NumSideSets();
	int num_elem = model_ALL.NumElements();
	
	/* partial model file */
	ExodusT model(cout);
	model.Create(partial_file, title, nothing, nothing, num_dim, num_node,
		num_elem, num_blks, num_ns, num_ss);
	
	/* coordinates */
	const dArray2DT& coords_ALL = model_ALL.Coordinates();
	dArray2DT coords(node_map.Length(), coords_ALL.MinorDim());		
	coords.RowCollect(node_map, coords_ALL);
	model.WriteCoordinates(coords);
	coords.Free();
		
	/* element sets */
	const ArrayT<StringT>& elem_ID = model_ALL.ElementGroupIDs();
	for (int j = 0; j < elem_ID.Length(); j++)
	{
		/* read global block */
		const iArray2DT& set_ALL = model_ALL.ElementGroup(elem_ID[j]);

		/* collect connectivities within the partition */
		const iArrayT& element_map = partition.ElementMap(elem_ID[j]);
		iArray2DT set(element_map.Length(), set_ALL.MinorDim());
		set.RowCollect(element_map, set_ALL);

		/* map to local scope */
		partition.SetNodeScope(PartitionT::kLocal, set);

		/* write to file */
		set++;
		GeometryT::CodeT geometry_code = model_ALL.ElementGroupGeometry(elem_ID[j]);
		int id = atoi(elem_ID[j]);
		model.WriteConnectivities(id, geometry_code, set);
	}
		
	/* node sets */
	const ArrayT<StringT>& node_ID = model_ALL.NodeSetIDs();
	for (int k = 0; k < node_ID.Length(); k++)
	{
		/* whole node set */
		const iArrayT& nodeset_ALL = model_ALL.NodeSet(node_ID[k]);
				
		/* partition node set */
		iArrayT local_indices;
		partition.ReturnPartitionNodes(nodeset_ALL, local_indices);
		
		/* non-empty set */
		if (1 || local_indices.Length() > 0)
		{
			iArrayT nodeset(local_indices.Length());
			nodeset.Collect(local_indices, nodeset_ALL);
			
			/* map to local numbering */
			partition.SetNodeScope(PartitionT::kLocal, nodeset);
				
			/* add */
			int id = atoi(node_ID[k]);
			nodeset++;
			model.WriteNodeSet(id, nodeset);
		}
	}

	/* side sets */		
	const ArrayT<StringT>& side_ID = model_ALL.SideSetIDs();
	for (int l = 0; l < side_ID.Length(); l++)
	{
		/* whole side set */
		iArray2DT sideset_ALL = model_ALL.SideSet(side_ID[l]);
		const StringT& element_set_ID = model_ALL.SideSetGroupID(side_ID[l]);

		/* partition side set */
		iArray2DT sideset;
		if (sideset_ALL.MajorDim() > 0)
		{
			iArrayT elements_ALL(sideset_ALL.MajorDim());
			sideset_ALL.ColumnCopy(0, elements_ALL);
//			elements_ALL--;
			sideset_ALL.SetColumn(0, elements_ALL);
				
			iArrayT local_indices;
			partition.ReturnPartitionElements(element_set_ID, elements_ALL, local_indices);
						
			/* non-empty set */
			if (local_indices.Length() > 0)
			{
				sideset.Dimension(local_indices.Length(), sideset_ALL.MinorDim());
				sideset.RowCollect(local_indices, sideset_ALL);

				iArrayT elements(sideset.MajorDim());
				sideset.ColumnCopy(0, elements);
				partition.SetElementScope(PartitionT::kLocal, element_set_ID, elements);
//				elements++;
				sideset.SetColumn(0, elements);
			}

			/* add */
			sideset++;
			int ss_id = atoi(side_ID[l]);
			int el_id = atoi(element_set_ID);
			model.WriteSideSet(ss_id, el_id, sideset);
		}			
	}
}

void FEExecutionManagerT::EchoPartialGeometry_TahoeII(const PartitionT& partition,
	ModelManagerT& model_ALL, const StringT& partial_file) const
{
	/* partition */
	int part = partition.ID();
	
	/* open model file */
	bool extern_file = true;
	ModelFileT model;
	model.OpenWrite(partial_file, extern_file);

	/* title */
	StringT title;
	title.Append("partition ", part);
	if (model.PutTitle(title) != ModelFileT::kOK) ExceptionT::GeneralFail();

	/* nodal coordinates */
	const iArrayT& node_map = partition.NodeMap();
	const dArray2DT& coords_ALL = model_ALL.Coordinates();
	dArray2DT coords(node_map.Length(), coords_ALL.MinorDim());		
	coords.RowCollect(node_map, coords_ALL);
	if (model.PutCoordinates(coords) != ModelFileT::kOK) ExceptionT::GeneralFail();
	coords.Free();	
		
	/* element sets */
	const ArrayT<StringT>& elem_ID = model_ALL.ElementGroupIDs();
	for (int j = 0; j < elem_ID.Length(); j++)
	{
		/* read global block */
		const iArray2DT& set_ALL = model_ALL.ElementGroup(elem_ID[j]);

		/* collect connectivities within the partition */
		const iArrayT& element_map = partition.ElementMap(elem_ID[j]);
		iArray2DT set(element_map.Length(), set_ALL.MinorDim());
		set.RowCollect(element_map, set_ALL);
			
		/* map to local scope */
		partition.SetNodeScope(PartitionT::kLocal, set);

		/* add */
		int id = atoi(elem_ID[j]);
		set++;
		if (model.PutElementSet(id, set) != ModelFileT::kOK) ExceptionT::GeneralFail();
	}
		
	/* node sets */
	const ArrayT<StringT>& node_ID = model_ALL.NodeSetIDs();
	for (int k = 0; k < node_ID.Length(); k++)
	{
		/* whole node set */
		const iArrayT& nodeset_ALL = model_ALL.NodeSet(node_ID[k]);
				
		/* partition node set */
		iArrayT local_indices;
		partition.ReturnPartitionNodes(nodeset_ALL, local_indices);
		iArrayT nodeset(local_indices.Length());
		nodeset.Collect(local_indices, nodeset_ALL);
			
		/* map to local numbering */
		partition.SetNodeScope(PartitionT::kLocal, nodeset);
			
		/* add */
		int id = atoi(node_ID[k]);
		nodeset++;
		if (model.PutNodeSet(id, nodeset) != ModelFileT::kOK) ExceptionT::GeneralFail();
	}

	/* side sets */		
	const ArrayT<StringT>& side_ID = model_ALL.SideSetIDs();
	for (int l = 0; l < side_ID.Length(); l++)
	{
		/* whole side set */
		iArray2DT sideset_ALL = model_ALL.SideSet(side_ID[l]);
		const StringT& element_set_ID = model_ALL.SideSetGroupID(side_ID[l]);

		/* partition side set */
		iArray2DT sideset;
		if (sideset_ALL.MajorDim() > 0)
		{
			iArrayT elements_ALL(sideset_ALL.MajorDim());
			sideset_ALL.ColumnCopy(0, elements_ALL);
//			elements_ALL--;
			sideset_ALL.SetColumn(0, elements_ALL);
				
			iArrayT local_indices;
			partition.ReturnPartitionElements(element_set_ID, elements_ALL, local_indices);
			sideset.Dimension(local_indices.Length(), sideset_ALL.MinorDim());
			sideset.RowCollect(local_indices, sideset_ALL);
				
			/* map to (block) local numbering */
			if (sideset.MajorDim() > 0)
			{
				iArrayT elements(sideset.MajorDim());
				sideset.ColumnCopy(0, elements);
				partition.SetElementScope(PartitionT::kLocal, element_set_ID, elements);
//				elements++;
				sideset.SetColumn(0, elements);
			}
		}
			
		/* add */
		sideset++;
		int ss_id = atoi(side_ID[l]);
		int el_id = atoi(element_set_ID);
		if (model.PutSideSet(ss_id, el_id, sideset) != ModelFileT::kOK) ExceptionT::GeneralFail();
	}

	/* close database */
	model.Close();
}

/* basic MP support */
int FEExecutionManagerT::Rank(void) const { return fComm.Rank(); }
int FEExecutionManagerT::Size(void) const { return fComm.Size(); }
