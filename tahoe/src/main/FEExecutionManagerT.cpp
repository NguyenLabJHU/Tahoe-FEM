/* $Id: FEExecutionManagerT.cpp,v 1.41.2.15 2003-06-11 18:56:03 hspark Exp $ */
/* created: paklein (09/21/1997) */
#include "FEExecutionManagerT.h"

#include <iostream.h>
#include <time.h>
#include <ctype.h>

#include "fstreamT.h"
#include "Environment.h"
#include "toolboxConstants.h"
#include "ExceptionT.h"

#if defined(__MWERKS__) && __option(profile)
#include <Profiler.h>
#endif

#include "fstreamT.h"
#include "FEManagerT.h"
#include "FEManagerT_mpi.h"
#include "IOManager_mpi.h"
#include "ModelManagerT.h"
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

/* needed for bridging calculations FEExecutionManagerT::RunBridging */
#ifdef BRIDGING_ELEMENT
#include "FEManagerT_bridging.h"
#include "FEManagerT_THK.h"
#include "TimeManagerT.h"
#include "NodeManagerT.h"
#include "dSPMatrixT.h"
#include "FieldT.h"
#include "IntegratorT.h"
#include "ElementBaseT.h"
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
}

/**********************************************************************
 * Protected
 **********************************************************************/

/* MUST be overloaded */
void FEExecutionManagerT::RunJob(ifstreamT& in, ostream& status)
{
	const char caller[] = "FEExecutionManagerT::RunJob";

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
				RunJob_serial(in, status);
			} else {
				cout << "\n RunJob_parallel: " << in.filename() << endl;
				RunJob_parallel(in, status);
			}
			break;
		}
		case kDecompose:
		{
			cout << "\n RunDecomp_serial: " << in.filename() << endl;
			if (fComm.Size() > 1) cout << " RunDecomp_serial: SERIAL ONLY" << endl;
				
			/* decompose using rank 0 */
			if (fComm.Rank() == 0) RunDecomp_serial(in, status);

			/* synch */
			fComm.Barrier();
			break;
		}
		case kJoin:
		{
			cout << "\n RunJoin_serial: " << in.filename() << endl;
			if (fComm.Size() > 1) cout << " RunJoin_serial: SERIAL ONLY" << endl;
			
			/* join using rank 0 */
			if (fComm.Rank() == 0) RunJoin_serial(in, status);

			/* synch */
			fComm.Barrier();
			break;
		}
		case kBridging:
		{
			cout << "\n RunBridging: " << in.filename() << endl;
			if (fComm.Size() > 1) ExceptionT::GeneralFail(caller, "RunBridging for SERIAL ONLY");

			RunBridging(in, status);
			break;
		}
		case kTHK:
		{
			cout << "\n RunTHK: " << in.filename() << endl;
			if (fComm.Size() > 1) ExceptionT::GeneralFail(caller, "RunTHK for SERIAL ONLY");

			RunTHK(in, status);
			break;
		}
		default:
			ExceptionT::GeneralFail("FEExecutionManagerT::RunJob", "unknown mode: %d", mode);
	}
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

/* Tahoe bridging calculation */
void FEExecutionManagerT::RunBridging(ifstreamT& in, ostream& status) const
{
	const char caller[] = "FEExecutionManagerT::RunBridging";
	StringT path;
	path.FilePath(in.filename());

	/* read atomistic and continuum source files */
	StringT atom_file, bridge_atom_file;
	in >> atom_file
	   >> bridge_atom_file;
	   
	StringT continuum_file, bridge_continuum_file;
	in >> continuum_file
	   >> bridge_continuum_file;

	/* streams for atomistic Tahoe */
	atom_file.ToNativePathName();
	atom_file.Prepend(path);
	ifstreamT atom_in('#', atom_file);
	if (atom_in.is_open())
		cout << " atomistic parameters file: " << atom_file << endl;
	else
		ExceptionT::BadInputValue(caller, "file not found: %s", atom_file.Pointer());
	StringT atom_out_file;
	atom_out_file.Root(atom_in.filename());
	atom_out_file.Append(".out");
	ofstreamT atom_out;
	atom_out.open(atom_out_file);

	bridge_atom_file.ToNativePathName();
	bridge_atom_file.Prepend(path);
	ifstreamT bridge_atom_in('#', bridge_atom_file);
	if (bridge_atom_in.is_open())
		cout << " atomistic bridging parameters file: " << bridge_atom_file << endl;
	else
		ExceptionT::BadInputValue(caller, "file not found: %s", bridge_atom_file.Pointer());

	/* streams for continuum Tahoe */
	continuum_file.ToNativePathName();
	continuum_file.Prepend(path);
	ifstreamT continuum_in('#', continuum_file);
	if (continuum_in.is_open())
		cout << " continuum parameters file: " << continuum_file << endl;
	else
		ExceptionT::BadInputValue(caller, "file not found: %s", continuum_file.Pointer());
	StringT continuum_out_file;
	continuum_out_file.Root(continuum_in.filename());
	continuum_out_file.Append(".out");
	ofstreamT continuum_out;
	continuum_out.open(continuum_out_file);

	bridge_continuum_file.ToNativePathName();
	bridge_continuum_file.Prepend(path);
	ifstreamT bridge_continuum_in('#', bridge_continuum_file);
	if (bridge_continuum_in.is_open())
		cout << " continuum bridging parameters file: " << bridge_continuum_file << endl;
	else
		ExceptionT::BadInputValue(caller, "file not found: %s", bridge_continuum_file.Pointer());

	/* brigding solver log file */
	StringT log_file;
	log_file.Root(in.filename());
	log_file.Append(".log");
	ofstreamT log_out;
	log_out.open(log_file);

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
		char job_char;
		atom_in >> job_char;
		
		/* initialize FEManager_THK using atom values */
		FEManagerT_THK atoms(atom_in, atom_out, fComm, bridge_atom_in);
		atoms.Initialize();
	    	continuum_in >> job_char;
		FEManagerT_bridging continuum(continuum_in, continuum_out, fComm, bridge_continuum_in);
		continuum.Initialize();
		t1 = clock();
		phase = 1;
                
		/* split here depending on whether integrators are explicit or implicit
		 * check only one integrator assuming they both are the same */
		const IntegratorT* mdintegrate = atoms.Integrator(0);
		IntegratorT::ImpExpFlagT impexp = mdintegrate->ImplicitExplicit();
       
		if (impexp == IntegratorT::kImplicit)
			RunStaticBridging(continuum, atoms, log_out);
		else if (impexp == IntegratorT::kExplicit)
			RunDynamicBridging(continuum, atoms, log_out);
		else
			ExceptionT::GeneralFail(caller, "unknown integrator type %d", impexp);
		
	
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
			status << " solution phase.";
		
		/* fix clock values */
		if (t1 == 0) t1 = clock();
		if (t2 == 0) t2 = clock();		

		atom_out << endl;
		continuum_out << endl;
	}
	/* stop day/date info */
	time_t stoptime;
	time(&stoptime);

	/* output timing */
	status << "\n     Filename: " << in.filename() << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   " Construction: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "     Solution: " << double(t2 - t1)/CLOCKS_PER_SEC << " sec.\n";
	status << "    Stop time: " << ctime(&stoptime);
	status << "\n End Execution\n" << endl;
}

void FEExecutionManagerT::RunStaticBridging(FEManagerT_bridging& continuum, FEManagerT_THK& atoms, ofstream& log_out) const
{
	const char caller[] = "FEExecutionManagerT::RunStaticBridging";
	    
	/* configure ghost nodes */
	int group = 0;
	int order1 = 0;
	StringT bridging_field = "displacement";
	bool active = false;
	atoms.InitGhostNodes();
	continuum.InitInterpolation(atoms.GhostNodes(), bridging_field, *atoms.NodeManager());
	continuum.InitProjection(atoms.NonGhostNodes(), bridging_field, *atoms.NodeManager(), active);

#if 0
	/* cross coupling matricies */
	int neq_A = atoms.NodeManager()->Field(bridging_field)->NumEquations();
	int neq_C = continuum.NodeManager()->Field(bridging_field)->NumEquations();
	dSPMatrixT K_AC(neq_A, neq_C, 0), K_G_NG;
	dSPMatrixT K_CA(neq_C, neq_A, 0), G_Interpolation;
	dArrayT F_A(neq_A), F_C(neq_C);
	continuum.InterpolationMatrix(bridging_field, G_Interpolation);
#endif

	/* time managers */
	TimeManagerT* atom_time = atoms.TimeManager();
	TimeManagerT* continuum_time = continuum.TimeManager();
 
	dArray2DT field_at_ghosts;
	atom_time->Top();
	continuum_time->Top();
	int d_width = OutputWidth(log_out, field_at_ghosts.Pointer());
	while (atom_time->NextSequence() && continuum_time->NextSequence())
	{	
		/* set to initial condition */
		atoms.InitialCondition();
		continuum.InitialCondition();

		/* loop over time increments */
		AutoArrayT<int> loop_count, atom_iter_count, continuum_iter_count;
		bool seq_OK = true;
		while (seq_OK && 
			atom_time->Step() &&
			continuum_time->Step()) //TEMP - same clock
		{
			log_out << "\n Step = " << atom_time->StepNumber() << '\n'
				<< " Time = " << atom_time->Time() << endl;
			
			/* running status flag */
			ExceptionT::CodeT error = ExceptionT::kNoError;		

			/* initialize step */
			if (error == ExceptionT::kNoError) error = atoms.InitStep();
			if (error == ExceptionT::kNoError) error = continuum.InitStep();
			
			/* solver phase status */
			const iArray2DT& atom_phase_status = atoms.SolverPhasesStatus();
			const iArray2DT& continuum_phase_status = continuum.SolverPhasesStatus();

#if 0
			/* set cross-coupling */
			atoms.Form_G_NG_Stiffness(bridging_field, K_G_NG);
			K_AC.MultAB(K_G_NG, G_Interpolation);
			K_CA.Transpose(K_AC);
#endif

			/* loop until both solved */
			int group_num = 0;
			double atoms_res, continuum_res, combined_res_0 = 0.0;
			int count = 0;
			int atom_last_iter, atom_iter, continuum_last_iter, continuum_iter;
			atom_last_iter = atom_iter = continuum_last_iter = continuum_iter = 0;
			while (count == 0 || (atom_iter > 0 || continuum_iter > 0)) //TEMP - assume just one phase

//			while (1 || error == ExceptionT::kNoError &&
//			(atom_phase_status(0, FEManagerT::kIteration) > 0 ||
//			continuum_phase_status(0, FEManagerT::kIteration) > 0)) //TEMP - assume just one phase
			{
				count++;

				/* solve atoms */
				if (1 || error == ExceptionT::kNoError) {
					atoms.ResetCumulativeUpdate(group);
					error = atoms.SolveStep();
				}

#if 0
				/* set cross-coupling */
				atoms.Form_G_NG_Stiffness(bridging_field, K_G_NG);
				K_AC.MultAB(K_G_NG, G_Interpolation);
				K_CA.Transpose(K_AC);
#endif
					
				/* apply solution to continuum */
				continuum.ProjectField(bridging_field, *atoms.NodeManager(), order1);
#if 0
				K_CA.Multx(atoms.CumulativeUpdate(group_num), F_C);
				F_C *= -1.0;
				continuum.SetExternalForce(group_num, F_C);
#endif
				continuum.FormRHS(group_num);
				continuum_res = continuum.Residual(group_num).Magnitude(); //serial
					
				/* solve continuum */
				if (1 || error == ExceptionT::kNoError) {
					continuum.ResetCumulativeUpdate(group);
					error = continuum.SolveStep();
				}
				
				/* apply solution to atoms */
				int order = 0;  // displacement only for static case
				continuum.InterpolateField(bridging_field, order, field_at_ghosts);
				atoms.SetFieldValues(bridging_field, atoms.GhostNodes(), order, field_at_ghosts);
#if 0
				K_AC.Multx(continuum.CumulativeUpdate(group_num), F_A);
				F_A *= -1.0;
				atoms.SetExternalForce(group_num, F_A);
#endif
				atoms.FormRHS(group_num);
				atoms_res = atoms.Residual(group_num).Magnitude(); //serial

				/* reset the reference errors */
				if (count == 1) {
					combined_res_0 = atoms_res + continuum_res;
					atoms.SetReferenceError(group_num, combined_res_0);
					continuum.SetReferenceError(group_num, combined_res_0);
				}
					
				/* log residual */
				double tot_rel_error = (fabs(combined_res_0) > kSmall) ? 
					(atoms_res + continuum_res)/combined_res_0 : 0.0;
				    log_out << setw(kIntWidth) << count << ": "
					<< setw(d_width) << atoms_res << " (A) | "
					<< setw(d_width) << continuum_res << " (C) | "
					<< setw(d_width) << tot_rel_error << endl;

				/* number of interations in last pass */
				int atom_total_iter = atom_phase_status(0, FEManagerT::kIteration);
				int continuum_total_iter = continuum_phase_status(0, FEManagerT::kIteration);
				atom_iter = atom_total_iter - atom_last_iter;
				continuum_iter = continuum_total_iter - continuum_last_iter;
				atom_last_iter = atom_total_iter;
				continuum_last_iter = continuum_total_iter;
			}
				
			loop_count.Append(count);
			atom_iter_count.Append(atom_last_iter);
			continuum_iter_count.Append(continuum_last_iter);
			
			/* close step */
			if (1 || error == ExceptionT::kNoError) error = atoms.CloseStep();
			if (1 || error == ExceptionT::kNoError) error = continuum.CloseStep();

			/* check for error */
			if (0)
//		if (error != ExceptionT::kNoError)
				ExceptionT::GeneralFail(caller, "hit error %d", error);
			//TEMP - no error recovery yet
		}
			
		cout << "\n Number of bridging iterations:\n";
		cout << setw(kIntWidth) << "step" 
			<< setw(kIntWidth) << "cycles" 
			<< setw(kIntWidth) << "a-its."
			<< setw(kIntWidth) << "c-its."<< '\n';
		for (int i = 0; i < loop_count.Length(); i++)
			cout << setw(kIntWidth) << i+1
				<< setw(kIntWidth) << loop_count[i]
				<< setw(kIntWidth) << atom_iter_count[i]
				<< setw(kIntWidth) << continuum_iter_count[i] << '\n';
	}
}

void FEExecutionManagerT::RunDynamicBridging(FEManagerT_bridging& continuum, FEManagerT_THK& atoms, ofstream& log_out) const
{
	const char caller[] = "FEExecutionManagerT::RunDynamicBridging";

	/* configure ghost nodes */
	int group = 0;
	int order1 = 0;	// For InterpolateField, 3 calls to obtain displacement/velocity/acceleration
	int order2 = 1;
	int order3 = 2;
	dArray2DT field_at_ghosts, totalu, fubig, fu, projectedu, boundghostdisp, boundghostvel, boundghostacc;
	dArray2DT thkforce, gaussdisp;
	dSPMatrixT ntf;
	iArrayT activefenodes;
	StringT bridging_field = "displacement";
	atoms.InitGhostNodes();
	bool makeinactive = false;	
	
	/* figure out boundary atoms for use with THK boundary conditions, 
	   ghost atoms for usage with MD force calculations */
	const iArrayT& boundaryghostatoms = atoms.InterpolationNodes();
	
	int numgatoms = (atoms.GhostNodes()).Length();	// total number of ghost atoms
	int numbatoms = boundaryghostatoms.Length() - numgatoms;	// total number of boundary atoms
	dArray2DT gadisp(numgatoms,2), gavel(numgatoms,2), gaacc(numgatoms,2);
	dArray2DT badisp(numbatoms,2), bavel(numbatoms,2), baacc(numbatoms,2);
	iArrayT allatoms(boundaryghostatoms.Length()), gatoms(numgatoms), batoms(numbatoms);
	allatoms.SetValueToPosition();
	batoms.CopyPart(0, allatoms, numgatoms, numbatoms);
	gatoms.CopyPart(0, allatoms, 0, numgatoms);
	continuum.InitInterpolation(boundaryghostatoms, bridging_field, *atoms.NodeManager());
	continuum.InitProjection(atoms.NonGhostNodes(), bridging_field, *atoms.NodeManager(), makeinactive);
	
	/* time managers */
	TimeManagerT* atom_time = atoms.TimeManager();
	TimeManagerT* continuum_time = continuum.TimeManager();

	atom_time->Top();
	continuum_time->Top();
	int d_width = OutputWidth(log_out, field_at_ghosts.Pointer());
	while (atom_time->NextSequence() && continuum_time->NextSequence())
	{	
		/* set to initial condition */
		atoms.InitialCondition();
	
		/* calculate fine scale part of MD displacement and total displacement u */
		continuum.InitialProject(bridging_field, *atoms.NodeManager(), projectedu, order1);
	
		/* solve for initial FEM force f(u) as function of fine scale + FEM */
		/* use projected totalu instead of totalu for initial FEM displacements */
		fubig = InternalForce(projectedu, atoms);
		
		/* calculate global interpolation matrix ntf */
		continuum.Ntf(ntf, atoms.NonGhostNodes(), activefenodes);
		//cout << "ntf = " << ntf << endl;
		
		/* compute FEM RHS force as product of ntf and fu */
		dArrayT fx(ntf.Rows()), fy(ntf.Rows()), tempx(ntf.Cols()), tempy(ntf.Cols());
		dArray2DT ntfproduct(ntf.Rows(), 2);
		fu.Dimension((atoms.NonGhostNodes()).Length(), 2);
		fu.RowCollect(atoms.NonGhostNodes(), fubig);
		fu.ColumnCopy(0,tempx);
		ntf.Multx(tempx,fx);
		ntfproduct.SetColumn(0,fx);
		fu.ColumnCopy(1,tempy);
		ntf.Multx(tempy,fy);
		ntfproduct.SetColumn(1,fy);
		
		/* Add FEM RHS force to RHS using SetExternalForce */
		continuum.SetExternalForce(bridging_field, ntfproduct, activefenodes);
		
		/* now d0, v0 and a0 are known after InitialCondition */
		continuum.InitialCondition();

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

		/* Write interpolated FEM values at MD ghost nodes into MD field - displacement only */
		atoms.SetFieldValues(bridging_field, atoms.GhostNodes(), order1, gadisp);
		
		/* store initial MD boundary displacement histories */
		thkforce = atoms.THKForce(badisp);

		/* figure out timestep ratio between fem and md simulations */
		int nfesteps = continuum_time->NumberOfSteps();
		double mddt = atom_time->TimeStep();
		double fedt = continuum_time->TimeStep();
		double d_ratio = fedt/mddt;		
		int ratio = int((2.0*d_ratio + 1.0)/2.0);
		
		/* running status flag */
		ExceptionT::CodeT error = ExceptionT::kNoError;	
		
		/* now loop over continuum and atoms after initialization */
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
	
				/* Write interpolated FEM values at MD ghost nodes into MD field - displacements only */
				atoms.SetFieldValues(bridging_field, atoms.GhostNodes(), order1, gadisp);				
			       	/* calculate THK force on boundary atoms, update displacement histories */
				thkforce = atoms.THKForce(badisp);
																			
				/* solve MD equations of motion */
				if (1 || error == ExceptionT::kNoError) {
						atoms.ResetCumulativeUpdate(group);
						error = atoms.SolveStep();
				}
				
				FieldT* cmd = atoms.NodeManager()->Field(bridging_field);
				dArray2DT cacc = (*cmd)[2];
				//cout << "md acc = " << cacc << endl;
				
				/* close  md step */
				if (1 || error == ExceptionT::kNoError) error = atoms.CloseStep();    
			}
			continuum_time->Step();
			
			/* initialize step */
			if (1 || error == ExceptionT::kNoError) error = continuum.InitStep();
            
			/* calculate total displacement u = FE + fine scale here using updated FEM displacement */
			continuum.BridgingFields(bridging_field, *atoms.NodeManager(), *continuum.NodeManager(), totalu);
			
			/* calculate FE internal force as function of total displacement u here */
			fubig = InternalForce(totalu, atoms);
			fu.RowCollect(atoms.NonGhostNodes(), fubig);
			fu.ColumnCopy(0,tempx);
			ntf.Multx(tempx,fx);
			ntfproduct.SetColumn(0,fx);
			fu.ColumnCopy(1,tempy);
			ntf.Multx(tempy,fy);
			ntfproduct.SetColumn(1,fy);	// SetExternalForce updated via pointer
			//cout << "fem force = " << ntfproduct << endl;
			
			/* solve FE equation of motion using internal force just calculated */
			if (1 || error == ExceptionT::kNoError) {
					continuum.ResetCumulativeUpdate(group);
					error = continuum.SolveStep();
			}
                 
			FieldT* fem = continuum.NodeManager()->Field(bridging_field);
			dArray2DT facc1 = (*fem)[2];
			//cout << "fem acc = " << facc1 << endl;
		
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
			
			/* Write interpolated FEM values at MD ghost nodes into MD field - displacement only */
			atoms.SetFieldValues(bridging_field, atoms.GhostNodes(), order1, gadisp);
			
			/* close fe step */
			if (1 || error == ExceptionT::kNoError) error = continuum.CloseStep();
                        
		}

		/* check for error */
		if (0)
			ExceptionT::GeneralFail(caller, "hit error %d", error);
                
	}
}

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
		FEManagerT_THK thk(in, out, fComm, dummy_bridging_input);
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
#endif

/* standard serial driver */
void FEExecutionManagerT::RunJob_serial(ifstreamT& in,
	ostream& status) const
{
	/* output stream */
	StringT outfilename;
	outfilename.Root(in.filename());
	outfilename.Append(".out");
	ofstreamT out;
	out.open(outfilename);

#ifdef __MWERKS__
	if (!out.is_open())
	{
		cout << "\n FEExecutionManagerT::RunJob_serial: could not open file: "
		     << outfilename << endl;
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
		FEManagerT analysis1(in, out, fComm);
		analysis1.Initialize();

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

/* generate decomposition files */
void FEExecutionManagerT::RunDecomp_serial(ifstreamT& in, ostream& status) const
{
	int size = 0;
	int index;
	if (!CommandLineOption("-decomp", index)) 
		ExceptionT::GeneralFail();
	else {
	
		/* look for size */
		if (fCommandLineOptions.Length() > index+1)
		{
			const char* opt = fCommandLineOptions[index+1];
			if (strlen(opt) > 1 && isdigit(opt[1]))
				size = atoi(opt+1); /* opt[0] = '-' */
		}
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
	if (size == 0) 
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

	/* set stream comment marker */
	in.set_marker('#');

	/* time markers */
	clock_t t0 = 0, t1 = 0;

	/* start day/date info */
	time_t starttime;
	time(&starttime);
	try
	{
		t0 = clock();
	
		/* get the model file name */
		StringT model_file, suffix;
		IOBaseT::FileTypeT format;
		GetModelFile(in, model_file, format);

		/* output map file */
		StringT map_file;
		map_file.Root(in.filename());
		map_file.Append(".n", size);
		map_file.Append(".io.map");

		/* set output map and and generate decomposition */
		Decompose(in, size, method, model_file, format, map_file);
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
	status << "\n     Filename: " << in.filename() << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   "Decomposition: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "    Stop time: " << ctime(&stoptime);	
	status << "\n End Execution\n" << endl;
}

/* join parallel results files */
void FEExecutionManagerT::RunJoin_serial(ifstreamT& in, ostream& status) const
{
	/* set stream comment marker */
	in.set_marker('#');

	/* time markers */
	clock_t t0 = 0, t1 = 0;

	/* start day/date info */
	time_t starttime;
	time(&starttime);
	try
	{
		t0 = clock();

		/* to read file parameters */
		ofstreamT out;
		FEManagerT fe_man(in, out, fComm);
		fe_man.Initialize(FEManagerT::kParametersOnly);
		
		/* model file parameters */
		StringT model_file = fe_man.ModelManager()->DatabaseName();
		IOBaseT::FileTypeT model_format = fe_man.ModelManager()->DatabaseFormat();
		IOBaseT::FileTypeT results_format = fe_man.OutputFormat();
		if (results_format == IOBaseT::kTahoe ||
		    results_format == IOBaseT::kTahoeII)
			results_format = IOBaseT::kTahoeResults;

		int size = fComm.Size();
		int index;
		if (!CommandLineOption("-join", index)) 
			ExceptionT::GeneralFail();
		else {
			if (fCommandLineOptions.Length() > index+1)
			{
				const char* opt = fCommandLineOptions[index+1];
				if (strlen(opt) > 1 && isdigit(opt[1]))
					size = atoi(opt+1);
			}
		}

		/* prompt for decomp size */
		if (size == 1)
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
		StringT version = fe_man.Version();
		StringT title   = fe_man.Title();
		StringT input   = in.filename();
		OutputBaseT* output = IOBaseT::NewOutput(program, version, title, input, 
			results_format, cout);
		
		/* construct joiner */
		JoinOutputT output_joiner(in.filename(), model_file, model_format, 
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
	status << "\n     Filename: " << in.filename() << '\n';
	status <<   "   Start time: " << ctime(&starttime);
	status <<   "         Join: " << double(t1 - t0)/CLOCKS_PER_SEC << " sec.\n";
	status <<   "    Stop time: " << ctime(&stoptime);	
	status << "\n End Execution\n" << endl;
}

/* testing for distributed execution */
void FEExecutionManagerT::RunJob_parallel(ifstreamT& in, ostream& status) const
{
	const char caller[] = "::RunJob_parallel";

	/* set stream comment marker */
	in.set_marker('#');

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
	out_file.Root(in.filename());
	out_file.Append(".p", rank);
	out_file.Append(".out");
	ofstreamT out;
	out.open(out_file);

	int token; // for run time check sums
	try {
	t0 = clock();
	
	/* get the model file name */
	StringT model_file, suffix;
	IOBaseT::FileTypeT format;
	GetModelFile(in, model_file, format);
	
	/* output map file */
	StringT map_file;
	map_file.Root(in.filename());
	map_file.Append(".n", size);
	map_file.Append(".io.map");

	/* set output map and generate decomposition */
	token = 1;
	if (rank == 0)
	{
		/* look for method */
		int index = 0;
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

		if (NeedDecomposition(model_file, size) || (!CommandLineOption("-split_io") && NeedOutputMap(in, map_file, size)))
		  {
		/* prompt if not found */
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
	
		/* run decomp */
		try { Decompose(in, size, method, model_file, format, map_file); }
		catch (ExceptionT::CodeT code)
		{
			cout << "\n " << caller << ": exception on decomposition: " << code << endl;
			token = 0;
		}
		  }
	}

	/* synch and check status */
	if (fComm.Sum(token) != size) ExceptionT::GeneralFail();

	/* read partition information */
	PartitionT partition;
	StringT part_file;
	part_file.Root(model_file);
	part_file.Append(".n", size);
	part_file.Append(".part", rank);
	ifstreamT part_in(in.comment_marker(), part_file);
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
	StringT partial_file;
	suffix.Suffix(model_file);
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
	ifstreamT in_loc(in.comment_marker(), in.filename());
	if (!fJobCharPutBack)
	{
		char filetypechar;
		in_loc >> filetypechar;
	}
	FEManagerT_mpi FEman(in_loc, out, fComm, &partition, FEManagerT_mpi::kRun);
	try { FEman.Initialize(); }
	catch (ExceptionT::CodeT code)
	{
		status << "\n \"" << in_loc.filename() << "\" exit on exception " << code << " during the\n";
		status << " construction phase. Check the input file for errors." << endl;
		
		/* echo some lines from input */
		if (code == ExceptionT::kBadInputValue) Rewind(in_loc, status);
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
		/* read output map */
		iArrayT output_map;
		ReadOutputMap(in, map_file, output_map);

		/* set-up local IO */
		IOMan = new IOManager_mpi(in, fComm, output_map, *(FEman.OutputManager()), FEman.Partition(), model_file, format);
		if (!IOMan) throw ExceptionT::kOutOfMemory;
		
		/* set external IO */
		FEman.SetExternalIO(IOMan);
		}
		
		catch (ExceptionT::CodeT code)
		{
			token = 0;
			status << "\n \"" << in.filename() << "\" exit on exception " << code 
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
		status << "\n \"" << in.filename() << "\" exit on exception " << code << " during the\n";
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
	status << "\n     Filename: " << in.filename() << '\n';
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

/* extract the model file name from the stream */
void FEExecutionManagerT::GetModelFile(ifstreamT& in, StringT& model_file,
	IOBaseT::FileTypeT& format) const
{
	/* partially construct FE manager */
	ifstreamT in_temp(in.comment_marker(), in.filename());
	if (!fJobCharPutBack)
	{
		char filetypechar;
		in_temp >> filetypechar;
	}

	ofstreamT out;
	FEManagerT fe_temp(in_temp, out, fComm);
	fe_temp.Initialize(FEManagerT::kParametersOnly);

	ModelManagerT* model = fe_temp.ModelManager();
	format = model->DatabaseFormat();
	model_file = model->DatabaseName();
}

void FEExecutionManagerT::Decompose(ifstreamT& in, int size,
	int decomp_type, const StringT& model_file, IOBaseT::FileTypeT format, 
	const StringT& output_map_file) const
{	
	/* dispatch */
	switch (decomp_type)
	{
		case PartitionT::kGraph:
			Decompose_graph(in, size, model_file, format, output_map_file);
			break;

		case PartitionT::kAtom:
			Decompose_atom(in, size, model_file, format, output_map_file);
			break;
			
		case PartitionT::kSpatial:
			cout << "\n FEExecutionManagerT::Decompose: spatial decomposition not implemented yet" << endl;
			break;
						
		default:
			cout << "\n FEExecutionManagerT::Decompose: unrecognized method: " << decomp_type << endl;
	}
}

void FEExecutionManagerT::Decompose_atom(ifstreamT& in, int size,
	const StringT& model_file, IOBaseT::FileTypeT format, 
	const StringT& output_map_file) const
{
#pragma unused(in)
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

	/* output map file? - not needed for PartitionT::DecompTypeT == kAtom */
#pragma unused(output_map_file)
}

void FEExecutionManagerT::Decompose_spatial(ifstreamT& in, int size,
	const StringT& model_file, IOBaseT::FileTypeT format, 
	const StringT& output_map_file) const
{
#pragma unused(in)
#pragma unused(size)
#pragma unused(model_file)
#pragma unused(format)
#pragma unused(output_map_file)
	cout << "\n FEExecutionManagerT::Decompose_spatial: not implemented" << endl;
}

/* graph-based decomposition */
void FEExecutionManagerT::Decompose_graph(ifstreamT& in, int size,
	const StringT& model_file, IOBaseT::FileTypeT format, 
	const StringT& output_map_file) const
{
	bool split_io = CommandLineOption("-split_io");
	bool need_output_map = NeedOutputMap(in, output_map_file, size) && !split_io;
	bool need_decomp = NeedDecomposition(model_file, size);
	if (need_output_map || need_decomp)
	{
		/* echo stream */
		StringT decomp_file;
		decomp_file.Root(in.filename());
		decomp_file.Append(".out");
		ofstreamT decomp_out;
		decomp_out.open(decomp_file);
	
		ifstreamT in_decomp(in.comment_marker(), in.filename());
		if (!fJobCharPutBack)
		{
			char filetypechar;
			in_decomp >> filetypechar;
		}

		/* construct global problem */
		FEManagerT_mpi global_FEman(in_decomp, decomp_out, fComm, NULL, FEManagerT_mpi::kDecompose);
		try { global_FEman.Initialize(FEManagerT::kAllButSolver); }
		catch (ExceptionT::CodeT code)
		{
			cout << "\n FEExecutionManagerT::Decompose: exception during construction: "
		         << ExceptionT::ToString(code) << endl;
			throw code;
		}

		/* set output map */
		if (need_output_map)
		{
			cout << "\n Generating output map: " << output_map_file << endl;
			IOManager* global_IOman = global_FEman.OutputManager();
			iArrayT output_map;
			SetOutputMap(global_IOman->ElementSets(), output_map, size);
	
			/* write map file */
			ofstreamT map_out(output_map_file);
			map_out << "# number of processors\n";
			map_out << size << '\n';
			map_out << "# number of output sets\n";
			map_out << output_map.Length() << '\n';
			map_out << "# set to processor output map\n";
			map_out << output_map.wrap(8) << '\n';
			cout << " Generating output map: " << output_map_file << ": DONE"<< endl;
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
			try
			{
				cout << "\n Decomposing: " << model_file << endl;
				global_FEman.Decompose(partition, graph, true, method);
				cout << " Decomposing: " << model_file << ": DONE"<< endl;
			}
			catch (ExceptionT::CodeT code)
			{
				cout << "\n FEExecutionManagerT::Decompose: exception during decomposition: "
			         << code << endl;
				throw code;
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
				dArrayT inex(nnd); inex = 0;
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
					ExceptionT::GeneralFail("FEExecutionManagerT::Decompose", 
						"unexpected node number shift: %d", shift);

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
				ifstreamT part_in(in.comment_marker(), part_file);
				
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
				ExceptionT::GeneralFail("FEExecutionManagerT::Decompose_spatial", 
					"error opening file: %s", (const char*) model_file);
	
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
					try { EchoPartialGeometry(partition[i], model_ALL, partial_file, format); }
					catch (ExceptionT::CodeT error)
					{
						cout << "\n ::Decompose: exception writing file: " << partial_file << endl;
						throw error;
					}
				}
			}
		}
	}
	else
		cout << "\n ::Decompose: decomposition files exist" << endl;
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

/* returns true if a new output map is needed */
bool FEExecutionManagerT::NeedOutputMap(ifstreamT& in, const StringT& map_file,
	int size) const
{
	/* open map file */
	ifstreamT map_in(in.comment_marker(), map_file);
	if (map_in.is_open())
	{
		int num_parts;
		map_in >> num_parts;
		if (num_parts != size) return true;
	
		int num_sets;
		map_in >> num_sets;
		iArrayT map(num_sets);
		map_in >> map;
		
		int min, max;
		map.MinMax(min, max);
		if (min < 0 || max >= size)
			return true;
		else
			return false;
	}
	else
		return true;
}

void FEExecutionManagerT::ReadOutputMap(ifstreamT& in, const StringT& map_file,
	iArrayT& map) const
{
	/* map file */
	ifstreamT map_in(in.comment_marker(), map_file);
	if (!map_in.is_open())
		ExceptionT::GeneralFail("FEExecutionManagerT::ReadOutputMap", 
			"could not open io map file: %s", (const char*) map_file);

	/* read map */
	int size, num_sets;
	map_in >> size >> num_sets;
	map.Dimension(num_sets);
	map_in >> map;

	/* check */
	int min, max;
	map.MinMax(min, max);
	if (min < 0 || max >= size)
	{
		cout << "\n FEExecutionManagerT::ReadOutputMap: map error\n";
		cout << map.wrap(5) << '\n';
		cout.flush();
		ExceptionT::GeneralFail();
	}
}

/* set output map based on length of map */
void FEExecutionManagerT::SetOutputMap(const ArrayT<OutputSetT*>& output_sets,
	iArrayT& output_map, int size) const
{
	/* initialize map - processor number for each element block */
	output_map.Dimension(output_sets.Length());
	output_map = 0;
	
	iArrayT output_counts(size);
	output_counts = 0;
	
	/* output set size */
	iArrayT set_size(output_sets.Length());
	for (int i = 0; i < output_sets.Length(); i++)
	{
		const OutputSetT& set = *(output_sets[i]);
	
		int size = 0;
		size += set.NumNodes()*set.NodeOutputLabels().Length();
		size += set.NumElements()*set.ElementOutputLabels().Length();
		
		set_size[i] = size;
	}

	/* map output sets to processors */
	for (int j = 0; j < output_sets.Length(); j++)
	{
		int max_at;
		set_size.Max(max_at);
	
		int min_at;
		output_counts.Min(min_at);
	
		output_map[max_at] = min_at;
		output_counts[min_at] += set_size[max_at];
		set_size[max_at] = 0;
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
