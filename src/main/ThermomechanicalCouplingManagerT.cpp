/* $Id $ */
#include "ThermomechanicalCouplingManagerT.h"

#ifdef BRIDGING_ELEMENT

#include "ifstreamT.h"
#include "SolverT.h"
#include "FEManagerT_THK.h"
#include "NodeManagerT.h"
#include "OutputSetT.h"
#include "TimeManagerT.h"
#include "ParticlePairT.h"
#include "FieldT.h"
#include "dSPMatrixT.h"
#include "ParameterContainerT.h"
#include "ParameterUtils.h"
#include "CommunicatorT.h"
#include "CommManagerT.h"
#include "IOManager.h"
#include "ModelManagerT.h"

using namespace Tahoe;

const double fkB = 0.00008617385; // Boltzmann's constant in [eV/K], see ScaledVelocityNodesT.cpp


/* constructor */
ThermomechanicalCouplingManagerT::ThermomechanicalCouplingManagerT(const StringT& input_file, ofstreamT& output, CommunicatorT& comm,
	const ArrayT<StringT>& argv, TaskT task):
	MultiManagerT(input_file, output, comm, argv, task)
{
	SetName("tahoe_thermomechanical_coupling");
	
}

ThermomechanicalCouplingManagerT::~ThermomechanicalCouplingManagerT(void)
{
	/* clean up communicators */
	if (fCoarseComm != &fComm) delete fCoarseComm;
	if (fFineComm != &fComm) delete fFineComm;
}
	
/* (re-)set system to initial conditions */
ExceptionT::CodeT ThermomechanicalCouplingManagerT::InitialCondition(void)
{
  	return MultiManagerT::InitialCondition();
}

/* (re-)set the equation number for the given group */
void ThermomechanicalCouplingManagerT::Solve(void)
{
	const char caller[] = "ThermomechanicalCouplingManagerT::Solve";

	NodeManagerT& fine_node_manager = *(fFine->NodeManager());
	int nsd = fine_node_manager.NumSD();
	int nat = fine_node_manager.NumNodes();
	cout << " number of atoms : " << nat << "\n";
	NodeManagerT& coarse_node_manager = *(fCoarse->NodeManager());
	int nnd = coarse_node_manager.NumNodes();

	/* time managers */
	TimeManagerT* atom_time = fFine->TimeManager();
	TimeManagerT* continuum_time = fCoarse->TimeManager();

	int nfesteps = continuum_time->NumberOfSteps();
	double mddt = atom_time->TimeStep();
	double fedt = continuum_time->TimeStep();
	double d_ratio = fedt/mddt;		
	cout << "time step ratio : " << d_ratio 
	<< ", continuum: " << fedt << ", atom: " << mddt << "\n";
	int ratio = int((2.0*d_ratio + 1.0)/2.0);

	/* set to initial condition */
	ExceptionT::CodeT error = InitialCondition();

	dArrayT ke(nat);
	dArray2DT  ave_kinetic_energy(nat,1), projected_temperature(nnd,1);

	int sample_interval = 100;
	int num_samples = 0;

	/* loop over time increments */
	for (int i = 0; i < nfesteps; i++)	
	{
		ave_kinetic_energy = 0.0;
		num_samples = 0;
		for (int j = 0; j < ratio; j++)	// MD update first
		{
			atom_time->Step();	
								
			/* initialize step */
			if (error == ExceptionT::kNoError) 
				error = fFine->InitStep();
				
			/* solve the current time step */
			if (error == ExceptionT::kNoError) 
				error = fFine->SolveStep();
			
			/* close the current time step */
			if (error == ExceptionT::kNoError)
				error = fFine->CloseStep();

// use WriteOutput()
        		//if (((j % sample_interval) == 0)
        		if ((atom_time->WriteOutput())
			&& (fFineField->Order() > 0)) // got velocities!
        		{
				num_samples++;
				fParticles->AtomicKineticEnergies(ke);
				cout << num_samples << " adding velocities @ " << atom_time->StepNumber() << " " <<  atom_time->Time() << " T : " 
<< ke.Average()/( 0.5*nsd*fkB ) << ", KE : " << ke.Sum() << 	"\n";
                		double *akep = ave_kinetic_energy.Pointer();
                		double *kep = ke.Pointer();
                		for (int ii = 0; ii < ave_kinetic_energy.MajorDim(); ii++,kep++,akep++) 
					*akep += *kep; 
			}
		}

		//continuum_time->Step();
		//fCoarse->ProjectField(fFineField->FieldName(), *fFine->NodeManager(), 0); // hard set "order" to 0
		ave_kinetic_energy /= 0.5*nsd*fkB*num_samples; // average
		fCoarse->Project(ave_kinetic_energy, projected_temperature);

#if 1
                /* screen dump */
		cout << "time : " << atom_time->Time()  ; 
                cout << " average kinetic energy \n";
		const dArray2DT& fine_coords = fine_node_manager.InitialCoordinates();
                double *ptr = ave_kinetic_energy.Pointer();
                for (int i = 0; i < nat; i++, ptr++) { // assume real atoms 
                cout << i << " " ;
		for (int ii = 0; ii < nsd; ii++)  {
			cout << fine_coords(i,ii) << " " ;
		}
		cout <<*ptr<<"\n"; }

		cout << "time : " << atom_time->Time()  ; 
                cout << " projection of ave. kinetic energy \n";
		const dArray2DT& coarse_coords = coarse_node_manager.InitialCoordinates();
                double *ppt = projected_temperature.Pointer();
                for (int i = 0; i < projected_temperature.MajorDim(); i++, ppt++) {
                cout << i << " " ;
		for (int ii = 0; ii < nsd; ii++)  {
			cout << coarse_coords(i,ii) << " " ;
		}
		cout <<*ppt<<"\n"; }
#endif

		dArray2DT e_values;
		fCoarse->WriteOutput(fCoarseOutputID, projected_temperature, e_values);
	}
}

#if 0
/* describe the parameters needed by the interface */
void ThermomechanicalCouplingManagerT::DefineParameters(ParameterListT& list) const
{
	        /* inherited */
	        ParameterInterfaceT::DefineParameters(list);

		        ParameterT beta(fBeta, "beta");
			        beta.AddLimit(0.0, LimitT::LowerInclusive);
				        list.AddParameter(beta);
}
#endif


/* accept parameter list */
void ThermomechanicalCouplingManagerT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "ThermomechanicalCouplingManagerT::TakeParameterList";

	/* inherited */
//	MultiManagerT::TakeParameterList(list);

	/* inherited - don't call direct base class method */
	ParameterInterfaceT::TakeParameterList(list);

	/* path to parameters file */
	StringT path;
	path.FilePath(fInputFile);
	TaskT task = kRun;

	/* parse/validate continuum input */
	StringT continuum_input = list.GetParameter("continuum_input");
	continuum_input.ToNativePathName();
	continuum_input.Prepend(path);
	ParameterListT continuum_params;
	ParseInput(continuum_input, continuum_params, true, true, true, fArgv);
			
	/* construct continuum solver */
	if (fCoarseComm->Size() != 1)
		ExceptionT::GeneralFail(caller, "parallel execution error");
	if (Size() > 1) /* change file name so output files are unique */  {
		StringT suffix;
		suffix.Suffix(continuum_input);
		continuum_input.Root();
		continuum_input.Append(".p", Rank());
		continuum_input.Append(suffix);
	}
	StringT continuum_output_file;
	continuum_output_file.Root(continuum_input);
	continuum_output_file.Append(".out");
	fCoarseOut.open(continuum_output_file);
	fCoarse = TB_DYNAMIC_CAST(FEManagerT_bridging*, FEManagerT::New(continuum_params.Name(), continuum_input, fCoarseOut, *fCoarseComm, fArgv, task));
	if (!fCoarse) ExceptionT::GeneralFail(caller, "could not construct continuum solver");
	fCoarse->TakeParameterList(continuum_params);

	/* parse/validate atomistic input */
	StringT atom_input = list.GetParameter("atom_input");
	atom_input.ToNativePathName();
	atom_input.Prepend(path);
	ParameterListT atom_params;
	ParseInput(atom_input, atom_params, true, true, true, fArgv);

	/* construct atomistic solver */
	if (Size() != fFineComm->Size())
		ExceptionT::GeneralFail(caller, "parallel execution error");
	StringT atom_output_file;
	atom_output_file.Root(atom_input);
	if (Size() > 1) atom_output_file.Append(".p", Rank());
	atom_output_file.Append(".out");
	fFineOut.open(atom_output_file);
	fFine = TB_DYNAMIC_CAST(FEManagerT_bridging*, FEManagerT::New(atom_params.Name(), atom_input, fFineOut, *fFineComm, fArgv, task));
	if (!fFine) ExceptionT::GeneralFail(caller, "could not construct atomistic solver");
	fFine->TakeParameterList(atom_params);

	/* check consistency between time managers */
	TimeManagerT* atom_time = fFine->TimeManager();
	TimeManagerT* continuum_time = fCoarse->TimeManager();

	/* use parameters from coarse scale solver */
	fTimeManager = fCoarse->TimeManager();
	fOutputFormat = fCoarse->OutputFormat();

	/* don't compute initial conditions */
	fFine->SetComputeInitialCondition(false);
	fCoarse->SetComputeInitialCondition(false);

	/* resolve bridging fields */
	const StringT& bridging_field = list.GetParameter("bridging_field");
	fFineField = fFine->NodeManager()->Field(bridging_field);
	if (!fFineField) ExceptionT::GeneralFail(caller, "could not resolve fine scale \"%s\" field", bridging_field.Pointer());
	fCoarseField = fCoarse->NodeManager()->Field(bridging_field);
	if (!fFineField) ExceptionT::GeneralFail(caller, "could not resolve coarse scale \"%s\" field", bridging_field.Pointer());

	/* resolve integrator types */
//	if (fFineField->Integrator().ImplicitExplicit() != fCoarseField->Integrator().ImplicitExplicit())
//		ExceptionT::GeneralFail(caller, "time integrator mismatch");
	fImpExp = fFineField->Integrator().ImplicitExplicit();

	/* collect the ghost atom ID list */
	ArrayT<StringT> ghost_atom_ID;
	const ParameterListT* ghosts = list.List("ghost_atom_ID_list");
	if (ghosts)	StringListT::Extract(*ghosts, ghost_atom_ID);

	/* configure projection/interpolation */
	NodeManagerT& fine_node_manager = *(fFine->NodeManager());	
	int group = 0;
	int order1 = 0;
	bool make_inactive = true;
	fFine->InitGhostNodes(fFineField->FieldName(), ghost_atom_ID, fCoarse->ProjectImagePoints());
	fCoarse->InitInterpolation(fFineField->FieldName(), fFine->GhostNodes(), fine_node_manager.InitialCoordinates());
	fCoarse->InitProjection(fFineField->FieldName(), *(fFine->CommManager()), fFine->NonGhostNodes(), fine_node_manager, make_inactive);

	/* send coarse/fine output through the fFine output */
	int ndof = fFine->NodeManager()->NumDOF(group);
	ArrayT<StringT> labels(2*ndof);
	const char* coarse_labels[] = {"UC_X", "UC_Y", "UC_Z"};
	const char* fine_labels[] = {"UF_X", "UF_Y", "UF_Z"};
	int dex = 0;
	for (int i = 0; i < ndof; i++) labels[dex++] = coarse_labels[i];
	for (int i = 0; i < ndof; i++) labels[dex++] = fine_labels[i];
	const iArrayT& non_ghost_nodes = fFine->NonGhostNodes();
	fAtomConnectivities.Alias(non_ghost_nodes.Length(), 1, non_ghost_nodes.Pointer());
	OutputSetT output_set(GeometryT::kPoint, fAtomConnectivities, labels, false);
	fOutputID = fFine->RegisterOutput(output_set);

	/* construct solver */
	int n1 = fFine->NumGroups();
	int n2 = fCoarse->NumGroups();
//	if (n1 != n2) ExceptionT::GeneralFail(caller, "number of groups must match: %d != %d", n1, n2);
	const ParameterListT* multi_solver = list.ListChoice(*this, "multi_solver_choice");
	if (multi_solver)
	{
		/* construct */
		SolverT* solver = SolverT::New(*this, multi_solver->Name(), 0);
		if (!solver) ExceptionT::GeneralFail(caller, "could not construct solver \"%s\"",
			multi_solver->Name().Pointer());
		solver->TakeParameterList(*multi_solver);
		
		/* store */
		fSolvers.Dimension(1);
		fSolvers[0] = solver;
	}
	SetSolver();

	/* default solver phases */
	fMaxSolverLoops = 1;
	fSolverPhases.Dimension(1, 3);
	fSolverPhases(0,0) = 0;
	fSolverPhases(0,1) =-1;
	fSolverPhases(0,2) =-1;
	fSolverPhasesStatus.Dimension(fSolverPhases.MajorDim(), kNumStatusFlags);
	fSolverPhasesStatus = 0;

	/* terms to include in the equilibrium equations */
	fFineToCoarse = list.GetParameter("fine_to_coarse");
	fCoarseToFine = list.GetParameter("coarse_to_fine");
	if (fCoarseToFine) /* enforce zero bond density in projected cells */
		fCoarse->DeactivateFollowerCells();

	/* needed to solve overlap */
	const dArray2DT& fine_init_coords = fine_node_manager.InitialCoordinates();
	const ParticlePairT* particle_pair = fFine->ParticlePair();
	if (!particle_pair) ExceptionT::GeneralFail(caller, "could not resolve ParticlePairT");
	
	/* overlap correction method */
	StringT overlap_path;
	const ParameterT* overlap_file = list.Parameter("overlap_file");
	if (overlap_file) {
		overlap_path = *overlap_file;
		overlap_path.ToNativePathName();

		StringT path;
		path.FilePath(fInputFile);
		overlap_path.Prepend(path);
	} /* default name */
	else {
		overlap_path.Root(fInputFile);
		overlap_path.Append(".overlap");
	}
	
	const ParameterListT& overlap = list.GetListChoice(*this, "overlap_correction_choice");
	if (overlap.Name() == "prescribe_overlap")
	{
		double density = overlap.GetParameter("bond_density");
		cout << "\n " << caller << ": \"prescribe_overlap\" not implemented. p = 1.0" << endl;
	}
	else if (overlap.Name() == "by-bond_multiplier")
	{
		/* extract parameters */
		double smoothing = overlap.GetParameter("smoothing");
		double bind_1 = overlap.GetParameter("bind_to_1.0");
		double reg = overlap.GetParameter("constraint_regularization");
		int nip = overlap.GetParameter("density_nip");
		double init_bound_width = overlap.GetParameter("init_bound_width");
		
		/* compute overlap correction */
		double bound_0 = init_bound_width/2.0;
		fCoarse->CorrectOverlap_2(particle_pair->Neighbors(), fine_init_coords, 
			overlap_path, smoothing, bind_1, reg, bound_0, nip);
	}
	else if (overlap.Name() == "by-bond_penalty")
	{
		/* extract parameters */
		double smoothing = overlap.GetParameter("smoothing");
		double bind_1 = overlap.GetParameter("bind_to_1.0");
		double bound_tol = overlap.GetParameter("bound_tolerance");
		double stiffness_jump = overlap.GetParameter("stiffness_jump");
		int nip = overlap.GetParameter("density_nip");
		
		/* compute overlap correction */
		fCoarse->CorrectOverlap_22(particle_pair->Neighbors(), fine_init_coords, 
			overlap_path, smoothing, bind_1, bound_tol, stiffness_jump, nip);
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized overlap correction method \"%s\"",
			overlap.Name().Pointer());



        /* output for kinetic temperature  */
	ModelManagerT* model = fCoarse->ModelManager();
        ArrayT<StringT> vlabels(1);
        const char* label[] = {"kT"};
        vlabels[0] = label[0];
	/* specify output - "free set" */
	const char* id = "1";
	StringT ID ; // HACK
	ID = fCoarse->ElementGroup(0)->ElementBlockID(1);
        OutputSetT coarse_output_set(model->ElementGroupGeometry(ID), model->ElementGroup(ID), vlabels);
        fCoarseOutputID = fCoarse->RegisterOutput(coarse_output_set);

	/* find particle class */
        for (int i = 0; i < fFine->NumElementGroups(); i++)
        {
           ElementBaseT* element_base = fFine->ElementGroup(i); 
           ParticleT* particle = dynamic_cast<ParticleT*>(element_base);
           if (particle) { fParticles = particle;   }
        }
}

/**********************************************************************
 * Protected
 **********************************************************************/

#endif /* BRIDGING_ELEMENT */
