#include "FEDEManagerT.h"
#include "ModelManagerT.h"
#include "TimeManagerT.h"
#include "NodeManagerT.h"
#include "SolverT.h"

using namespace Tahoe;

/* constructor */
FEDEManagerT::FEDEManagerT(const StringT& input_file, ofstreamT& output, CommunicatorT& comm,
			   const ArrayT<StringT>& argv, TaskT task):
    FEManagerT(input_file, output, comm, argv, task)
{
    SetName("tahoe_DEM_coupling");

    dem::g_exceptioninf.open("dem_exception");
    if(! dem::g_exceptioninf) { cout<<"stream error!"<<endl; exit(-1); }
    dem::g_exceptioninf.setf(ios::scientific, ios::floatfield);
}

/* destructor */
FEDEManagerT::~FEDEManagerT(void)
{
    dem::g_exceptioninf.close();    
}

/* accept parameter list */
void FEDEManagerT::TakeParameterList(const ParameterListT& list)
{
    const char caller[] = "FEDEManagerT::TakeParameterList";
    
    /* inherited */
    FEManagerT::TakeParameterList(list);
    
    /* read DEM data */
    fDEManager.TakeParameter("dem.dat", fGhostElemSet);

    /* calculate parent domain coordinates for ghost particles */
    fDEManager.MapToParentDomain(ModelManager(), fGhostElemSet);
}

/* perform DEM simulation and calculate ghost forces */
void FEDEManagerT::DemComputeStep(void)
{

    if (fTimeManager->WriteOutput()) {

	/* print FE mesh info */
	fDEManager.PrintFE(NodeManager());

	/* print DE particle info */
	fDEManager.PrintDE();
    }

    /* perform DEM simulation*/
    fDEManager.Run();

    /* calculate ghost forces */
    fDEManager.GhostForce(fGhostFBC, ModelManager(), fGhostElemSet);
}

/* update ghost particles coordinates and print FE new mesh */
void FEDEManagerT::GhostDisplaceStep(void)
{
    fDEManager.GhostDisplace(NodeManager(), fGhostElemSet);
}

void FEDEManagerT::Solve(void)
{
    const char caller[] = "FEDEManagerT::Solve";
    
    /* set to initial condition */
    ExceptionT::CodeT error = InitialCondition();
    
    /* loop over time increments */
    while (error == ExceptionT::kNoError && fTimeManager->Step())
    {
	/* initialize the current time step */
	if (error == ExceptionT::kNoError) 
	    error = InitStep();
	
	/* perform DEM simulation and calculate ghost forces */
	DemComputeStep();
	
	/* solve the current time step */
	if (error == ExceptionT::kNoError) 
	    error = SolveStep();

	/* update ghost particles coordinates based on FE mesh deformation */
	GhostDisplaceStep();
	
	/* close the current time step */
	if (error == ExceptionT::kNoError) 
	    error = CloseStep();
	
	/* handle errors */
	switch (error)
	{
	    case ExceptionT::kNoError:
		/* nothing to do */
		break;
	    case ExceptionT::kGeneralFail:
	    case ExceptionT::kBadJacobianDet:
	    {
		cout << '\n' << caller << ": trying to recover from error: " << ExceptionT::ToString(error) << endl;
		
		/* reset system configuration */
		error = ResetStep();
		
		/* cut time step */
		if (error == ExceptionT::kNoError)
		    if (!DecreaseLoadStep())
			error = ExceptionT::kGeneralFail;
		
		break;
	    }
	    default:
		cout << '\n' << caller <<  ": no recovery for error: " << ExceptionT::ToString(error) << endl;
	}
    }
}


ExceptionT::CodeT FEDEManagerT::SolveStep(void)
{
    ExceptionT::CodeT error = ExceptionT::kNoError;
    try {
	
	/* check group flag */
	if (fCurrentGroup != -1) throw ExceptionT::kGeneralFail;
	
	int loop_count = 0;
	bool all_pass = false;
	SolverT::SolutionStatusT status = SolverT::kContinue;
	
	/* clear status */
	fSolverPhasesStatus = 0;
	
	while (!all_pass && 
	       (fMaxSolverLoops == -1 || loop_count < fMaxSolverLoops) &&
	       status != SolverT::kFailed)
	{
	    
	    /* one solver after the next */
	    all_pass = true;
	    for (int i = 0; status != SolverT::kFailed && i < fSolverPhases.MajorDim(); i++)
	    {
		/* group parameters */
		fCurrentGroup = fSolverPhases(i,0);
		int iter = fSolverPhases(i,1);
		int pass = fSolverPhases(i,2);
		
		/* call solver */
		status = fSolvers[fCurrentGroup]->Solve(iter, *this, fGhostFBC);
		
		/* check result */
		fSolverPhasesStatus(i, kGroup) = fCurrentGroup;
		int last_iter = fSolverPhasesStatus(i, kIteration);
		fSolverPhasesStatus(i, kIteration) = fSolvers[fCurrentGroup]->IterationNumber();
		if (status == SolverT::kFailed) 
		{
		    all_pass = false;
		    fSolverPhasesStatus(i, kPass) = -1;					
		}
		else if (status == SolverT::kConverged && 
			 (pass == -1 || (fSolverPhasesStatus(i, kIteration) - last_iter) <= pass)
			 || fSolverPhasesStatus(i, kIteration) > 1000)
		{
		    all_pass = all_pass && true; /* must all be true */
		    fSolverPhasesStatus(i, kPass) = 1;
		}
		else
		{
		    all_pass = false;
		    fSolverPhasesStatus(i, kPass) = 0;
		}
	    }
	    
	    /* count */
	    loop_count++;
	    
	    /* write status */
	    if (fSolverPhases.MajorDim() > 1)
	    {
		cout << "\n Solver status: pass " << loop_count << '\n';
		cout << setw(kIntWidth) << "phase"
		     << setw(kIntWidth) << "solver"
		     << setw(kIntWidth) << "its."
		     << setw(kIntWidth) << "pass" << '\n';
		fSolverPhasesStatus.WriteNumbered(cout);
		cout << endl;
	    }			
	}
	
	/* solver looping not successful */
	if (!all_pass) {
	    cout << "\n Solver loop not converged after " << loop_count << " iterations" << endl;
	    error = ExceptionT::kGeneralFail;
	} else if (fSolverPhases.MajorDim() > 1) {
	    cout << "\n Solver loop converged in " << loop_count << " iterations" << endl;
	}
    }	
    catch (ExceptionT::CodeT exc) {
	cout << "\n FEManagerT::SolveStep: caught exception: " 
	     << ExceptionT::ToString(exc) << endl;
	error = exc;
    }
    
    /* reset group flag */
    fCurrentGroup = -1;
    
    /* done */
    return error;
}

void FEDEManagerT::FormRHS(int group, ArrayT<FBC_CardT>& fGhostFBC)
{
    /* state */
    SetStatus(GlobalT::kFormRHS);
    
    /* unlock assembly into RHS */
    fSolvers[group]->UnlockRHS();
    
    /* nodal force contribution - F(t) */
    fNodeManager->FormRHS(group, fGhostFBC);
    
    /* lock assembly into RHS */
    fSolvers[group]->LockRHS();
}
