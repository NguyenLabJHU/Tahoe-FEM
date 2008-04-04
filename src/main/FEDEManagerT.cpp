
#include "FEDEManagerT.h"

// include whenever necessary
#include "TimeManagerT.h"

using namespace Tahoe;

/* constructor */
FEDEManagerT::FEDEManagerT(const StringT& input_file, ofstreamT& output, CommunicatorT& comm,
	const ArrayT<StringT>& argv, TaskT task):
	FEManagerT(input_file, output, comm, argv, task)
{
	SetName("tahoe_FEDE_coupling");
}

/* destructor */
FEDEManagerT::~FEDEManagerT(void)
{

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

		/* apply force from DE ghost particles */
		ApplyDEForce(); 
	
		/* solve the current time step */
		if (error == ExceptionT::kNoError) 
			error = SolveStep();
			
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

void FEDEManagerT::ApplyDEForce(void)
{
    ApplyNodeForce(8, 2, -10.0); // (node number, dof, force)
}

void FEDEManagerT::ApplyNodeForce(int node_num, int dof, double force)
{

}
