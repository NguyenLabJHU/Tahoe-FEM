/* $Id: Aztec_fe.cpp,v 1.12 2005-04-05 16:07:07 paklein Exp $ */
/* created: paklein (08/01/1998) */
#include "Aztec_fe.h"

/* library support options */
#ifdef __AZTEC__

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>

#include "toolboxConstants.h"
#include "ExceptionT.h"
#include "az_aztec.h"
#include "ifstreamT.h"

#include "MSRBuilderT.h"
#include "iArray2DT.h"

using namespace Tahoe;

/* constructor */
Aztec_fe::Aztec_fe(const ParameterListT& parameters, ostream& msg, CommunicatorT& comm):
	AztecBaseT(msg, comm),
	fMSRBuilder(NULL),
	fMSRSet(0)
{
	/* keep non-default options and parameters */
	fAztecParams.TakeParameterList(parameters);

	/* construct MSR data builder */
	fMSRBuilder = new MSRBuilderT(false);
	if (!fMSRBuilder) throw ExceptionT::kOutOfMemory;
}

/* destructor */
Aztec_fe::~Aztec_fe(void) { delete fMSRBuilder; }

/* set solver options */
void Aztec_fe::SetAztecOptions(void)
{
	/* inherited (default settings) */
	AztecBaseT::SetAztecOptions();

	/* non-default options and parameters */
	fAztecParams.SetAztecOptions(options, params);
}

/* clear values in the matrix */
void Aztec_fe::Clear(void) { fval = 0.0; }

/* add to structure - active equations only (# > 0)
* NOTE: structure not set until #CALL#. equation data
* must persist outside of Aztec until (at least) then */
void Aztec_fe::AddEquationSet(const iArray2DT& eqnos)
{
	/* send to MSR data builder */
	fMSRBuilder->AddGroup(eqnos);
}

void Aztec_fe::AddEquationSet(const RaggedArray2DT<int>& eqnos)
{
	/* send to MSR data builder */
	fMSRBuilder->AddGroup(eqnos);
}

/* solve system based on data passed in for the rhs
* and return solution in result */
void Aztec_fe::Solve(dArrayT& rhs2result)
{
	/* general function */
	Solve(finitguess, rhs2result);
}

void Aztec_fe::Solve(const dArrayT& initguess, dArrayT& rhs2result)
{
	/* checks */
	if (!fMSRSet) throw ExceptionT::kGeneralFail;
	if (rhs2result.Length() != fupdate.Length()) throw ExceptionT::kSizeMismatch;

	/* allocate initial guess */
	finitguess.Dimension(InitGuessLength());
	
	/* set initial guess */	
	if (&initguess == &finitguess)
	{
		finitguess = 1.0; /* no initial guess */
	
		//TEMP - give steepest descent direction
		//       as the initial guess
		//finitguess.CopyPart(0, rhs2result, 0, rhs2result.Length());
		//finitguess.UnitVector();
		//finitguess *= -1.0;
	}
	else
	{
		/* check dimension */
		if (initguess.Length() != rhs2result.Length()) throw ExceptionT::kSizeMismatch;
	
		/* copy guess */
		finitguess.CopyPart(0, initguess, 0, initguess.Length());
	}

	/* set rhs */
	dArrayT rhs(RHSLength());
	rhs.CopyPart(0, rhs2result, 0, rhs2result.Length());

	/* Aztec solution driver */
	SolveDriver(rhs.Pointer(), finitguess.Pointer());
		
	/* copy result */
	rhs2result.CopyPart(0, finitguess, 0, rhs2result.Length());
		
	/* free memory from initial guess */
	finitguess.Free();
	
	/* check output status */
	cout << "\n number of iterations: " << int(status[AZ_its]) << '\n';
	cout <<   "   termination status: " << int(status[AZ_why]) << '\n';
	if (int(status[AZ_why]) != AZ_normal)
	  cout << "\n Aztec_fe::Solve: WARNING: exit status was not normal (0)\n" << endl;
	else
	  cout << endl;
}

/* statistics */
int Aztec_fe::NumNonZeroValues(void) const { return fval.Length() - 1; }

/*************************************************************************
* Private
*************************************************************************/

/* configure the update, bindx, and values array and return
* their length. if the columns indices for each row in bindx
* are sorted, is_sorted returns 1 and 0 otherwise */
int Aztec_fe::SetMSRData(int** update, int** bindx, double** val,
	int& is_sorted)
{
	/* set update vector - global numbering */
	fupdate.Dimension(N_update);
	int* pupdate = fupdate.Pointer();
	int n_update = Start_update; //OFFSET
	for (int i = 0; i < N_update; i++)
		(*pupdate++) = n_update++;

	/* set MSR structure data */
	fMSRBuilder->SetMSRData(fupdate, fbindx);
	is_sorted = 1; // MSRBuilderT sorts bindx data

	/* shift from equation numbers to global rows */
	fupdate--; //OFFSET

//NOTE - should be working with global equation numbers
#if 0	
	/* offset to global equation numbers */
	int offset = Start_update - 1;
	if (offset > 0)
	{
		fupdate += offset;
		for (int i = fbindx[0]; i < fbindx.Length(); i++)
			fbindx[i] += offset;
	}
#endif

//DEBUG
#if 0
cout << "\n Aztec_fe::SetMSRData: MSR data written to message file" << endl;
fMSRBuilder->WriteMSRData(fMessage, fupdate, fbindx);
#endif
	
	/* allocate the matrix and initialize to 0.0 */
	fval.Dimension(fbindx.Length());
	fval = 0.0;
	
	/* set pointers */
	*update = fupdate.Pointer();
	*bindx  = fbindx.Pointer();
	*val    = fval.Pointer();
	
	/* set flag */
	fMSRSet = 1;
	fMSRBuilder->ClearGroups();
	
	return fbindx.Length() - 1;
}

/* copy MSR data to RCV */
void Aztec_fe::GenerateRCV(iArrayT& r, iArrayT& c, dArrayT& v)
{
	/* overall dimension */
	int num_vals = fval.Length() - 1; // MSR format has 1 unused value
	r.Dimension(num_vals);
	c.Dimension(num_vals);
	v.Dimension(num_vals);

	/* start of off-diagonal data (MSR) */
	int*    pcol = fbindx.Pointer(N_update + 1);
	double* pval = fval.Pointer(N_update + 1);

	/* output rows in ascending column order */
	int shift = Start_update - 1; //OFFSET
	int count = 0;
	for (int row = 0; row < N_update; row++)
	{
		/* the diagonal */
		r[count] = row + shift;
		c[count] = row + shift;
		v[count] = fval[row];
		count++;

		int numvals = fbindx[row+1] - fbindx[row]; /* not incl. diagonal */
		for (int i = 0; i < numvals; i++)
		{
			r[count] = row + shift;
			c[count] = *pcol++;
			v[count] = *pval++;
			count++;
		}
	}
	
	/* check */
	if (count != num_vals)
	{
		cout << "\n SPOOLESMatrixT::GenerateRCV: translation error:\n"
		     <<   "   expected number of values = " << num_vals << '\n'
		     <<   "            number of values = " << count << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/* library support options */
#endif /* __AZTEC__ */
