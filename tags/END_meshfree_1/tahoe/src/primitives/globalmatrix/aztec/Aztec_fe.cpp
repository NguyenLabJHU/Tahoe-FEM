/* $Id: Aztec_fe.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (08/01/1998)                                          */
/* base class for interface for using Aztec with fe++                     */

#include "Aztec_fe.h"

/* library support options */
#ifdef __AZTEC__

#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>

#include "Constants.h"
#include "ExceptionCodes.h"
#include "az_aztec.h"

#include "fstreamT.h"
#include "MSRBuilderT.h"
#include "AztecReaderT.h"
#include "iArray2DT.h"

//DEBUG
#if 0
#include "StringT.h"
#include <fstream.h>
#endif

/* constructor */
Aztec_fe::Aztec_fe(ifstreamT& in):
	fMSRBuilder(NULL),
	fMSRSet(0)
{
	/* read non-default options and parameters */
	ReadOptionsParameters(in);

	/* construct MSR data builder */
	fMSRBuilder = new MSRBuilderT(false);
	if (!fMSRBuilder) throw eOutOfMemory;
}

/* destructor */
Aztec_fe::~Aztec_fe(void) { delete fMSRBuilder; }

/* set solver options */
void Aztec_fe::SetAztecOptions(void)
{
	/* inherited (default settings) */
	AztecBaseT::SetAztecOptions();

	/* non-default options */
	for (int i = 0; i < AZ_options.Length(); i++)
		options[AZ_options_dex[i]] = AZ_options[i];

	/* non-default parameters */
	for (int j = 0; j < AZ_params.Length(); j++)
		params[AZ_params_dex[j]] = AZ_params[j];
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
//DEBUG
#if 0
StringT name = "az.rcv";
name.Append(".p", proc_config[AZ_node]);
name.Append(".out");
ofstream out(name);
iArrayT r, c;
dArrayT v;
GenerateRCV(r, c, v);
for (int i = 0; i < r.Length(); i++)
	out << setw(kIntWidth) << r[i]
	    << setw(kIntWidth) << c[i]
		<< setw(kDoubleWidth) << v[i] << endl;
cout << "\n  Aztec_fe::Solve: EXIT" << endl;
throw;
#endif
//DEBUG

	/* checks */
	if (!fMSRSet) throw eGeneralFail;
	if (rhs2result.Length() != fupdate.Length()) throw eSizeMismatch;

	/* allocate initial guess */
	finitguess.Allocate(InitGuessLength());
	
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
		if (initguess.Length() != rhs2result.Length()) throw eSizeMismatch;
	
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
	{
		cout << "\n Aztec_fe::Solve: solver failed to converge" << endl;
		throw eBadJacobianDet;
	}
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
	fupdate.Allocate(N_update);
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
	
	/* allocate the matrix and initialize to 0.0 */
	fval.Allocate(fbindx.Length());
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

/* read (non-default) Aztec solver options and parameters */
void Aztec_fe::ReadOptionsParameters(ifstreamT& in)
{
	/* input interpreter */
	AztecReaderT az_reader;

	/* options */
	int num_options;
	in >> num_options;
	if (num_options < 0 || num_options > AZ_OPTIONS_SIZE)
		throw eBadInputValue;
//TEMP - num_options strictly can't be AZ_OPTIONS_SIZE, should
//       get this limit from AztecReaderT
	
	AZ_options_dex.Allocate(num_options);
	AZ_options.Allocate(num_options);
	for (int i = 0; i < num_options; i++)
	{
		/* read */
		int index, value;
		az_reader.ReadOption(in, index, value);
		
		/* store */
		AZ_options_dex[i] = index;
		AZ_options[i]     = value;
	}
	
	/* parameters */
	int num_params;
	in >> num_params;
	if (num_params < 0 || num_params > AZ_PARAMS_SIZE)
		throw eBadInputValue;
//TEMP - num_params strictly can't be AZ_PARAMS_SIZE, should
//       get this limit from AztecReaderT

	AZ_params_dex.Allocate(num_params);
	AZ_params.Allocate(num_params);
	for (int j = 0; j < num_params; j++)
	{
		/* read */
		int    index;
		double value;
		az_reader.ReadParameter(in, index, value);
		
		/* store */
		AZ_params_dex[j] = index;
		AZ_params[j]     = value;
	}
}

/* copy MSR data to RCV */
void Aztec_fe::GenerateRCV(iArrayT& r, iArrayT& c, dArrayT& v)
{
	/* overall dimension */
	int num_vals = fval.Length() - 1; // MSR format has 1 unused value
	r.Allocate(num_vals);
	c.Allocate(num_vals);
	v.Allocate(num_vals);

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
		throw eGeneralFail;
	}
}

/* library support options */
#endif /* __AZTEC__ */
