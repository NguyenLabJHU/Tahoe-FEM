/* $Id: AztecBaseT.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (07/28/1998)                                          */
/* Base class for Aztec iterative solver                                  */

#include "AztecBaseT.h"

/* library support options */
#ifdef __AZTEC__

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream.h>
#include <iomanip.h>

#include "ExceptionCodes.h"
#include "Constants.h"
#include "az_aztec.h"

#ifdef __NEW_THROWS__
#include <new.h>
#endif

/* constructor */
AztecBaseT::AztecBaseT(void): N_update(0), update_index(NULL),
	update_bin(NULL), srow_dex(NULL), srow_val(NULL), external(NULL),
	extern_index(NULL), rpntr(NULL), cpntr(NULL), indx(NULL), bpntr(NULL),
	data_org(NULL),
	bindx_transform(NULL)
{
	/* allocate parameter arrays */
	proc_config = new int[AZ_PROC_SIZE];
	options     = new int[AZ_OPTIONS_SIZE];
	params      = new double[AZ_PARAMS_SIZE];
	status      = new double[AZ_STATUS_SIZE];
	
	/* check */
	if (!proc_config || !options || !params || !status)
		throw eOutOfMemory;
		
	/* get number of processors and the name of this processor */
	AZ_processor_info(proc_config);
}

/* destructor */
AztecBaseT::~AztecBaseT(void)
{
	/* free parameter arrays */
	delete[] proc_config;
	delete[] options;
	delete[] params;
	delete[] status;

	/* allocated by Aztec.lib */
	FreeAztecMemory();
	
	/* free quick find data */
	delete[] update_bin;
	delete[] srow_dex;
	delete[] srow_val;
	
	/* bindx in transormed format */
	if (bindx_transform != bindx)
		delete[] bindx_transform;
}

/* set solver options */
void AztecBaseT::SetAztecOptions(void)
{
	/* choose among AZTEC options (see User's Guide) */
	AZ_defaults(options, params);

	options[AZ_solver]   = AZ_cgs;
	options[AZ_scaling]  = AZ_none;
	options[AZ_precond]  = AZ_ls;
	options[AZ_conv]     = AZ_r0;
	options[AZ_output]   = 1;
	options[AZ_pre_calc] = AZ_calc;
	options[AZ_max_iter] = 500;
	options[AZ_poly_ord] = 5;
	options[AZ_overlap]  = AZ_none;
	options[AZ_kspace]   = 30;
	options[AZ_orthog]   = AZ_modified;
	options[AZ_aux_vec]  = AZ_resid;

	params[AZ_tol]       = 4.00e-9;
	params[AZ_drop]      = 0.0;
}

void AztecBaseT::WriteAztecOptions(ostream& out) const
{
	out << "\n Aztec options (see \"az_aztec_defs.h\"):\n";
	out << " Solver type . . . . . . . . . . . . . . . . . . = "  << options[AZ_solver]   << '\n';
	out << " Scaling option. . . . . . . . . . . . . . . . . = "  << options[AZ_scaling]  << '\n';
	out << " Preconditioner. . . . . . . . . . . . . . . . . = "  << options[AZ_precond]  << '\n';
	out << " Residual expression . . . . . . . . . . . . . . = "  << options[AZ_conv]     << '\n';
	out << " Output information flags/print increment. . . . = "  << options[AZ_output]   << '\n';
	out << " Factorization re-use option . . . . . . . . . . = "  << options[AZ_pre_calc] << '\n';
	out << " Maximum number of iterations. . . . . . . . . . = "  << options[AZ_max_iter] << '\n';
	out << " Polynomial order/steps for precondtioner. . . . = "  << options[AZ_poly_ord] << '\n';
	out << " Overlap option for decomposition algorithms . . = "  << options[AZ_overlap]  << '\n';
	out << " GMRES Krylov subspace size. . . . . . . . . . . = "  << options[AZ_kspace]   << '\n';
	out << " GMRES orthogonalization scheme. . . . . . . . . = "  << options[AZ_orthog]   << '\n';
	out << " Auxiliary vector option . . . . . . . . . . . . = "  << options[AZ_aux_vec]  << '\n';

	out << "\n Aztec parameters (see \"az_aztec_defs.h\"):\n";
	out << " Convergence tolerance . . . . . . . . . . . . . = "  << params[AZ_tol]       << '\n';
	out << " LU preconditioner drop tolerance. . . . . . . . = "  << params[AZ_drop]      << endl;
}

/* set Aztec data structures - calls SetMSRData() to set the
* MSR data and then initializes other needed arrays, sets and
* verifies solver options */
void AztecBaseT::Initialize(int num_eq, int start_eq)
{
	Start_update = start_eq;
	N_update = num_eq;

	/* checks */
	if (Start_update < 1) throw eBadInputValue;
	if (N_update < 1) throw eBadInputValue;

	/* free existing memory */
	FreeAztecMemory();

	/* derived classes set MSR data */
	int is_sorted;
	int numterms = SetMSRData(&update, &bindx, &val, is_sorted);
	
	/* bindx data must be sorted */
	if (!is_sorted) Sort_bindx();

	/* keep copy */
	if (bindx_transform != bindx) delete[] bindx_transform;
	if (proc_config[AZ_N_procs] == 1)
		bindx_transform = bindx;
	else
	{
#ifdef __NEW_THROWS__
		try { bindx_transform = new int[numterms+1]; }
		catch (bad_alloc) { bindx_transform = NULL; }
#else
		bindx_transform = new int[numterms+1];
#endif
		if (!bindx_transform)
		{
			cout << "\n AztecBaseT::Initialize: out of memory" << endl;
			throw eOutOfMemory;
		}
		memcpy(bindx_transform, bindx, (numterms+1)*sizeof(int));
	}

	/* convert matrix to a local distributed MSR */
	AZ_transform(proc_config, &external, bindx_transform, val, update, &update_index,
		&extern_index, &data_org, N_update, indx, bpntr, rpntr, &cpntr,
		AZ_MSR_MATRIX);

	/* free memory not used for serial execution */
	if (proc_config[AZ_N_procs] == 1) FreeAztec_MP_Memory();

	/* initial quick find data */
	SetUpQuickFind();

	/* set solver options */
	SetAztecOptions();

	/* verify Aztec parameters */
	int error_code = AZ_check_input(data_org, options, params, proc_config);
	if (error_code)
	{
		AZ_print_error(error_code);
		throw eGeneralFail;
	}	
}

/* assemble row values into global structure using the global column
* indices given in coldex - status is 1 if successful, 0 otherwise */
void AztecBaseT::AssembleRow(int row, int numvals, const int* col_dex,
	const double* rowvals, int& status)
{
	/* find row index */
	int rowkey = AZ_quick_find(row, update, N_update, QF_shift, update_bin);
	if (rowkey == -1)
		status = 0;
	else
	{			
		/* copy */
		memcpy(srow_dex, col_dex, numvals*sizeof(int));
		memcpy(srow_val, rowvals, numvals*sizeof(double));
	
		/* sort copies */
		AZ_sort(srow_dex, numvals, NULL, srow_val);
		
		/* assemble */
		int*    pasmdex = srow_dex;
		double* pasmval = srow_val;

		int  rowlength = bindx[rowkey+1] - bindx[rowkey] + 1; // incl. diagonal
		int*      pdex = bindx + bindx[rowkey];
		double*   pval = val   + bindx[rowkey];
			
		int rowcount = 0;
		int valcount = 0;
		while (rowcount < rowlength && valcount < numvals)
		{
			/* diagonal value */
			if (*pasmdex == row)
			{
				/* assemble */
				val[rowkey] += *pasmval;
				
				/* increment assembled */
				pasmdex++;
				pasmval++;
				valcount++;
			}
			/* matching index */
			else if (*pasmdex == *pdex)
			{
				/* assemble */
				*pval += *pasmval;
					
				/* increment assembled */
				pasmdex++;
				pasmval++;
				valcount++;
			}
			/* no matching index */
			else
			{
				/* increment matrix */
				pdex++;
				pval++;
				rowcount++;
			}
		}
			
		/* check completion */
		if (rowcount == rowlength)
			status = 0;
		else
			status = 1;
	}
}

/* assemble diagonal values into the global structure
* status is 1 if successful, 0 otherwise */
void AztecBaseT::AssembleDiagonals(int numvals, const int* rows,
	const double* vals, int& status)
{
	/* assemble values */
	status = 1;
	for (int i = 0; i < numvals && status; i++)
	{
		int rowkey = AZ_quick_find(*rows, update, N_update, QF_shift, update_bin);
		if (rowkey == -1)
			status = 0;
		else	
		{
			/* assemble */
			val[rowkey] += *vals;		
		
			/* next */
			rows++;
			vals++;
		}		
	}
}

/* write non-zero values to stream as {row,col,value} */
void AztecBaseT::PrintNonZero(ostream& out) const
{
	/* off-diagonal data */
	int*    pcol = bindx + N_update + 1;
	double* pval = val   + N_update + 1;

	/* output rows in ascending column order */
	for (int row = 0; row < N_update; row++)
	{
		int numvals    = bindx[row+1] - bindx[row]; /* inc. diagonal */
		int beforediag = 1;
		int valcount   = 0;
		
		while (valcount < numvals || beforediag)
		{
			/* output diagonal value */
			if (beforediag && (*pcol > row || valcount == numvals))
			{
				out << setw(kIntWidth)    << row + 1;
				out << setw(kIntWidth)    << row + 1;
				out << setw(kDoubleWidth) << val[row] << '\n';
				beforediag = 0;
			}
			/* off-diagonal values */
			else
			{
				out << setw(kIntWidth)    << row + 1;
				out << setw(kIntWidth)    << *pcol + 1;
				out << setw(kDoubleWidth) << *pval << '\n';
				
				/* increment */
				pcol++;
				pval++;
				valcount++;
			}
		}
	}
}

/*************************************************************************
* Protected
*************************************************************************/

/* return correct dimension */
int AztecBaseT::InitGuessLength(void) const
{
	return data_org[AZ_N_internal] +
		   data_org[AZ_N_border  ] +
		   data_org[AZ_N_external];
}

int AztecBaseT::RHSLength(void) const
{
	return data_org[AZ_N_internal] +
	       data_org[AZ_N_border  ];
}

/* check parameters and solver with given rhs and initguess */
void AztecBaseT::SolveDriver(double* rhs, double* initguess)
{
	/* Aztec solution driver */
AZ_solve(initguess, rhs, options, params, indx, bindx_transform,
	rpntr, cpntr, bpntr, val, data_org, status, proc_config);

	/* free internal memory */
	AZ_free_memory(AZ_SYS);
}

/*************************************************************************
* Private
*************************************************************************/

/* allocate memory for quick find */
void AztecBaseT::SetUpQuickFind(void)
{
	/* quick find bin */
	delete[] update_bin;
	update_bin = new int[2 + (N_update + 4)/4]; /* oversize */
	if (!update_bin) throw eOutOfMemory;

	/* initialize shift and bin */
	AZ_init_quick_find(update, N_update, &QF_shift, update_bin);

	/* find max row length */
	int maxlength = 0;
	for (int i = 0; i < N_update; i++)
	{
		int rowln = bindx[i+1] - bindx[i];
		maxlength = (rowln > maxlength) ? rowln : maxlength;
	}
	
	/* add space for the diagonal */
	maxlength += 1;

	/* allocate space for sorted row data */
	srow_dex = new int[maxlength];
	if (!srow_dex) throw eOutOfMemory;
	
	srow_val = new double[maxlength];
	if (!srow_val) throw eOutOfMemory;
}

/* free memory allocated by Aztec.lib */
void AztecBaseT::FreeAztecMemory(void)
{	
	/* allocated by Aztec using calloc or malloc */
free((void *) indx);
	free((void *) rpntr);
	free((void *) cpntr);
	free((void *) bpntr);
	free((void *) data_org);

	/* practice safe set */
	indx         = NULL;
	rpntr        = NULL;
	cpntr        = NULL;
	bpntr        = NULL;
	data_org     = NULL;

	/* memory for MP execution */
	FreeAztec_MP_Memory();
}

/* free memory which is only used in multi-processor execution */
void AztecBaseT::FreeAztec_MP_Memory(void)
{
	/* allocated by Aztec using calloc or malloc */
free((void *) update_index);
free((void *) external);
free((void *) extern_index);

	/* practice safe set */
	update_index = NULL;
	external     = NULL;
	extern_index = NULL;
}

/* sort active columns in bindx into ascending order */
void AztecBaseT::Sort_bindx(void)
{
	for (int i = 0; i < N_update; i++)
	{
		/* data length */
		int numvals = bindx[i+1] - bindx[i];

		/* sort copies */
		AZ_sort(bindx + bindx[i], numvals, NULL, NULL);
	}
}

/* library support options */
#endif  /* __AZTEC__ */
