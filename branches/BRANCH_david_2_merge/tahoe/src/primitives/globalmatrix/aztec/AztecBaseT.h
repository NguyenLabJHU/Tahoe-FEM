/* $Id: AztecBaseT.h,v 1.5 2003-02-28 02:07:17 paklein Exp $ */
/* created: paklein (07/28/1998) */

#ifndef _AZTEC_BASE_T_H_
#define _AZTEC_BASE_T_H_

#include "Environment.h"

/* library support options */
#ifdef __AZTEC__

/* forward declarations */
#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class CommunicatorT;

/** low-level interface to the Aztec iterative solver */
class AztecBaseT
{
public:

	/** constuctor */
	AztecBaseT(ostream& msg, CommunicatorT& comm);
	
	/* destructor */
	virtual ~AztecBaseT(void);

	/* set structure - n_update is number of local update rows */
	void Initialize(int num_eq, int start_eq);

	/* set solver options */
	virtual void SetAztecOptions(void);
	void WriteAztecOptions(ostream& out) const;

	/* assemble row values into global structure using the global column
	 * indices given in coldex - status is 1 if successful, 0 otherwise */
	void AssembleRow(int row, int numvals, const int* col_dex,
		const double* rowvals, int& status);

	/* assemble diagonal values into the global structure
	 * status is 1 if successful, 0 otherwise */
	void AssembleDiagonals(int numvals, const int* rows, const double* vals,
		int& status);

	/* write non-zero values to stream as {row,col,value} */
	void PrintNonZero(ostream& out) const;

protected:

	/* return correct dimension */
	int InitGuessLength(void) const;
	int RHSLength(void) const;

	/* check parameters and solver with given rhs and initguess */
	void SolveDriver(double* rhs, double* initguess);

private:

	/* configure the update, bindx, and values array and return
	 * their length. if the columns indices for each row in bindx
	 * are sorted, is_sorted returns 1 and 0 otherwise */
	virtual int SetMSRData(int** update, int** bindx, double** val,
		int& is_sorted) = 0;
	
	/* allocate memory for quick find */
	void SetUpQuickFind(void);

	/* free memory allocated by Aztec.lib */
	void FreeAztecMemory(void);	
	
	/* free memory which is only used in multi-processor execution */
	void FreeAztec_MP_Memory(void);
	
	/* sort active columns in bindx into ascending order */
	void Sort_bindx(void);

protected:

	/** output stream for messages */
	ostream& fMessage;

	/** multiprocessor communicator */
	CommunicatorT& fCommunicator;

	/* number of update rows */
	int Start_update; //1,...
	int N_update;

	/* Aztec work arrays */

	/* See Aztec User's Guide for more information   */
	/* on the variables that follow.                 */

	/* parameter arrays */
	int*    proc_config; // Processor information:
	int*    options;     // Array used to select solver options.
	double* params;      // User selected solver paramters.
	double* status;      // Information returned from AZ_solve()
// indicating success or failure.

	/* dynamically allocated arrays */
	int* data_org;       // Array to specify data layout
	int* external;       // vector elements needed by this node.
	int* update_index;   // ordering of update[] and external[]
	int* extern_index;   // locally on this processor. For example
// update_index[i] gives the index
// location of the vector element which
// has the global index 'update[i]'.

	/* quick find data */
	int     QF_shift;
	int*    update_bin;
	int*    srow_dex; // space for sorted row indices
	double* srow_val; // space for sorted row values

	/* DVBR index arrays - not used */
	int* rpntr;
	int* cpntr;
	int* indx;
	int* bpntr;
	
private:	

	/* used but not managed */
	int*    update; // global indices updated on this processor
	int*    bindx;  // MSR structure data
	double* val;    // matrix value array
		// derived classes should manage the memory and initialization
		// of these arrays. initialization occurs in SetMSRData()

	int* bindx_transform; // version set by AZ_transform
};

} // namespace Tahoe 
#endif /*__AZTEC__ */
#endif /* _AZTEC_BASE_T_H_ */
