/* $Id: CommManagerT.h,v 1.1 2002-12-05 08:31:13 paklein Exp $ */
#ifndef _COMM_MANAGER_T_H_
#define _COMM_MANAGER_T_H_

/* direct members */
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class PartitionT;
class CommunicatorT;
class ModelManagerT;
class iArrayT;

/** manage processor to processor transactions. Manages partition information.
 * Creates ghosts nodes. Manages communication lists. Manipulates the
 * ModelManagerT to create the environment for the local processor. */
class CommManagerT
{
public:

	/** constructor */
	CommManagerT(CommunicatorT& comm, int nsd);

	/** set partition information */
	void SetPartition(PartitionT* partition) { fPartition = partition; };

	/** set the model manager */
	void SetModelManager(ModelManagerT* model_manager) { fModelManager = model_manager; };

	/** \name setting periodic boundaries */
	/*@{*/
	void SetPeriodicBoundaries(int i, double x_i_min, double x_i_max);
	void ClearPeriodicBoundaries(int i);
	/*@}*/

	/** configure the current local coordinate list and register it with the
	 * model manager */
	void Configure(double range);

private:

	/** determine the local coordinate bounds 
	 * \param coords coordinate list
	 * \param local rows of coords to use in determining bounds
	 * \param bounds returns with the bounds of the given points */
	void GetBounds(const dArray2DT& coords, const iArrayT& local, dArray2DT& bounds) const;

private:

	/** communicator */
	CommunicatorT& fComm;

	/** \name periodic boundaries */
	/*@{*/
	/** flags to indicate if periodic boundary conditions are imposed */
	ArrayT<bool> fIsPeriodic;
	
	/** rows give the lower and upper periodic bounds for that coordinate */
	dArray2DT fPeriodicBoundaries;
	/*@}*/
	
	/** processor bounds */
	dArray2DT fBounds;

	/** partition information */
	PartitionT* fPartition;

	/** the model manager */
	ModelManagerT* fModelManager;
};

} /* namespace Tahoe */

#endif /* _COMM_MANAGER_T_H_ */
