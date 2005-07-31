/* $Id: IOManager_mpi.h,v 1.1.1.1 2001-01-29 08:20:21 paklein Exp $ */
/* created: paklein (03/14/2000)                                          */

#ifndef _IOMANAGER_MPI_H_
#define _IOMANAGER_MPI_H_

/* base class */
#include "IOManager.h"

/* direct members */
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "MapSetT.h"
#ifdef __MPI__
#include "mpi.h"
#endif

/* forward declarations */
class PartitionT;

class IOManager_mpi: public IOManager
{
public:

	/* constructor */
	IOManager_mpi(ifstreamT& in, const iArrayT& io_map,
		const IOManager& local_IO, const PartitionT& partition,
		const StringT& global_geom_file, IOBaseT::FileTypeT format);

	/* destructor */
	virtual ~IOManager_mpi(void);

#ifdef __MPI__
	/* distribute/assemble/write output */
	virtual void WriteOutput(double time, int ID, const dArray2DT& n_values,
		const dArray2DT& e_values);
#endif

private:

	/* communicate output counts */
	void SetCommunication(const IOManager& local_IO);

	/* return the global node numbers of the set nodes residing
	 * on the current partition */
	void GlobalSetNodes(const iArrayT& local_set_nodes, iArrayT& nodes);

	/* determine map from local nodes into global array, such that:
*
	 *             global[lg_map[i]] = local[i]
	 */
	void SetInverseMap(const iArrayT& global, iArrayT& inv_global,
		int& shift, int fill) const;
	void SetAssemblyMap(const iArrayT& inv_global, int shift,
		const iArrayT& local, iArrayT& lg_map) const;		

	/* MPI information */
	int Rank(void) const;
	int Size(void) const;

#ifdef __MPI__
	/* clear all outstanding requests - returns 1 of all OK */
	int Clear(ArrayT<MPI_Request>& requests);
#endif

	/* check that assembly maps are compact and complete */
	void CheckAssemblyMaps(void);

	/* load global geometry */
	void ReadOutputGeometry(const StringT& global_geom_file,
		IOBaseT::FileTypeT format);

private:

	/* ID to processor map */
	const iArrayT fIO_map;

	/* partition info */
	const PartitionT& fPartition;

	/* global model data */
	dArray2DT fCoordinates; //NOTE: right now storing ALL of them
	iArrayT   fNodeMap;
	ArrayT<iArray2DT> fConnectivities; //NOTE: only reading

	/* communication maps */
	iArray2DT fNodeCounts;    // [element sets] x [number of processors]
	iArray2DT fElementCounts; // [element sets] x [number of processors]
	
	/* maps (for each output set) from processor to global position */
	ArrayT<MapSetT> fMapSets;

	/* number of outgoing per set */
	iArrayT fOutNodeCounts; // assumes values for border nodes will be at the end
	  	
#ifdef __MPI__
	/* requests for non-blocking communication */
	ArrayT<MPI_Request> fRecvRequest;
	ArrayT<MPI_Request> fSendRequest;
#endif
};

#endif /* _IOMANAGER_MPI_H_ */
