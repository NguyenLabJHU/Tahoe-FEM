/* $Id: IOManager_mpi.h,v 1.16 2003-10-28 07:47:53 paklein Exp $ */
/* created: paklein (03/14/2000) */

#ifndef _IOMANAGER_MPI_H_
#define _IOMANAGER_MPI_H_

/* base class */
#include "IOManager.h"

/* direct members */
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "MapSetT.h"
#include "CommunicatorT.h"

namespace Tahoe {

/* forward declarations */
class PartitionT;
class ModelManagerT;

/** class to support runtime joining of parallel output. The generation of
 * assembly maps differs here from the approach used in JoinOutputT in that
 * in this implementation, only information about the local processor is
 * read and stored. The remaining information is communicated during initialization.
 * The JoinOutputT class reads information about every partition in order to
 * generate the communication and assembly maps. */
class IOManager_mpi: public IOManager
{
public:

	/** constructor 
	 * \param model_file total model database */
	IOManager_mpi(ifstreamT& in, CommunicatorT& comm, const IOManager& local_IO, 
		const PartitionT& partition, const StringT& model_file, IOBaseT::FileTypeT format);

	/** destructor */
	virtual ~IOManager_mpi(void);

	/** distribute/assemble/write output */
	virtual void WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values);

	/** temporarily re-route output to a database with the given filename */
	virtual void DivertOutput(const StringT& outfile);

private:

	/** write the assembly maps. Used for debugging */
	void WriteMaps(ostream& out) const;

	/** communicate output counts */
	void SetCommunication(const IOManager& local_IO);

	/** load balance output across processors 
	 * \param local_IO IOManager for local output 
	 * \param output_map returns with processor handling each output set */
	void SetOutputMap(const IOManager& local_IO, iArrayT& output_map) const;

	/** return the global node numbers of the set nodes residing
	 * on the current partition */
	void GlobalSetNodes(const iArrayT& local_set_nodes, iArrayT& nodes);

	/** determine map from local nodes into global array, such that:
	 *
	 *             global[lg_map[i]] = local[i]
	 */
	void SetInverseMap(const iArrayT& global, iArrayT& inv_global,
		int& shift, int fill) const;
	void SetAssemblyMap(const iArrayT& inv_global, int shift,
		const iArrayT& local, iArrayT& lg_map) const;		

	/** determine the assembly map. Incorporate the partial assemble map
	 * into the global assemble map.
	 * \param set output set index to assemble
	 * \param block_ID block ID for the partial map
	 * \param block_map global block ID as a function of partial index
	 * \param map global output set ID as a function of set index. This array
	 *        is allocate/initialized if passed in empty. Subsequently the
	 *        array length must match the number of elements in the output set */
	void BuildElementAssemblyMap(int set, const StringT& block_ID, const iArrayT& block_map, 
		iArrayT& map) const;

	/** check that assembly maps are compact and complete */
	void CheckAssemblyMaps(void);

	/** read global geometry. Reads \e all of the nodal coordinates, but only
	 * the connectivities for output sets written by the processor.  
	 * \param model_file path to the total model file
	 * \param element_sets processor local output sets
	 * \param format model database format */
	void ReadOutputGeometry(const StringT& model_file,
		const ArrayT<OutputSetT*>& element_sets, IOBaseT::FileTypeT format);

private:

	/** MP communicator */
	CommunicatorT& fComm;

	/** output set index to output processor map 
	 * output_processor[output_set_index] */
	iArrayT fIO_map;

	/** partition info */
	const PartitionT& fPartition;
	
	/** global model geometry */
	ModelManagerT* fOutputGeometry;

	/* communication maps */
	iArray2DT fNodeCounts;    // [element sets] x [number of processors]
	iArray2DT fElementCounts; // [element sets] x [number of processors]
	
	/** maps (for each output set) from processor to global node number */
	ArrayT<MapSetT> fMapSets;

	/** number of outgoing nodes per set. This is the sum of the internal and
	 * border nodes. These ID's are determined assuming the local ID's are
	 * assigned in order as _i, _b, _e. An analogous array for elements is
	 * not needed because all elements are resident, hence the ID's (by
	 * block) can be taken directly from the elements maps stored in the
	 * local PartitionT. Seems to be redundant with IOManager_mpi::fNodeCounts */
	iArrayT fOutNodeCounts;	
};

} // namespace Tahoe 
#endif /* _IOMANAGER_MPI_H_ */
