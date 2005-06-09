/* $Id: JoinOutputT.h,v 1.6 2005-06-09 00:03:09 paklein Exp $ */
/* created: paklein (03/24/2000) */

#ifndef _JOIN_OUTPUT_T_H_
#define _JOIN_OUTPUT_T_H_

/* direct members */
#include "IOBaseT.h"
#include "PartitionT.h"
#include "MapSetT.h"
#include "dArray2DT.h"
#include "StringT.h"

namespace Tahoe {

/* forward declarations */
class OutputBaseT;
class ModelManagerT;

/** class to join results from a parallel calculation */
class JoinOutputT
{
public:

	/** constructor */
	JoinOutputT(const StringT& param_file, const StringT& model_file,
		IOBaseT::FileTypeT model_file_type, IOBaseT::FileTypeT results_file_type, 
		OutputBaseT* output, int size);

	/** destructor */
	~JoinOutputT(void);

	/** do join */
	void Join(void);

private:

	/* set output */
	void SetOutput(void);

	/* set assembly maps */
	void SetMaps(void);

	/* resident partition for each node */
	void SetNodePartitionMap(iArrayT& node_partition);

	/* return the global node numbers of the set nodes residing
	 * in the partition */
	void PartitionSetNodes(int partition, const iArrayT& node_part_map,
		const iArrayT& set_nodes, iArrayT& nodes) const;

	/* determine map from local nodes into global array, such that:
	 *
	 *             global[lg_map[i]] = local[i]
	 */
	void SetInverseMap(const iArrayT& global, iArrayT& inv_global,
		int& shift, int fill) const;
	void SetAssemblyMap(const iArrayT& inv_global, int shift,
		const iArrayT& local, iArrayT& lg_map) const;		

	/* check that assembly maps are compact and complete */
	void CheckAssemblyMaps(void);

	/* generate output file name */
	void ResultFileName(int part, int group, bool changing, int print_step, StringT& name) const;

	/* returns the number of output steps */
	int NumOutputSteps(int group) const;
	
	/* retrieve output labels */
	void OutputLabels(int group, bool changing, ArrayT<StringT>& node_labels,
		ArrayT<StringT>& element_labels) const;

private:

	/** parameters file name */
	const StringT fJobFile;

	/** model database manager */
	ModelManagerT* fModel;

	/** geometry database file type */
	IOBaseT::FileTypeT fResultsFileType;
	
	/** partition data */
	ArrayT<PartitionT> fPartitions;
	
	/** output formatter */
	OutputBaseT* fOutput;

	/** maps (for each output set) from processor to global position */
	ArrayT<MapSetT> fMapSets;	
};

} // namespace Tahoe 
#endif /* _JOIN_OUTPUT_T_H_ */
