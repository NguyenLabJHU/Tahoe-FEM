/* $Id: JoinOutputT.h,v 1.1 2002-01-09 12:36:51 paklein Exp $ */
/* created: paklein (03/24/2000)                                          */

#ifndef _JOIN_OUTPUT_T_H_
#define _JOIN_OUTPUT_T_H_

/* direct members */
#include "IOBaseT.h"
#include "PartitionT.h"
#include "MapSetT.h"
#include "dArray2DT.h"
#include "StringT.h"

/* forward declarations */
class IOManager;

class JoinOutputT
{
public:

	/* constructor */
	JoinOutputT(ifstreamT& in, const StringT& model_file,
		const StringT& global_model_file, IOBaseT::FileTypeT file_type,
		int size);

	/* destructor */
	~JoinOutputT(void);

	/* do join */
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
	void ResultFileName(int part, int group, StringT& name) const;

	/* returns the number of output steps */
	int NumOutputSteps(int group) const;
	
	/* retrieve output labels */
	void OutputLabels(int group, ArrayT<StringT>& node_labels,
		ArrayT<StringT>& element_labels) const;

private:

	/* file name info */
	const StringT fJobFile;
	const StringT fModelFile;
	const StringT fGlobalModelFile;
	
	/* data format */
	IOBaseT::FileTypeT fFileType;
	
	/* partition data */
	ArrayT<PartitionT> fPartitions;
	
	/* I/O manager */
	IOManager* fIO;
	
	/* model data */
	dArray2DT fCoordinates;
	ArrayT<iArray2DT> fConnects;
	
	/* maps (for each output set) from processor to global position */
	ArrayT<MapSetT> fMapSets;	
};

#endif /* _JOIN_OUTPUT_T_H_ */
