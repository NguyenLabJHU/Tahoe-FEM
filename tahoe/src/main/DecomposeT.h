/* $Id: DecomposeT.h,v 1.2 2004-09-28 15:35:37 paklein Exp $ */
#ifndef _FE_DECOMPOSE_T_H_
#define _FE_DECOMPOSE_T_H_

/* direct members */
#include "IOBaseT.h"

namespace Tahoe
{

/* forward declarations */
class StringT;
class PartitionT;
class ModelManagerT;
class CommunicatorT;
class GraphT;
class FEManagerT;

/** class to handle decomposition of models for parallel execution */
class DecomposeT
{
public:

	/** constructor */
	DecomposeT(void);

	/** decomposes model if */
	void CheckDecompose(const StringT& input_file, int size, int decomp_type, CommunicatorT& comm,
		const StringT& model_file, IOBaseT::FileTypeT format, const ArrayT<StringT>& commandlineoptions) const;
	
	/** returns true if a new decomposition is needed */
	bool NeedDecomposition(const StringT& model_file, int size) const;
	
	/** returns true if the global output model file is not found */
	bool NeedModelFile(const StringT& model_file, IOBaseT::FileTypeT format) const;
	
	/** write partial file for the given format */
	void EchoPartialGeometry(const PartitionT& partition, ModelManagerT& model_ALL,
		const StringT& partial_file, IOBaseT::FileTypeT format) const;
private:

	/** \name decomposition methods */ 
	/*@{*/
	/** graph-based decomposition. Partition model based on the connectivites
 	 * in the model files and those generated at run time. The actual
  	 * decomposition is calculated by a FEManagerT. */
	void Decompose_graph(const StringT& input_file, int size, CommunicatorT& comm, 
		const StringT& model_file, IOBaseT::FileTypeT format, const ArrayT<StringT>& commandlineoptions) const;

	/** "atom" decomposition. Partition model by dividing global list
 	 * of coordinates into sequential, nearly equal length lists. The
 	 * number of atoms per partition is \f$ \frac{N}{n_p} \f$ for
 	 * all partitions except the last, which also includes any remainder.
 	 * \f$ N \f$ is the total number nodes and \f$ n_p \f$ is the number
 	 * of partitions. The partition for a given node is then given by
 	 \f[
 		 p_i = floor \left( \frac{i n_p}{N} \right).
 	 \f]
 	 */
	void Decompose_atom(const StringT& input_file, int size, const StringT& model_file,
		IOBaseT::FileTypeT format, const ArrayT<StringT>& commandlineoptions) const;

	/** spatial decomposition. Partition model based on a grid. */
	void Decompose_spatial(const StringT& input_file, int size, const StringT& model_file,
		IOBaseT::FileTypeT format) const;
	/*@}*/
	
	/** \name graph decompositon methods */
	/*@{*/
	/** domain decomposition (graph is returned) */
	void DoDecompose_graph(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int method, const FEManagerT& feman) const;
	void DoDecompose_graph_1(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int method, const FEManagerT& feman) const;
	void DoDecompose_graph_2(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int method, const FEManagerT& feman) const;		
	/*@}*/

	/** \name write partial geometry files 
 	* \param partition partition information for the part to be written
 	* \param model_ALL ModelManagerT accessing the total model database
 	* \param partial_file path to the partial geometry file
 	*/
	/*@{*/
	/** write partial model file in ExodusII format */
	void EchoPartialGeometry_ExodusII(const PartitionT& partition,
		ModelManagerT& model_ALL, const StringT& partial_file) const;

	/** write partial model file in TahoeII format */
	void EchoPartialGeometry_TahoeII(const PartitionT& partition,
		ModelManagerT& model_ALL, const StringT& partial_file) const;
	/*@}*/
	
};
}
#endif 