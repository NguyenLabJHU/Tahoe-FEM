/* $Id: FEDecomposeT.h,v 1.1.2.2 2004-09-15 02:14:15 d-farrell2 Exp $ */
/* created: d-farrell2 (08/03/2004) */
#ifndef _FE_DECOMPOSE_T_H_
#define _FE_DECOMPOSE_T_H_


// Includes
/* program parameters */
#include "GlobalT.h"

/* base class */
#include "iConsoleObjectT.h"
#include "ParameterInterfaceT.h"

/* direct members */
#include "IOBaseT.h"
#include "iArrayT.h"
#include "StringT.h"
#include "ElementListT.h"
#include "IOBaseT.h"
#include "iArray2DT.h"
/* direct members, formerly in FEManagerT_mpi.h, DEF 28 July 04 */
#include "PartitionT.h"
#include "dArray2DT.h"
#include "ofstreamT.h"
#include "SolverT.h"
#include "IOManager.h"

#include "ios_fwd_decl.h"

namespace Tahoe
{
	/* forward declarations */
	template <class TYPE> class ArrayT;
	class iArrayT;
	class OutputSetT;
	class IOManager;
	class FEManagerT;
	class PartitionT;
	class ModelManagerT;
	class dArray2DT;
	class dArray2DT;
	class StringT;
	class ParameterListT;
	/* forward declarations, formerly in FEManagerT_mpi.h, DEF 28 July 04 */
	class IOManager_mpi;
	class PartitionT;
	
	class FEDecomposeT
	{
	public:
		// Constructor(?) DEF 3 Aug 04
		FEDecomposeT();
		
		// Destructor
		~FEDecomposeT();
			
		// calls one of decomposition methods below based on user input
		void CheckDecompose(const StringT& input_file, int size, int decomp_type, CommunicatorT& comm,
		const StringT& model_file, IOBaseT::FileTypeT format, const ArrayT<StringT>& commandlineoptions) const; // was FEExecutionManagerT::Decompose(..)
		
		/* domain decomposition (graph is returned), formerly in FEManagerT_mpi.h, DEF 28 July 04 */
		void Decompose(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int method, const FEManagerT& feman) const;
		
		/** returns true if a new decomposition is needed */
		bool NeedDecomposition(const StringT& model_file, int size) const;
		
		/** returns true if the global output model file is not found */
		bool NeedModelFile(const StringT& model_file, IOBaseT::FileTypeT format) const;
		
		/** write partial file for the given format */
		void EchoPartialGeometry(const PartitionT& partition, ModelManagerT& model_ALL,
			const StringT& partial_file, IOBaseT::FileTypeT format) const;
	private:
		// These formerly in FEExecutionManagerT.h DEF 3 Aug 04
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
		
		/** \name decomposition methods, were in FEManagerT.h, 3 Aug 04 */ 
		/*@{*/
		void DoDecompose_1(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int method, const FEManagerT& feman) const;
		void DoDecompose_2(ArrayT<PartitionT>& partition, GraphT& graph, bool verbose, int method, const FEManagerT& feman) const;		
		/*@}*/
		
		/** returns true if the global output model file is not found */
		//bool NeedModelFile(const StringT& model_file, IOBaseT::FileTypeT format) const;
		
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