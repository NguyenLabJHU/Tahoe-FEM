/* $Id: FEExecutionManagerT.h,v 1.24 2004-04-09 02:03:11 hspark Exp $ */
/* created: paklein (09/21/1997) */
#ifndef _FE_EXECMAN_T_H_
#define _FE_EXECMAN_T_H_

/* element configuration header */
#include "ElementsConfig.h"

/* base class */
#include "ExecutionManagerT.h"

/* direct members */
#include "IOBaseT.h"
#include "ofstreamT.h"

namespace Tahoe {

/* forward declarations */
template <class TYPE> class ArrayT;
class iArrayT;
class OutputSetT;
class IOManager;
class FEManagerT;
class PartitionT;
class ModelManagerT;
class FEManagerT_bridging;
class dArray2DT;
#ifdef __DEVELOPMENT__
class FEManagerT_THK;
#endif
class dArray2DT;
class StringT;
class EAMFCC3D;

/** class to handle file driven finite element simulations */
class FEExecutionManagerT: public ExecutionManagerT
{
public:

	/** constructor */
	FEExecutionManagerT(int argc, char* argv[], char job_char, char batch_char,
		CommunicatorT& comm);

protected:

	/** add the command line option to the list. \returns true if the option was
	 * added, false otherwise. */
	virtual bool AddCommandLineOption(const char* str);

	/** remove the command line option to the list. \returns true if the option was
	 * removed, false otherwise. */
	virtual bool RemoveCommandLineOption(const char* str);

	/** overloaded */
	virtual void RunJob(ifstreamT& in, ostream& status);

	/** \name basic MP support */
	/*@{*/
	int Rank(void) const;
	int Size(void) const;
	/*@}*/

private:

	/** \name execution modes */
	/*@{*/
	/** enum for execution modes */
	enum ModeT {
        kJob = 0,
  kDecompose = 1,
       kJoin = 2,
   kBridging = 3,
        kTHK = 4,
        kDTD = 5
	};
	
	/** standard serial driver */
	void RunJob_serial(ifstreamT& in, ostream& status) const;
	
	/** parallel driver */
	void RunJob_parallel(ifstreamT& in, ostream& status) const;

	/** generate decomposition files */
	void RunDecomp_serial(ifstreamT& in, ostream& status, CommunicatorT& comm, int size = -1) const;

	/** join parallel results files */
	void RunJoin_serial(ifstreamT& in, ostream& status, CommunicatorT& comm, int size = -1) const;

	/** multi-Tahoe, bridging scale test */
	void RunBridging(ifstreamT& in, ostream& status) const;

#ifdef __DEVELOPMENT__
	/** time history kernel tests */
	void RunTHK(ifstreamT& in, ostream& status) const;
#endif

	/** dump current parameter description file */
	void RunWriteDescription(int doc_type) const;
	/*@}*/

#ifdef BRIDGING_ELEMENT
	/** \name bridging scale with different integrators */
	/*@{*/
	/** quasistatic multi-Tahoe bridging scale using a staggered solution strategy */
	void RunStaticBridging_staggered(FEManagerT_bridging& continuum, FEManagerT_bridging& atoms,
		ofstream& log_out) const;

	/** quasistatic multi-Tahoe bridging scale using a monolithic solution strategy */
	void RunStaticBridging_monolithic(ifstreamT& in, FEManagerT_bridging& continuum, FEManagerT_bridging& atoms,
		ofstream& log_out) const;
        
#ifdef __DEVELOPMENT__
	/** dynamic multi-Tahoe bridging scale */
	void RunDynamicBridging(FEManagerT_bridging& continuum, FEManagerT_THK& atoms,
		ofstream& log_out, ifstreamT& in) const;
#endif
				
	/** calculate MD internal force as a function of total displacement u */
	const dArray2DT& InternalForce(dArray2DT& totalu, FEManagerT_bridging& atoms) const;
	/*@}*/
#endif

	/** print message on exception */
	void Rewind(ifstreamT& in, ostream& status) const;

	/** extract the model file name from the stream */
	void GetModelFile(ifstreamT& in, StringT& model_file,
		IOBaseT::FileTypeT& format) const;

	/** \name generate decomposition data */
	/*@{*/
	/** name calls one of decomposition methods below based on user input */
	void Decompose(ifstreamT& in, int size, int decomp_type, CommunicatorT& comm,
		const StringT& model_file, IOBaseT::FileTypeT format) const;

	/** graph-based decomposition. Partition model based on the connectivites
	 * in the model files and those generated at run time. The actual
	 * decomposition is calculated by a FEManagerT_mpi. */
	void Decompose_graph(ifstreamT& in, int size, CommunicatorT& comm, 
		const StringT& model_file, IOBaseT::FileTypeT format) const;

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
	void Decompose_atom(ifstreamT& in, int size, const StringT& model_file,
		IOBaseT::FileTypeT format) const;

	/** spatial decomposition. Partition model based on a grid. */
	void Decompose_spatial(ifstreamT& in, int size, const StringT& model_file,
		IOBaseT::FileTypeT format) const;
	/*@}*/

	/** returns true if a new decomposition is needed */
	bool NeedDecomposition(const StringT& model_file, int size) const;

	/** returns true if the global output model file is not found */
	bool NeedModelFile(const StringT& model_file, IOBaseT::FileTypeT format) const;

	/** construct and return the local IOManager */
	IOManager* NewLocalIOManager(const FEManagerT* global_FEman,
		const iArrayT& output_map) const;
		
	/** \name write partial geometry files 
	 * \param partition partition information for the part to be written
	 * \param model_ALL ModelManagerT accessing the total model database
	 * \param partial_file path to the partial geometry file
	 */
	/*@{*/
	/** write partial file for the given format */
	void EchoPartialGeometry(const PartitionT& partition, ModelManagerT& model_ALL,
		const StringT& partial_file, IOBaseT::FileTypeT format) const;

	/** write partial model file in ExodusII format */
	void EchoPartialGeometry_ExodusII(const PartitionT& partition,
		ModelManagerT& model_ALL, const StringT& partial_file) const;

	/** write partial model file in TahoeII format */
	void EchoPartialGeometry_TahoeII(const PartitionT& partition,
		ModelManagerT& model_ALL, const StringT& partial_file) const;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _FE_EXECMAN_T_H_ */
