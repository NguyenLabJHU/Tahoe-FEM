/* $Id: FS_SCNIMF_AxiT.h,v 1.4 2004-08-04 22:00:23 cjkimme Exp $ */
#ifndef _FS_SCNIMF_AXI_T_H_
#define _FS_SCNIMF_AXI_T_H_

/* base class */
#include "SCNIMFT.h"

/* direct members */
#include "FSMatSupportT.h"

namespace Tahoe {


/** base class for particle types */
class FS_SCNIMF_AxiT: public SCNIMFT
{
public:

	/** constructor */
	FS_SCNIMF_AxiT(const ElementSupportT& support, const FieldT& field);
	FS_SCNIMF_AxiT(const ElementSupportT& support);

	/** destructor */
	~FS_SCNIMF_AxiT(void);

	/** write output. ParticleT::WriteOutput only writes search grid statistics.
	 * Sub-classes are responsible for writing data for each particle, given the
	 * variables names returned by ParticleT::GenerateOutputLabels. */
	virtual void WriteOutput(void);

	/* compute specified output parameter and send for smoothing */
	virtual void SendOutput(int kincode);

	/** trigger reconfiguration */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/** \name restart functions */
	/*@{*/
	/** write restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::ReadRestart implementation. */
	virtual void WriteRestart(ostream& out) const;

	/** read restart data to the output stream. Should be paired with
	 * the corresponding ElementBaseT::WriteRestart implementation. */
	virtual void ReadRestart(istream& in);
	/*@}*/

	/** */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
	
	/** */
	virtual void RHSDriver(void);

protected: /* for derived classes only */
	
	virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;
	
	virtual MaterialListT* NewMaterialList(const StringT& name, int size);
	
	/** translate internal storage of bVector to Strain-Displacement matrix */	
	void bVectorToMatrix(double *bVector, dMatrixT& BJ);
	
	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels);

	/** compute B matrices for strain smoothing/nodal integration */
	virtual void ComputeBMatrices(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	/*@}*/

protected:

	FSMatSupportT* fFSMatSupport;
	int fNeedsOffset;
	
private:
	ArrayT< LinkedListT<double> > circumferentialWorkSpace;
	RaggedArray2DT<double> circumferential_B;
	
};

} /* namespace Tahoe */

#endif /* _FD_SCNIMF_AXI_T_H_ */


