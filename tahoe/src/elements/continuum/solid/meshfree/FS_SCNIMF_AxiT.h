/* $Id: FS_SCNIMF_AxiT.h,v 1.7 2004-10-28 20:30:53 gjwagne Exp $ */
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

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/	
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected: /* for derived classes only */
	
	virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;
	
	virtual MaterialListT* NewMaterialList(const StringT& name, int size);
	
	/** translate internal storage of bVector to Strain-Displacement matrix */	
	void bVectorToMatrix(double *bVector, dMatrixT& BJ);
	
	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels);

	/** assemble particle mass matrix into LHS of global equation system */
	virtual void AssembleParticleMass(const double rho);

	/** compute B matrices for strain smoothing/nodal integration */
	virtual void ComputeBMatrices(void);

protected:

	FSMatSupportT* fFSMatSupport;
	int fNeedsOffset;
	
private:
	ArrayT< LinkedListT<double> > circumferentialWorkSpace;
	RaggedArray2DT<double> circumferential_B;

	/* deformation gradients passed to the materials */
	ArrayT<dMatrixT> fF_list;
	ArrayT<dMatrixT> fF_last_list;
};

} /* namespace Tahoe */

#endif /* _FD_SCNIMF_AXI_T_H_ */


