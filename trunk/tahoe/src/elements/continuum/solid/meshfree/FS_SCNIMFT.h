/* $Id: FS_SCNIMFT.h,v 1.8 2005-03-06 04:02:35 cjkimme Exp $ */
#ifndef _FS_SCNIMF_T_H_
#define _FS_SCNIMF_T_H_

/* base class */
#include "SCNIMFT.h"

/* direct members */
#include "FSMatSupportT.h"

namespace Tahoe {

/** base class for particle types */
class FS_SCNIMFT: public SCNIMFT
{
public:

	/** constructor */
	FS_SCNIMFT(const ElementSupportT& support, const FieldT& field);
	FS_SCNIMFT(const ElementSupportT& support);

	/** destructor */
	~FS_SCNIMFT(void);

	/** write output. ParticleT::WriteOutput only writes search grid statistics.
	 * Sub-classes are responsible for writing data for each particle, given the
	 * variables names returned by ParticleT::GenerateOutputLabels. */
	virtual void WriteOutput(void);

	/** trigger reconfiguration */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/** compute stiffness matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
	
	/** compute residual force */
	virtual void RHSDriver(void);

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
	
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected: /* for derived classes only */

	virtual dMatrixT& TransformModuli(const dMatrixT& moduli, const dMatrixT& F, dMatrixT& Csig);
	
	virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;
	virtual MaterialListT* NewMaterialList(const StringT& name, int size);
	
	/** translate internal storage of bVector to Strain-Displacement matrix */	
	void bVectorToMatrix(double *bVector, dMatrixT& BJ);

	FSMatSupportT* fFSMatSupport;

	/** \name deformation gradients passed to the constitutive models */
	/*@{*/
	ArrayT<dMatrixT> fF_list;
	ArrayT<dMatrixT> fF_last_list;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _FD_SCNIMF_T_H_ */


