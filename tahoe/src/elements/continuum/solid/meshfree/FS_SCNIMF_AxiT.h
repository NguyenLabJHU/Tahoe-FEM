/* $Id: FS_SCNIMF_AxiT.h,v 1.11 2005-03-01 08:26:29 paklein Exp $ */
#ifndef _FS_SCNIMF_AXI_T_H_
#define _FS_SCNIMF_AXI_T_H_

/* base class */
#include "FS_SCNIMFT.h"

/* direct members */
#include "FSMatSupportT.h"

namespace Tahoe {


/** base class for particle types */
class FS_SCNIMF_AxiT: public FS_SCNIMFT
{
public:

	/** constructor */
	FS_SCNIMF_AxiT(const ElementSupportT& support);

	/** write output. ParticleT::WriteOutput only writes search grid statistics.
	 * Sub-classes are responsible for writing data for each particle, given the
	 * variables names returned by ParticleT::GenerateOutputLabels. */
	virtual void WriteOutput(void);

	/** trigger reconfiguration */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

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

	/** return number of values for each output variable */
	virtual void SetOutputCount(const iArrayT& flags, iArrayT& counts) const;
	
	virtual void CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const;
			
	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels);

	/** assemble particle mass matrix into LHS of global equation system */
	virtual void AssembleParticleMass(const double rho);
	
};

} /* namespace Tahoe */

#endif /* _FD_SCNIMF_AXI_T_H_ */


