/* $Id: SS_SCNIMFT.h,v 1.3 2004-07-29 23:42:06 cjkimme Exp $ */
#ifndef _SS_SCNIMF_T_H_
#define _SS_SCNIMF_T_H_

/* base class */
#include "SCNIMFT.h"

/* direct members */
#include "SSMatSupportT.h"

namespace Tahoe {

/** base class for particle types */
class SS_SCNIMFT: public SCNIMFT
{
public:

	/** constructor */
	SS_SCNIMFT(const ElementSupportT& support, const FieldT& field);
	SS_SCNIMFT(const ElementSupportT& support);

	/** destructor */
	~SS_SCNIMFT(void);

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
	
	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels);

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
	
	//virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected: /* for derived classes only */
	
	virtual void ReadMaterialData(void);
	virtual MaterialListT* NewMaterialList(const StringT&, int size);
	
	/** translate internal storage of bVector to Strain-Displacement matrix */	
	void bVectorToMatrix(double *bVector, dMatrixT& BJ);
	

	SSMatSupportT* fSSMatSupport;
	
	/** offset to material needs */
	int fNeedsOffset; 
  

};

} /* namespace Tahoe */

#endif /* _SS_SCNIMF_T_H_ */


