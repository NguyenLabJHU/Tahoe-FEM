/* $Id */
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
	
protected: /* for derived classes only */
	
	virtual void ReadMaterialData(ifstreamT& in);
	
	virtual void WriteMaterialData(ostream& out) const;
	
	virtual MaterialListT* NewMaterialList(int nsd, int size);
	
	/** translate internal storage of bVector to Strain-Displacement matrix */	
	void bVectorToMatrix(double *bVector, dMatrixT& BJ);
	
	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
		SubListT& sub_sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

protected:

	SSMatSupportT* fSSMatSupport;
	
	/** offset to material needs */
	int fNeedsOffset; 
  

};

} /* namespace Tahoe */

#endif /* _SS_SCNIMF_T_H_ */


