/* $Id: ParticlePairT.h,v 1.17 2004-06-26 06:06:26 paklein Exp $ */
#ifndef _PARTICLE_PAIR_T_H_
#define _PARTICLE_PAIR_T_H_

/* base class */
#include "ParticleT.h"

/* direct members */
#include "RaggedArray2DT.h"
#include "VariArrayT.h"
#include "ofstreamT.h"

namespace Tahoe {

/* forward declarations */
class PairPropertyT;
class StringT;

/** base class for particle types */
class ParticlePairT: public ParticleT
{
public:

	/** constructor */
	ParticlePairT(const ElementSupportT& support, const FieldT& field);
	ParticlePairT(const ElementSupportT& support);

	/** collecting element group equation numbers */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** class initialization. Among other things, element work space
	 * is allocated and connectivities are read. */
	virtual void Initialize(void);

	/** \name connectivities.
	 * See ElementBaseT::ConnectsX and ElementBaseT::ConnectsU for more
	 * information about what these are used for */
	/*@{*/
	/** collecting element geometry connectivities */
	virtual void ConnectsX(AutoArrayT<const iArray2DT*>& connects) const;

	/** collecting element field connectivities */
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	             AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;
	/*@}*/

	/** write output. Writes output data for all tags listed in
	 * ParticleT::fGlobalTag. The values per node are those specified
	 * by ParticlePairT::GenerateOutputLabels. */
	virtual void WriteOutput(void);

	/** compute the part of the stiffness matrix associated with rows of the given 
	 * particles. This is the mixed part of the stiffness matrix involving free
	 * particles and ghost particles which have prescribed motion. */
	virtual void FormStiffness(const InverseMapT& col_to_col_eq_row_map,
		const iArray2DT& col_eq, dSPMatrixT& stiffness);

	/** access to the neighbor pair list */
	const RaggedArray2DT<int>& Neighbors(void) const { return fNeighbors; };

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

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the LHS matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	void RHSDriver1D(void);
	void RHSDriver2D(void);
	void RHSDriver3D(void);
	/*@}*/
	
	/** set neighborlists and any other system configuration information
	 * based on the current information. Uses ParticleT::GenerateNeighborList
	 * to determine the neighborlists. */
	virtual void SetConfiguration(void);

	/** construct the list of properties from the given input stream */
	virtual void EchoProperties(ifstreamT& in, ofstreamT& out);

	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels) const;


	/*nearest neighbor list*/
	RaggedArray2DT<int> NearestNeighbors;
	RaggedArray2DT<int> RefNearestNeighbors;

	/** return a new pair property or NULL if the name is invalid */
	PairPropertyT* New_PairProperty(const StringT& name, bool throw_on_fail) const;

private:

	/** particle pair-properties list */
	ArrayT<PairPropertyT*> fPairProperties;

	/** neighbor lists */
	RaggedArray2DT<int> fNeighbors;


	/** equation numbers */
	RaggedArray2DT<int> fEqnos;

	/** \name workspace for ParticlePairT::RHSDriver. Used to accumulate the force for
	 * a single row of ParticlePairT::fNeighbors. */
	/*@{*/
	dArrayT fForce_list;
	VariArrayT<double> fForce_list_man;

	/** constant matrix needed to compute the stiffness */
	dMatrixT fOneOne;
	
	/** new variables for file I/O */
	bool fopen;
	ofstreamT fout, fout2;
	StringT fsummary_file, fsummary_file2;
	/*@}*/

	
	
};

} /* namespace Tahoe */

#endif /* _PARTICLE_PAIR_T_H_ */
