/* $Id: ParticlePairT.h,v 1.7 2002-11-30 16:37:25 paklein Exp $ */
#ifndef _PARTICLE_PAIR_T_H_
#define _PARTICLE_PAIR_T_H_

/* base class */
#include "ParticleT.h"

/* direct members */
#include "RaggedArray2DT.h"
#include "VariArrayT.h"

namespace Tahoe {

/* forward declarations */
class PairPropertyT;

/** base class for particle types */
class ParticlePairT: public ParticleT
{
public:

	/** constructor */
	ParticlePairT(const ElementSupportT& support, const FieldT& field);

	/** destructor */
	~ParticlePairT(void);

	/** collecting element group equation numbers */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

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

protected:

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the LHS matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	void RHSDriver_1(void);
	void RHSDriver_2(void);
	void RHSDriver_3(void);
	void RHSDriver3D_3(void);
	/*@}*/
	
	/** set neighborlists and any other system configuration information
	 * based on the current information. Uses ParticleT::GenerateNeighborList
	 * to determine the neighborlists. */
	virtual void SetConfiguration(void);

	/** construct the list of properties from the given input stream */
	virtual void EchoProperties(ifstreamT& in, ofstreamT& out);

	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels) const;

private:

	/** neighbor cut-off distance */
	double fNeighborDistance;

	/** particle properties list */
	ArrayT<PairPropertyT*> fProperties;

	/** neighbor lists */
	RaggedArray2DT<int> fNeighbors;

	/** equation numbers */
	RaggedArray2DT<int> fEqnos;

	/** \name workspace for ParticlePairT::RHSDriver. Used to accumulate the force for
	 * a single row of ParticlePairT::fNeighbors. */
	/*@{*/
	dArrayT fForce_list;
	VariArrayT<double> fForce_list_man;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _PARTICLE_PAIR_T_H_ */
