/* $Id: ParticlePairT.h,v 1.3 2002-11-22 01:49:45 paklein Exp $ */
#ifndef _PARTICLE_PAIR_T_H_
#define _PARTICLE_PAIR_T_H_

/* base class */
#include "ParticleT.h"

/* direct members */
#include "RaggedArray2DT.h"

namespace Tahoe {

/** base class for particle types */
class ParticlePairT: public ParticleT
{
public:

	/** constructor */
	ParticlePairT(const ElementSupportT& support, const FieldT& field);

	/** initialization */
	virtual void Initialize(void);

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

	/** trigger reconfiguration */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/** close current time increment. Since this is called only once per
	 * time step, this is used to increment counters. */
	virtual void CloseStep(void);

protected:

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the LHS matrix */
	virtual void LHSDriver(void);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	/*@}*/
	
	/** set neighborlists. Recalculates the neighborlists at intervals defined
	 * by ParticlePairT::fReNeighborIncr.
	 * \param force if true, forces recalculation of neighbors regardless of
	 *        of the state of the counters
	 * \return true if the configuration has changed */
	bool SetConfiguration(bool force = false);
	
private:

	/** neighbor cut-off distance */
	double fNeighborDistance;

	/** number of steps between reseting neighbor lists */
	int fReNeighborIncr;

	/** neighbor lists */
	RaggedArray2DT<int> fNeighbors;

	/** equation numbers */
	RaggedArray2DT<int> fEqnos;

	/** \name run time information. Incremented in */
	/*@{*/
	int fReNeighborCounter;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _PARTICLE_PAIR_T_H_ */
