/* $Id: nIntegratorT.h,v 1.2 2001-08-27 17:12:13 paklein Exp $ */
/* created: paklein (10/14/1996) */

#ifndef _N_CONTROLLERT_H_
#define _N_CONTROLLERT_H_

/* direct members */
#include "KBC_CardT.h"
#include "ArrayT.h"

/* forward declarations */
class dArray2DT;
class iArray2DT;
class dArrayT;
class iArrayT;

/** defines the interface for time integrators of field data. All
 * field arrays must be set with nIntegratorT::RegisterField
 * before calling predictors, correctors, or boundary condition
 * calculations. Note that the predictor is applied to all
 * degrees of freedom, while the correctors are only applied to
 * the \a active degrees of freedom */
class nIntegratorT
{
public:

	/** constructor.
	 * \param order number of time derivatives needed for the field.
	 *        order must be >= 0. The zero-th order field is always
	 *        assumed to exist. */
	nIntegratorT(int order);

	/* destructor */
	virtual ~nIntegratorT(void);

	/** pseudo-boundary conditions for external nodes. For parallel
	 * calculations, external nodes, or "ghost" nodes, are assigned
	 * a boundary condition state in order to remove their associated
	 * degrees of freeom from the local "active" set. */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const = 0;

	/** register field array.
	 * \param field array of field variables
	 * \param order specifies which order time derivative is being
	 *        registered. */
	void RegisterField(dArray2DT& field, int order);

	/** predictor. Maps ALL degrees of freedom forward. */
	virtual void Predictor(void) = 0;

	/** corrector. Maps only the ACTIVE degrees of freedom forward.
	 * \param eqnos equations for the degrees of freedom of every node
	 * \param update vector of updates to the active degrees of freedom 
	 * \param eq_start lowest equation number to consider \a active.
	 * \param eq_stop highest equation number to consider \a active. */
	virtual void Corrector(const iArray2DT& eqnos, const dArrayT& update,
		int eq_start, int eq_stop) = 0;

	/** apply corrector to active equations with a node number map.
	 * \param map list of nodes corresponding to the rows of eqnos and update 
	 * \param eqnos equations for the degrees of freedom of every node
	 * \param update updates to the nodal field data
	 * \param eq_start lowest equation number to consider \a active.
	 * \param eq_stop highest equation number to consider \a active. */
	virtual void MappedCorrector(const iArrayT& map, const iArray2DT& eqnos,
		const dArray2DT& update, int eq_start, int eq_stop) = 0;

	/** return the field array needed by nIntegratorT::MappedCorrector. */
	virtual const dArray2DT& MappedCorrectorField(void) const = 0;
	
	/** prescribe the field and derivatives consistent BC's */
	virtual void ConsistentKBC(const KBC_CardT& KBC) = 0;

protected:  	
	
	/** recalculate time stepping constants */
	virtual void nComputeParameters(void) = 0;
	
protected:

	/** field data and derivatives */
	ArrayT<dArray2DT*> fU;
};

/* inlines */

#endif /* _N_CONTROLLERT_H_ */
