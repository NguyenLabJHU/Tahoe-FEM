#ifndef _N_GEAR_04_H_
#define _N_GEAR_04_H_

/* base class */
#include "Gear4.h"
#include "nIntegratorT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

namespace Tahoe {

/** Node controller for an explicit 6th order accurate Gear time integration
 * algorithm. The BasicFieldT used with this integrator must have
 * BasicFieldT::Order of 6. */
class nGear4: public virtual Gear4, public nIntegratorT
{
public:

	/** constructor */
	nGear4(void);

	/** consistent BC's */
	virtual void ConsistentKBC(BasicFieldT& field, const KBC_CardT& KBC);

	/** pseudo-boundary conditions for external nodes */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const;

	/** predictor. Maps ALL degrees of freedom forward. */
	virtual void Predictor(BasicFieldT& field);

	/** corrector. Maps ALL degrees of freedom forward. */
	virtual void Corrector(BasicFieldT& field, const dArray2DT& update);

	/** corrector - map ACTIVE. See nIntegratorT::Corrector for more
	 * documentation */
	virtual void Corrector(BasicFieldT& field, const dArrayT& update, 
		int eq_start, int num_eq);

	/** corrector with node number map - map ACTIVE. See 
	 * nIntegratorT::MappedCorrector for more documentation */
	virtual void MappedCorrector(BasicFieldT& field, const iArrayT& map, 
		const iArray2DT& flags, const dArray2DT& update);

	/** return the field array needed by nIntegratorT::MappedCorrector. */
	virtual const dArray2DT& MappedCorrectorField(BasicFieldT& field) const;

private:

	/** \name Gear constants */
	/*@{*/
	double F02; 
	double F12;
	double F22;
	double F32;
	/*@}*/
	
	/** \name timestep-scaled Gear constants */
	/*@{*/
	double f02;
	double f12;
	double f22;
	double f32;
	/*@}*/
	
	/** \name Taylor expansion factors */
	/*@{*/
	double fdt2; /**< \f$ \frac{\Delta t^2}{2!} \f$ */
	double fdt3; /**< \f$ \frac{\Delta t^3}{3!} \f$ */
	/*@}*/
	
	/** \name Volume "field" from AndersenPressureT */
	/*@{*/
	dArrayT v_field;
	/*@}*/

 protected:

	/** recalculate time stepping constants */
	virtual void nComputeParameters(void);

};

} // namespace Tahoe

#endif /* _N_GEAR_06_H_ */
