/* $Id: nExplicitCD.h,v 1.2 2001-08-27 17:12:13 paklein Exp $ */
/* created: paklein (03/23/1997) */

#ifndef _N_EXP_CD_H_
#define _N_EXP_CD_H_

/* base class */
#include "ControllerT.h"
#include "nControllerT.h"

/** Node controller for an explicit 2nd order accurate, central 
 * difference time-stepping algorithm. */
class nExplicitCD: public virtual ControllerT, public nControllerT
{
public:

	/** constructor */
	nExplicitCD(void);

	/** consistent BC's - updates predictors and acceleration only */
	virtual void ConsistentKBC(const KBC_CardT& KBC);
	
	/** predictor - map ALL */
	virtual void Predictor(void);

	/** corrector - map ACTIVE. See nControllerT::Corrector for more
	 * documentation */
	virtual void Corrector(const iArray2DT& eqnos, const dArrayT& update,
		int eq_start, int eq_stop);

	/** corrector with node number map - map ACTIVE. See 
	 * nControllerT::MappedCorrector for more documentation */
	virtual void MappedCorrector(const iArrayT& map, const iArray2DT& eqnos,
		const dArray2DT& update, int eq_start, int eq_stop);

	/** return the field array needed by nControllerT::MappedCorrector. */
	virtual const dArray2DT& MappedCorrectorField(void) const;

	/** pseudo-boundary conditions for external nodes */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const;

protected:  	
	
	/* recalculate time stepping constants */
	virtual void nComputeParameters(void);
	
private:

	/* predictor */
	double	dpred_v;
	double	dpred_a;
	double	vpred_a;
	
	/* corrector */  	
	double	vcorr_a; //also used for consistent BC
	  	  	
};

#endif /* _N_EXP_CD_H_ */
