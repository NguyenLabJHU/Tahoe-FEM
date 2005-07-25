/* $Id: nLinearHHTalpha.h,v 1.11.22.1 2005-07-25 02:37:16 paklein Exp $ */
/* created: paklein (10/14/1996) */
#ifndef _N_LINEARHHT_A_H_
#define _N_LINEARHHT_A_H_

/* base class */
#include "HHTalpha.h"
#include "nIntegratorT.h"

/* direct members */
#include "dArray2DT.h"

namespace Tahoe {

/** HHT alpha integration for linear systems */
class nLinearHHTalpha: virtual public HHTalpha, public nIntegratorT
{
public:

	/** constructor */
	nLinearHHTalpha(double alpha);

	/** declare field dimensions to the integrator */
	virtual void Dimension(const BasicFieldT& field);

	/** prescribe the field and derivatives consistent BC's 
	 * \param field field to which boundary conditions will be applied
	 * \param KBC boundary condition specification
	 * \param nodes pointer to a list of nodes if KBC card has KBC_CardT::ModeT set to KBC_CardT::kSet
	 *        or NULL if it is set to KBC_CardT::kNode. */
	virtual void ConsistentKBC(BasicFieldT& field, const KBC_CardT& KBC, const iArrayT* nodes = NULL);

	/** pseudo-boundary conditions for external nodes */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const;

	/** predictor. Maps ALL degrees of freedom forward Unless specified otherwise */
	virtual void Predictor(BasicFieldT& field, int fieldstart = 0, int fieldend = -1);

	/** corrector. Maps ALL degrees of freedom forward. */
	virtual void Corrector(BasicFieldT& field, const dArray2DT& update, int fieldstart = 0, int fieldend = -1);

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

protected:  	
	
	/* recalculate time stepping constants */
	virtual void nComputeParameters(void);

private:

	/* values from t_n (needed for HHT-alpha) */
	dArray2DT	dn;
	dArray2DT	vn;

	/* predictor */
	double	dpred_v;
	double	dpred_a;
	double	vpred_a;
	
	/* corrector */
	double	dcorr_d;
	double	dcorr_dpred;
	double	dcorr_a;
	
	double	vcorr_v;
	double	vcorr_vpred;
	double	vcorr_a;
	  	
	/* consistent BC */
	double	dalpha_a;
	double	valpha_a;

	//TEMP - only support one field
	const BasicFieldT* fField;
};

} // namespace Tahoe 
#endif /* _N_LINEARHHT_A_H_ */
