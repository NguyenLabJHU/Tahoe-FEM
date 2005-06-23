/* $Id: SIERRA_HypoElasticT.h,v 1.3 2004-07-27 03:16:04 paklein Exp $ */
#ifndef _SIERRA_HYPO_ELASTIC_T_H_
#define _SIERRA_HYPO_ELASTIC_T_H_ 

/* base class */
#include "SIERRA_Material_BaseT.h"

namespace Tahoe {

/** hypoelastic material written through the Sierra materials interface. */
class SIERRA_HypoElasticT: public SIERRA_Material_BaseT
{
public:

	/** constructor */
	SIERRA_HypoElasticT(void);

	/** returns the strain energy density for the specified strain. This
	 * material write the strain energy density into the state variable
	 * array and returns the value here. */
	virtual double StrainEnergyDensity(void);

protected:

	/** state variable variable indicies */
	enum StateIndexT{
		kStrainEnergyDensity = 0,
		  kNumStateVariables = 1
	};

	/** \name required by sub-classes */
	/*@{*/
	/** call the SIERRA registration function */
	virtual void Register_SIERRA_Material(void) const;

	/** set material output */
	virtual void SetOutputVariables(iArrayT& variable_index,
		ArrayT<StringT>& output_labels) const;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _SIERRA_HYPO_ELASTIC_T_H_ */
