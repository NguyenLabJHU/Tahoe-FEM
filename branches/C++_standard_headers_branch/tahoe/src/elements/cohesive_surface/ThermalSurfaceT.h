/* $Id: ThermalSurfaceT.h,v 1.6 2004-07-15 08:25:57 paklein Exp $ */

#ifndef _THERMAL_SURFACE_T_H_
#define _THERMAL_SURFACE_T_H_

/* base class */
#include "CSEBaseT.h"

/* direct members */
#include "pArrayT.h"

namespace Tahoe {

/* forward declarations */
class C1FunctionT;

/** surface with thermal resistance */
class ThermalSurfaceT: public CSEBaseT
{
public:

	/** constructor */
	ThermalSurfaceT(const ElementSupportT& support);

	/** form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** tangent matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT);

	/** force vector */
	virtual void RHSDriver(void);

	/** compute output values */
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
		const iArrayT& e_codes, dArray2DT& e_values);

private:

	/** conduction parameters */
	enum ParametersT {
		kK_0 = 0, /**< initial conduction */
		kd_c = 1  /**< opening to zero conductivity */
	};
	
private:

	/** \name conduction parameters
	 * This is temporary. There parameters will be generalized later. */
	/*@{*/
	dArray2DT fConduction;
	/*@}*/
	
	/** nodal temperatures */
	LocalArrayT fLocTemperatures;
};

} // namespace Tahoe 
#endif /* _THERMAL_SURFACE_T_H_ */
