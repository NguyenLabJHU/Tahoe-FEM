/* $Id: ShapeTools.h,v 1.1 2004-09-02 18:25:08 rdorgan Exp $ */ 
#ifndef _SHAPE_TOOLS_H_
#define _SHAPE_TOOLS_H_

#include <iostream.h>
#include <ctype.h>

#include "ifstreamT.h"
#include "StringT.h"
#include "ExceptionT.h"
#include "toolboxConstants.h"

#include "LocalArrayT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

namespace Tahoe {

class ShapeTools
{
 public:
	/** constructor */
	ShapeTools();

	/** destructor */
	virtual ~ShapeTools(void);

	/** set all local parameters. call immediately after constructor */
	virtual void Initialize(void) = 0;

	/** shape function derivatives.
	 * compute the derivatives of the shape functions with respect
	 * to the given coordinates by the chain rule for all integration
	 * points at once.
	 * \param coords nodal coordinates */
	virtual void SetDerivatives() = 0;

	/** \name shape function values */
	/*@{*/
	/** shape functions defining the field at the given integration point */
	virtual const double* IPShapeU(int ip) const = 0;

	/** gradient of the shape functions defining the field at the given integration point */
	virtual const dArray2DT& IPDShapeU(int ip) const = 0;
	
	/** Laplacian of shape functions defining the field at the given integration point */
	virtual const dArray2DT& IPDDShapeU(int ip) const = 0;
	/*@}*/

 private: 
	/** compute local shape functions and derivatives at IPs for 1D. 
	 * The shape functions and their derivatives are evaluated for one 
	 * of the pre-defined integration rules. */
	virtual void SetLocalShape() = 0;

	/** compute the jacobian of the nodal values.
	 * uses externally provided shape function derivatives.
	 * \param nodal values at the nodes: [nnd] x [ndim]
	 * \param LDNaX shape function derivatives: [ndim] x [nnd]
	 * \param jacobian resulting jacobian: [ndim] x [ndim] */
	virtual void Jacobian(const LocalArrayT& nodal, const dArray2DT& LDNaX, dMatrixT& jac) = 0;
   
 private:
};

} // namespace Tahoe 
#endif /* _SHAPE_TOOLS_H_ */
