/* $Id: C0_QuadT.h,v 1.1 2004-09-02 18:25:08 rdorgan Exp $ */ 
#ifndef _C0_QUAD_T_H_
#define _C0_QUAD_T_H_

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

/* base class */
#include "ShapeTools.h"
#include "GeometryT.h"

namespace Tahoe {

class C0_QuadT: public ShapeTools, public GeometryT
{
 public:
	/** constructor */
	C0_QuadT(GeometryT::CodeT geometry_code, int numip, int numnd_u, int numdof_u, const LocalArrayT& coords);

	/** set all local parameters. call immediately after constructor */
	void Initialize(void);

	/** shape function derivatives.
	 * compute the derivatives of the shape functions with respect
	 * to the given coordinates by the chain rule for all integration
	 * points at once.
	 * \param coords nodal coordinates */
	void SetDerivatives();

	/** \name shape function values */
	/*@{*/
	/** shape functions defining the field at the given integration point */
	const double* IPShapeU(int ip) const { return fGNaU(ip); };

	/** gradient of the shape functions defining the field at the given integration point */
	const dArray2DT& IPDShapeU(int ip) const { return fGDNaU[ip]; };
	
	/** Laplacian of shape functions defining the field at the given integration point */
	const dArray2DT& IPDDShapeU(int ip) const { return fGDDNaU[ip]; };
	/*@}*/

 private: 
	/** compute local shape functions and derivatives at IPs for 1D. 
	 * The shape functions and their derivatives are evaluated for one 
	 * of the pre-defined integration rules. */
	void SetLocalShape();

	/** compute the jacobian of the nodal values.
	 * uses externally provided shape function derivatives.
	 * \param nodal values at the nodes: [nnd] x [ndim]
	 * \param LDNaX shape function derivatives: [ndim] x [nnd]
	 * \param jacobian resulting jacobian: [ndim] x [ndim] */
	void Jacobian(const LocalArrayT& nodal, const dArray2DT& LDNaX, dMatrixT& jac);
   
 private:
	/** local coordinates */
	const LocalArrayT& fCoords;

	/** \name dimensions */
	/*@{*/
	GeometryT::CodeT fGeometryCode; /**< geometry shape code */
	int fNumIP;                 /**< number of integration points */
	int fNumSD;                 /**< number of spatial dimensions */

	int fNumNodes_U;            /**< number of domain nodes */
	int fNumNodes_X;            /**< number of domain nodes */
	int fNumDOF_U;              /**< number of degrees of freedom for field */
	int fNumDOF_X;              /**< number of degrees of freedom for field */
	/*@}*/

	/** \name field shape functions */
	/*@{*/
	dArray2DT fLNaU;            /**< (#IP x #nodes#dof) (local, parent domain) */
	dArray2DT fGNaU;            /**< (#IP x #nodes#dof) (global, physical domain) */
	/*@}*/

	/** \name field shape function derivatives */
	/*@{*/
	ArrayT<dArray2DT> fLDNaU;   /**< [#IP x (#sd x #nodes#dof)] (local, parent domain) */
	ArrayT<dArray2DT> fGDNaU;   /**< [#IP x (#sd x #nodes#dof)] (global, physical domain) */
	/*@}*/

	/** \name field shape function Laplacians */
	/*@{*/
	ArrayT<dArray2DT> fLDDNaU;  /**< [#IP x (#sd x #nodes#dof)] (local, parent domain) */
	ArrayT<dArray2DT> fGDDNaU;  /**< [#IP x (#sd x #nodes#dof)] (global, physical domain) */
	/*@}*/

	/** \name geometry shape functions (parent domain) */
	/*@{*/
	dArray2DT fLNaX;            /**< (#IP x #nodes#dof) (local, parent domain) */
	/*@}*/

	/** \name geometry shape function derivatives */
	/*@{*/
	ArrayT<dArray2DT> fLDNaX;   /**< [#IP x (#sd x #nodes#dof)] (local, parent domain) */
	/*@}*/

	/** \name work space */
	/*@{*/
	//	dMatrixT fNodalExtrap;      /**< extrapolation matrix - least square smoothing (#nodes x #IP) */
	dArrayT  fWeights;          /**< integration weights */
	dArrayT  fArray1;           /**< (#sd) */
	dArrayT  fArray2;           /**< (#sd) */
	dMatrixT fJacobian;         /**< jacobian matrix (#sd x #sd) */
	dMatrixT fMatx1;            /**< (#sd x #sd) */
	/*@}*/
};

} // namespace Tahoe 
#endif /* _C0_QUAD_T_H_ */
