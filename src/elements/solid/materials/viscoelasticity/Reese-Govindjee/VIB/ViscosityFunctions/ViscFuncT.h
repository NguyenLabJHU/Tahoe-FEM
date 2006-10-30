/* $Id: ViscFuncT.h,v 1.1 2006-10-30 23:32:05 thao Exp $ */
/* created: paklein (12/04/1996) */

#ifndef _VISC_FUNC_T_H_
#define _VISC_FUNC_T_H_

#include "Environment.h"
#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;

/** interface for a twice differentiable function */
class ViscFuncT
{
public:

	/* function codes of derived classes */
	enum TypesT {kConstant = 0,
	            kVariVisc = 1,
		     kBiQuadratic = 2};

	/* constructor */
	ViscFuncT(void);

	/* destructor */
	virtual ~ViscFuncT(void);
	
	/* I/O */
	virtual void Print(ostream& out) const = 0;     	    	   	
	virtual void PrintName(ostream& out) const = 0;     	    	   	
	    	   	    	
	/* returning values */
	virtual double Function(double Jv, double Je) const = 0;
	virtual double DFuncDJv(double Jv, double Je) const = 0;
	virtual double DFuncDJe(double Jv, double Je) const = 0;

};
}
#endif /* _VISC_FUNC_T_H_ */
