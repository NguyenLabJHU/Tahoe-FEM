/* $Id: LinearDamageT.h,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (08/26/2000)                                          */

#ifndef _LINEAR_DAMAGE_T_H_
#define _LINEAR_DAMAGE_T_H_

/* base class */
#include "DecohesionLawT.h"

/* forward declarations */
class ifstreamT;

class LinearDamageT: public DecohesionLawT
{
public:

	/* constructor */
	LinearDamageT(ifstreamT& in, const dArrayT& init_traction,
		iArrayT& i_store, dArrayT& d_store);

	/* surface potential */
	virtual double Potential(const dArrayT& jump_u);
	
	/* traction vector given displacement jump vector */	
	virtual const dArrayT& Traction(const dArrayT& jump_u);

	/* potential stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump_u);

	/* surface status */
	virtual StatusT Status(const dArrayT& jump_u);

	/* print parameters to the output stream */
	virtual void PrintName(ostream& out) const;
	virtual void Print(ostream& out) const;

	/* storage dimensions */
	virtual int IntegerStorage(void) const; // per ip
	virtual int DoubleStorage(void) const;  // per ip

	/* initialize surface */
	virtual void InitializeFacet(void); // ip at a time

	/* update/reset internal variables */
	virtual void UpdateHistory(void); // ip at a time
	virtual void ResetHistory(void);  // ip at a time
	
private:

	/* traction at initialization */
	const dArrayT& fInitTraction;

	/* traction potential parameters */
	double fd_c_n;     // characteristic normal opening to failure
	double fd_c_t;     // characteristic tangential opening to failure
	
	/* penetration stiffness */
	double fpenalty; // stiffening multiplier
	double fK;       // penetration stiffness
};

#endif /* _LINEAR_DAMAGE_T_H_ */
