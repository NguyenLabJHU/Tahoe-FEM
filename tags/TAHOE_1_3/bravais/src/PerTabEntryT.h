// DEVELOPMENT
/* $Id: PerTabEntryT.h,v 1.3 2003-05-04 23:17:30 paklein Exp $ */

#ifndef _PER_TAB_ENTRY_T_H_
#define _PER_TAB_ENTRY_T_H_

#include "dArrayT.h"
#include "StringT.h"

using namespace Tahoe;

class PerTabEntryT {
private:
	double atomicmass;
	StringT latticetype;	
	StringT name;
	StringT symbol; 
	dArrayT latt_params[3];

public:
	PerTabEntryT();
	~PerTabEntryT();

	void SetName(const char * s);
	void SetSymbol(const char * s);
	void SetLattParams(double a);
	void SetLattParams(double a, double c);
	void SetLattParams(double a, double b, double c);
	void SetMass(double m);
	void SetLattType(const char * s);

	const StringT& GetName();
	const StringT& GetSymbol();
	const dArrayT& GetLattParams();
	const double& GetMass();
	const StringT& GetLattType();

};

#endif
