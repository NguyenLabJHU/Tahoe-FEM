/* $Id: C1FunctionT.cpp,v 1.3 2002-10-20 22:38:47 paklein Exp $ */
/* created: paklein (12/04/1996) */
#include "C1FunctionT.h"
#include "dArrayT.h"
#include <float.h>

using namespace Tahoe;

/* constructor */
C1FunctionT::C1FunctionT(void) { }

/* destructor */
C1FunctionT::~C1FunctionT(void) { }

/* returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime.  Default behavior is just to map the virtual functions
* above */
dArrayT& C1FunctionT::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw ExceptionT::kGeneralFail;
	
	double *pin  =  in.Pointer();
	double *pout = out.Pointer();
	int    length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = Function(*pin++);
		
	return out;	
}

dArrayT& C1FunctionT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw ExceptionT::kGeneralFail;
	
	double *pin  =  in.Pointer();
	double *pout = out.Pointer();
	int    length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = DFunction(*pin++);	
		
	return out;
}

dArrayT& C1FunctionT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw ExceptionT::kGeneralFail;
	
	double *pin  =  in.Pointer();
	double *pout = out.Pointer();
	int    length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = DDFunction(*pin++);	
		
	return out;
}

/* return 0th, 1st, and 2nd derivative in the respective
* fields of the dArrayT. Default behavior is just to call the
* virtual functions above */  	
void C1FunctionT::SetAll(double x, dArrayT& data) const
{
	data[0] = Function(x);
	data[1] = DFunction(x);
	data[2] = DDFunction(x);
}

/* function domain */
double C1FunctionT::DomainMin(void) const { return DBL_MIN; }
double C1FunctionT::DomainMax(void) const { return DBL_MAX; }
