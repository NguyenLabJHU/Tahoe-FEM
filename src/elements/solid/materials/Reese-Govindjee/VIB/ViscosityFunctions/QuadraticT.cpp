/* $Id: QuadraticT.cpp,v 1.1 2003-03-19 19:00:58 thao Exp $ */
/* created: paklein (03/25/1999)                                          */

#include "QuadraticT.h"
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructors */
QuadraticT::QuadraticT(double A, double B, double C): 
fA(A),
fB(B),
fC(C)
 { }

/* I/O */
void QuadraticT::Print(ostream& out) const
{
	/* parameters */
	out <<"A: "<<fA<< " ...... B: "<<fB<< " ...... C: "<<fC<<'\n';
}

void QuadraticT::PrintName(ostream& out) const
{
        out <<"Function ....... fA*(J-fB)^2+fC"<<'\n';
}

