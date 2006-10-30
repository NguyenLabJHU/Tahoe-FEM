/* $Id: Linear.cpp,v 1.1 2006-10-30 23:32:05 thao Exp $ */
/* created: paklein (03/25/1999)                                          */

#include "Linear.h"
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructors */
Linear::Linear(double A, double B): 
  fA(A),
  fB(B)
{ }

/* I/O */
void Linear::Print(ostream& out) const
{
	/* parameters */
	out <<"A: "<< fA << '\n';
	out <<"B: "<< fB << '\n';
}

void Linear::PrintName(ostream& out) const
{
        out << "Function .............  linear "<<'\n';
	out << "\tn(J) = no * J\n";
}


