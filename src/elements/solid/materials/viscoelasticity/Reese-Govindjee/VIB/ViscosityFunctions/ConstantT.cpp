/* $Id: ConstantT.cpp,v 1.1 2006-10-30 23:32:05 thao Exp $ */
/* created: paklein (03/25/1999)                                          */

#include "ConstantT.h"
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructors */
ConstantT::ConstantT(double A): fA(A){ }

/* I/O */
void ConstantT::Print(ostream& out) const
{
	/* parameters */
       out <<"\n      A = "<< fA << '\n';
}

void ConstantT::PrintName(ostream& out) const
{
        out << "Constant"<<'\n';
}

