/* $Id: ConstantT.cpp,v 1.2 2011-12-01 20:38:12 beichuan Exp $ */
/* created: paklein (03/25/1999)                                          */

#include "ConstantT.h"
#include <iostream>
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

