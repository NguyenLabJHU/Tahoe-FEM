/* $Id: BiQuadraticT.cpp,v 1.1 2006-10-30 23:32:05 thao Exp $ */
/* created: paklein (03/25/1999)                                          */

#include "BiQuadraticT.h"
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructors */
BiQuadraticT::BiQuadraticT(double A1, double A2, double B, double C): 
	fA1(A1),
	fA2(A2),
	fB(B),
	fC(C)
 { }

/* I/O */
void BiQuadraticT::Print(ostream& out) const
{
	/* parameters */
        out <<"\n       A1 = " << fA1; 
	out <<"\n       A2 = " << fA2; 
	out <<"\n        B = " << fB;  
        out <<"\n        C = " << fC  << '\n';
} 

void BiQuadraticT::PrintName(ostream& out) const
{
	out << "Bi-Quadratic \n";
}

