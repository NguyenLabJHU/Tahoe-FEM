/* $Id: VariViscT.cpp,v 1.2 2002-10-20 22:48:48 paklein Exp $ */
/* created: paklein (03/25/1999)                                          */

#include "VariViscT.h"
#include "ExceptionT.h"

using namespace Tahoe;

/* constructors */
VariViscT::VariViscT(double no,double d,double z,double Jcrit,double Jo): 
	fno(no),
	fdo(d),
	fz(z),
	fJcrit(Jcrit),
	fJo(Jo){ }

/* I/O */
void VariViscT::Print(ostream& out) const
{
	/* parameters */
        out <<"\n         no = " << fno; 
	out <<"\n         do = " << fdo; 
	out <<"\n          z = " << fz;  
        out <<"\n        Jcr = " << fJcrit  << '\n';
        out <<"\n         Jo = " << fJo  << '\n';
} 

void VariViscT::PrintName(ostream& out) const
{
	out <<"Exponentially decaying viscosity function\n";
	/*        out <<"n(lv, Je) = no exp( -(Jv-Jo)/d(Je) )\n";
		  out <<"d(Je) = do (1 + exp( (Je - Jcr) / z) )\n";   */
}

