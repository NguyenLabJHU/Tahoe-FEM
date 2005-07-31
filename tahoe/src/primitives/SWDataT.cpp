/* $Id: SWDataT.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (03/22/1997)                                          */
/* Container class for Stillinger-Weber potential parameters              */

#include "SWDataT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "Constants.h"
#include "ExceptionCodes.h"

#include "fstreamT.h"

/* constructor */
SWDataT::SWDataT(void): feps(0.0), fA(0.0), fdelta(0.0), fgamma(0.0),
	flambda(0.0), frcut(0.0), fa(0.0), fB(0.0)
{

}

SWDataT::SWDataT(ifstreamT& in)
{
	Read(in);
}

/* I/O operators */
void SWDataT::Read(ifstreamT& in)
{
	/* unit scaling */
	in >> feps;	if (feps <= 0.0) throw eBadInputValue;

	/* 2 body potential */
	in >> fA;		if (fA     <= 0.0) throw eBadInputValue;
	in >> fdelta;	if (fdelta <= 0.0) throw eBadInputValue;
	
	/* 3 body potential */
	in >> fgamma;	if (fgamma  <= 0.0) throw eBadInputValue;
	in >> flambda;	if (flambda <= 0.0) throw eBadInputValue;
	
	in >> frcut;	if (frcut <= 0.0) throw eBadInputValue;		
	in >> fa;		if (fa    <= 0.0) throw eBadInputValue;

	/* compute B factor */
	double a0 = pow(2.0,1.0/6.0);
	fB =-(fdelta*pow(a0,5))/(-(a0*fdelta) - 4.0*a0*a0 +
	           8.0*a0*frcut - 4.0*frcut*frcut);
}

void SWDataT::Write(ostream& out) const
{
	out << "\n Stillinger-Weber Parameters:\n";

	/* unit scaling */
	out << "    epsilon = " << feps << '\n';		

	/* 2 body potential */
	out << " 2 body terms:\n";
	out << "          A = " << fA << '\n';
	out << "          B = " << fB << "   **COMPUTED\n";
	out << "      delta = " << fdelta << '\n';
	
	/* 3 body potential */
	out << " 3 body terms:\n";
	out << "      gamma = " << fgamma << '\n';
	out << "     lambda = " << flambda << '\n';
	
	out << " Lattice scaling and cut-off terms:\n";
	out << "      r_cut = " << frcut << '\n';
	out << "         a0 = " << fa << '\n';
}
