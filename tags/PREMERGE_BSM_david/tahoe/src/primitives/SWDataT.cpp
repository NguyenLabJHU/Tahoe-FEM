/* $Id: SWDataT.cpp,v 1.6 2004-07-15 08:31:36 paklein Exp $ */
/* created: paklein (03/22/1997) */
#include "SWDataT.h"

#include "ifstreamT.h"
#include <math.h>
#include <iostream.h>
#include <iomanip.h>

using namespace Tahoe;

/* constructor */
SWDataT::SWDataT(void): 
	ParameterInterfaceT("Stillinger-Weber"),
	feps(0.0), fA(0.0), fdelta(0.0), fgamma(0.0),
	flambda(0.0), frcut(0.0), fa(0.0), fB(0.0)
{

}

SWDataT::SWDataT(ifstreamT& in):
	ParameterInterfaceT("Stillinger-Weber")
{
	Read(in);
}

/* I/O operators */
void SWDataT::Read(ifstreamT& in)
{
	/* unit scaling */
	in >> feps;	if (feps <= 0.0) throw ExceptionT::kBadInputValue;

	/* 2 body potential */
	in >> fA;		if (fA     <= 0.0) throw ExceptionT::kBadInputValue;
	in >> fdelta;	if (fdelta <= 0.0) throw ExceptionT::kBadInputValue;
	
	/* 3 body potential */
	in >> fgamma;	if (fgamma  <= 0.0) throw ExceptionT::kBadInputValue;
	in >> flambda;	if (flambda <= 0.0) throw ExceptionT::kBadInputValue;
	
	in >> frcut;	if (frcut <= 0.0) throw ExceptionT::kBadInputValue;		
	in >> fa;		if (fa    <= 0.0) throw ExceptionT::kBadInputValue;

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

/* describe the parameters needed by the interface */
void SWDataT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	list.AddParameter(feps, "epsilon");
	list.AddParameter(fA, "A");
	list.AddParameter(fdelta, "delta");
	list.AddParameter(fgamma, "gamma");
	list.AddParameter(flambda, "lambda");
	list.AddParameter(frcut, "r_cut");
	list.AddParameter(fa, "a");
}
 
/* accept parameter list */
void SWDataT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	feps    = list.GetParameter("epsilon");
	fA      = list.GetParameter("A");
	fdelta  = list.GetParameter("delta");
	fgamma  = list.GetParameter("gamma");
	flambda = list.GetParameter("lambda");
	frcut   = list.GetParameter("r_cut");
	fa      = list.GetParameter("a");

	/* compute B factor */
	double a0 = pow(2.0,1.0/6.0);
	fB =-(fdelta*pow(a0,5))/(-(a0*fdelta) - 4.0*a0*a0 +
	           8.0*a0*frcut - 4.0*frcut*frcut);
}
