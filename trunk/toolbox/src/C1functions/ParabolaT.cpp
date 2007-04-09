/* $Id: ParabolaT.cpp,v 1.9 2007-04-09 21:06:15 tdnguye Exp $ */
/* created: paklein (03/25/1999)                                          */

#include "ParabolaT.h"
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

/* constructors */

using namespace Tahoe;

ParabolaT::ParabolaT(double k, double B, double l0): fk(k), fl0(l0), fB(B) {SetName("parabola");}

ParabolaT::ParabolaT(void):
	fk(0.0), fB(0.0), fl0(0.0)
{
	SetName("parabola");
}
/* I/O */
void ParabolaT::Print(ostream& out) const
{
        /* parameters */
        out << " U'' . . . . . . . . . . . . . . . . . . . . . . = " << fk << '\n';
}

void ParabolaT::PrintName(ostream& out) const
{
        out << "    Quadratic function\n";
}

/* returning values in groups */
dArrayT& ParabolaT::MapFunction(const dArrayT& in, dArrayT& out) const
{
        /* dimension checks */
        if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

        const double* pl = in.Pointer();
        double* pU = out.Pointer();
        
        for (int i = 0; i < in.Length(); i++)
        {
	        double x = *pl-fl0;
                *pU++ = 0.5*fk*x*x-0.5*fk*fB;
                pl++;
        }

        return out;
}

dArrayT& ParabolaT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
        /* dimension checks */
        if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

        const double* pl  = in.Pointer();
        double* pdU = out.Pointer();
        
        for (int i = 0; i < in.Length(); i++)
	{
	        double x = *pl-fl0;
                *pdU++ = fk*x;
		pl++;
		
	}
        return out;
}

dArrayT& ParabolaT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
        /* dimension checks */
        if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

        out = fk;
        return out;
}

void ParabolaT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);

	list.AddParameter(fk, "k");
	list.AddParameter(fB, "b");
	list.AddParameter(fl0, "l0");
	
	/* set the description */
	list.SetDescription("f(x) = 0.5*k*(x-l0)^2 - 0.5*k*b");	
}

void ParabolaT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	fk = list.GetParameter("k");
	fB = list.GetParameter("b");
	fl0 = list.GetParameter("l0");
}

