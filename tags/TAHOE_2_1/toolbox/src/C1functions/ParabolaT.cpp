/* $Id: ParabolaT.cpp,v 1.8 2003-11-21 22:41:27 paklein Exp $ */
/* created: paklein (03/25/1999)                                          */

#include "ParabolaT.h"
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

/* constructors */

using namespace Tahoe;

ParabolaT::ParabolaT(double k, double B, double l0): fk(k), fl0(l0), fB(B) { }

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
