/* $Id: VIB.cpp,v 1.5 2002-07-02 19:55:55 cjkimme Exp $ */
/* created: paklein (10/30/1997) */
#include "VIB.h"

#include <math.h>
#include <iostream.h>

#include "Constants.h"
#include "ExceptionCodes.h"

#include "fstreamT.h"
#include "dSymMatrixT.h"

/* potential functions */
#include "LennardJones612.h"
#include "SmithFerrante.h"
#include "GaoKlein.h"
#include "GaoJi.h"
#include "GaoJi2.h"
#include "GaoVicky.h"
#include "ParabolaT.h"
#include "Triantafyllidis.h"

/* constructors */

using namespace Tahoe;

VIB::VIB(ifstreamT& in, int nsd, int numstress, int nummoduli):
	fNumSD(nsd),
	fNumStress(numstress),
	fNumModuli(nummoduli)
{
	/* set potential function */
	int potentialcode;
	in >> potentialcode;	
	switch(potentialcode)
	{
		case C1FunctionT::kLennardJones:
		{	
			double A;
			in >> A;
			fPotential = new LennardJones612(A);
			break;
		}	
		case C1FunctionT::kSmithFerrante:
		{
			double A, B;
			in >> A >> B;		
			fPotential = new SmithFerrante(A,B);
			break;
		}
		case C1FunctionT::kGaoKlein:
		{
			double A, B, C;
			in >> A >> B >> C;		
			fPotential = new GaoKlein(A,B,C);
			break;
		}
		case C1FunctionT::kQuadratic:
		{
			double A;
			in >> A;		
			fPotential = new ParabolaT(A);
			break;
		}
		case C1FunctionT::kTriantafyllidis:
		{
			double A;
			in >> A;		
			fPotential = new Triantafyllidis(A);
			break;
		}
		case C1FunctionT::kGaoJi:
		{
			double A, B, C;
			in >> A >> B >> C;		
			fPotential = new GaoJi(A,B,C);
			break;
		}
		case C1FunctionT::kGaoJi2:
		{
			double A, B, C;
			in >> A >> B >> C;
			fPotential = new GaoJi2(A,B,C);
			break;
		}
		case C1FunctionT::kGaoVicky:
		{
			double A, B, C, D;
			in >> A >> B >> C >>D;
			fPotential = new GaoVicky(A,B,C,D);
			break;
		}
		default:
		
			throw eBadInputValue;	
	}
	if (!fPotential) throw eOutOfMemory;
}

/* destructor */
VIB::~VIB(void)
{
	delete fPotential;
}

/* print parameters */
void VIB::Print(ostream& out) const
{
	out << " Number of spatial dimensions. . . . . . . . . . = " << fNumSD << '\n';
	
	/* potential parameters */
	fPotential->Print(out);
}

void VIB::PrintName(ostream& out) const
{
	out << "    Virtual Internal Bond\n";
	
	/* potential name */
	fPotential->PrintName(out);
}

/*************************************************************************
* Protected
*************************************************************************/

/* allocate memory for all the tables */
void VIB::Allocate(int numbonds)
{
	/* length table */
	fLengths.Allocate(numbonds);

	/* potential tables */
	fU.Allocate(numbonds);
	fdU.Allocate(numbonds);
	fddU.Allocate(numbonds);

	/* jacobian table */
	fjacobian.Allocate(numbonds);

	/* STRESS angle tables - by associated stress component */
	fStressTable.Allocate(fNumStress, numbonds);
	  	
	/* MODULI angle tables - using Cauchy symmetry */ 	
	fModuliTable.Allocate(fNumModuli, numbonds);	
}
