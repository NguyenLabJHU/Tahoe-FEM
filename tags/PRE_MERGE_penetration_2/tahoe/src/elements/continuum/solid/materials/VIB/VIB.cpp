/* $Id: VIB.cpp,v 1.13 2004-06-17 07:40:37 paklein Exp $ */
/* created: paklein (10/30/1997) */
#include "VIB.h"

#include <math.h>
#include <iostream.h>

#include "toolboxConstants.h"
#include "ExceptionT.h"

#include "ifstreamT.h"
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
#include "SF2.h"

using namespace Tahoe;

/* constructors */
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
			double A, B;
			in >> A >> B;
			fPotential = new LennardJones612(A,B);
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
		case C1FunctionT::kQuadraticPot:
		{
			double A, B;
			in >> A;		
	                in >> B;
			fPotential = new ParabolaT(A,B,1.0);
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
		case C1FunctionT::kSF2:
		{
			double A, B;
			in >> A >> B;		
			fPotential = new SF2(A,B);
			break;
		}
		default:		
			throw ExceptionT::kBadInputValue;	
	}
	if (!fPotential) throw ExceptionT::kOutOfMemory;
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
void VIB::Dimension(int numbonds)
{
	/* length table */
	fLengths.Dimension(numbonds);

	/* potential tables */
	fU.Dimension(numbonds);
	fdU.Dimension(numbonds);
	fddU.Dimension(numbonds);

	/* jacobian table */
	fjacobian.Dimension(numbonds);

	/* STRESS angle tables - by associated stress component */
	fStressTable.Dimension(fNumStress, numbonds);
	  	
	/* MODULI angle tables - using Cauchy symmetry */ 	
	fModuliTable.Dimension(fNumModuli, numbonds);	
}
