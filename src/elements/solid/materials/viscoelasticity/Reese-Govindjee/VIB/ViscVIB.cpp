/* $Id: ViscVIB.cpp,v 1.1 2006-10-30 23:32:05 thao Exp $ */
/* created: TDN (1/19/2000) */

#include <math.h>
#include <iostream.h>

#include "ViscVIB.h"
#include "toolboxConstants.h"
#include "ExceptionT.h"

#include "ifstreamT.h"

/* potential functions */
#include "SmithFerrante.h"
#include "SF2.h"
#include "ParabolaT.h"
#include "VariViscT.h"
#include "ConstantT.h"
#include "BiQuadraticT.h"

using namespace Tahoe;

//see 	static int NumValues(int nsd) in dSymMatrixT.h

ViscVIB::ViscVIB(ifstreamT& in, int nsd, int numstress, int nummoduli):
	fNumSD(nsd),
	fNumStress(numstress),
	fNumModuli(nummoduli)
{
	/*set potential function*/
	int potentialcode;
	in >> potentialcode;	
	switch(potentialcode)
	{
		case C1FunctionT::kSmithFerrante:
		{
			double AE, BE, AIN, BIN;
			in >> AE >> BE;		
			fPotential_E = new SmithFerrante(AE,BE);
			in >> AIN >> BIN;
			fPotential_I = new SmithFerrante(AIN,BIN);
			break;
		}
		case C1FunctionT::kSF2:
		{
			double AE, BE, AIN, BIN;
			in >> AE >> BE;		
			fPotential_E = new SF2(AE,BE);
			in >> AIN >> BIN;
			fPotential_I = new SF2(AIN,BIN);
			break;
		}
	        case C1FunctionT::kQuadraticPot:
		{
		        double AE, AIN, BE,BIN;
			in >> AE >> BE;        
			fPotential_E = new ParabolaT(AE, BE,1.0);
			in >> AIN >> BIN;
			fPotential_I = new ParabolaT(AIN, BIN,1.0);
			break;
		}
		default:
		
			throw ExceptionT::kBadInputValue;	
	}
	in >> potentialcode;	
	switch(potentialcode)
	{
		case ViscFuncT::kConstant:
		{
			double A;
			in >> A;		
			fShearVisc = new ConstantT(A);
			in >> A;		
			fBulkVisc = new ConstantT(A);
			break;
		}	
		case ViscFuncT::kVariVisc:
		{
		        double n0, d0,zed,Jcrit;
			in >> n0 >> d0 >> zed >> Jcrit;
			fShearVisc = new VariViscT(n0, d0, zed, Jcrit);
			in >> n0 >> d0 >> zed >> Jcrit;
			fBulkVisc = new VariViscT(n0, d0, zed, Jcrit);
			break;
		}
		case ViscFuncT::kBiQuadratic:
		{
		        double A1,A2,B,C;
			in >> A1 >> A2 >> B >> C;
			fShearVisc = new BiQuadraticT(A1, A2, B, C);
			in >> A1 >> A2 >> B >> C;
			fBulkVisc = new BiQuadraticT(A1, A2, B, C);
			break;
		}
        	default:
		
			throw ExceptionT::kBadInputValue;	
	}
	
	/*set viscosity function*/
	if (!fPotential_E || !fPotential_I || !fShearVisc || !fBulkVisc) 
		throw ExceptionT::kOutOfMemory;
}

ViscVIB::~ViscVIB(void)
{
	delete fPotential_E;
	delete fPotential_I;
	delete fShearVisc;
	delete fBulkVisc;
}
/* print parameters */
void ViscVIB::Print(ostream& out) const
{
	out << " Number of spatial dimensions. . . . . . . . . . = " << fNumSD << '\n';
	
	/* potential parameters */
	out << " Elastic potential";
	fPotential_E->Print(out);

	out << " Inelastic potential";
	fPotential_I->Print(out);

	out << " Shear viscosity function";
	fShearVisc->Print(out);

	out << " Bulk Viscosity function ";
	fBulkVisc->Print(out);

}

void ViscVIB::PrintName(ostream& out) const
{
	out << "    Viscoelastic Virtual Internal Bond\n";
	
	/* potential name */
	out << "    Elastic potential:           ";
	fPotential_E -> PrintName(out);

	out << "    Inelastic potential:         ";
	fPotential_I -> PrintName(out);
	
	out << "    Shear viscosity function:    ";
	fShearVisc -> PrintName(out);

	out << "     Bulk viscosity function:    ";
	fBulkVisc -> PrintName(out);
}


/*************************************************************************
 * Protected
 *************************************************************************/

/* allocate memory for all the tables */
void ViscVIB::Allocate(int numbonds)
{
  	/* length table */
	fLengths_E.Dimension(numbonds);
	fLengths_I.Dimension(numbonds);

	/* potential tables */
	fU_E.Dimension(numbonds);
	fdU_E.Dimension(numbonds);
	fddU_E.Dimension(numbonds);

	fU_I.Dimension(numbonds);
	fdU_I.Dimension(numbonds);
	fddU_I.Dimension(numbonds);

  	/* jacobian table */
	fjacobian.Dimension(numbonds);
  
  	/* STRESS angle tables - by associated stress component */
	fStressTable.Dimension(fNumStress, numbonds);
  	  	
  	/* MODULI angle tables - using Cauchy symmetry */ 	
	fModuliTable.Dimension(fNumModuli, numbonds);	
} 
