/* $Id: ViscVIB.h,v 1.1 2006-10-30 23:32:05 thao Exp $ */
/* created: TDN (1/19/2000) */
 
#ifndef _VISCVIB_H_
#define _VISCVIB_H_

/* direct members */
#include "dArrayT.h"
#include "dArray2DT.h"
#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class ViscFuncT;
class C1FunctionT;

/** Base class for isotropic Viscoelastic VIB solvers */
class ViscVIB
{
	public:

		/*constructor*/
		ViscVIB(ifstreamT& in, int nsd, int numstress, int nummoduli);
		
		/*destructor*/
		virtual ~ViscVIB(void);
		
		/*print parameters*/
		virtual void Print(ostream& out) const;
		virtual void PrintName(ostream& out) const;
		
	protected:
		
		/*allocate memory for all the tables*/
		void Allocate(int numbonds);
		
	protected:
	
		/*number of spatial dimensions*/
		int fNumSD;
		
		/*elastic potential function*/
		C1FunctionT* fPotential_E;
		
		/*inelastic potential function*/
		C1FunctionT* fPotential_I;
		
		/*viscosity function*/
		ViscFuncT* fShearVisc;
		ViscFuncT* fBulkVisc;

		/*length table*/
		dArrayT fLengths_E;
		dArrayT fLengths_I;    	
		
		/*potential tables*/
		dArrayT fU_E;
		dArrayT fdU_E;
		dArrayT fddU_E;
		
		dArrayT fU_I;
		dArrayT fdU_I;
		dArrayT fddU_I;
		
		/*jacobian table*/
		dArrayT fjacobian;
		
		/*Stress and Moduli angle tables - by associated principal components       */
		/*Note: Only one stress table is used to calculate the elastic and inelastic*
		 *parts of stress and modulus because the material in the intermediate and  *
		 *underformed reference frame are both isotropic.                           */
		  int fNumStress;
		  dArray2DT fStressTable;
		  
		  int fNumModuli;
		  dArray2DT fModuliTable;
};
}
#endif /*_VISCVIB_H_*/
