/* $Id: EAM_particle.cpp,v 1.2 2004-04-09 02:02:58 hspark Exp $ */
/* created: hspark(02/25/2004) */
#include "EAM_particle.h"
#include <iostream.h> //TEMP
#include "CBLatticeT.h"

/* EAM property types */
#include "ParadynEAMT.h"

using namespace Tahoe;

/* constructor */
EAM_particle::EAM_particle(CBLatticeT& lattice): fLattice(lattice),
	fCounts( fLattice.BondCounts() ), fBonds( fLattice.DeformedLengths() ),
	fNumSpatialDim( fLattice.NumberOfSpatialDim() ),
	fNumBonds( fLattice.NumberOfBonds() ),
	fModuliDim(dSymMatrixT::NumValues(fNumSpatialDim)),
	fBondTensor4(fModuliDim), fAmn(fNumBonds),
	fBondTensor2(fModuliDim), //fBondTensor2b(fModuliDim),
	fTensor2Table(fNumBonds,fModuliDim),
	fBond1(fNumBonds), fBond2(fNumBonds), fBond3(fNumBonds),
	fEAMProperty(NULL),
	fPairEnergy(NULL),
	fPairForce(NULL),
	fPairStiffness(NULL),
	fEmbedEnergy(NULL),
	fEmbedForce(NULL),
	fEmbedStiffness(NULL),
	fEDEnergy(NULL),
	fEDForce(NULL),
	fEDStiffness(NULL)
{
	/* dimension checks */
	if (fCounts.Length() != fNumBonds ||
	     fBonds.Length() != fNumBonds) throw ExceptionT::kGeneralFail;
}

/* Destructor */
EAM_particle::~EAM_particle(void) 
{
	delete fEAMProperty;
}

/* Set "glue" functions */
void EAM_particle::SetGlueFunctions(const StringT& param_file)
{
	/* construct EAM property - only the Paradyn EAM potentials are implemented */
	fEAMProperty = new ParadynEAMT(param_file);

	/* cache function pointers */
	fPairEnergy    = fEAMProperty->getPairEnergy();
	fPairForce     = fEAMProperty->getPairForce();
	fPairStiffness = fEAMProperty->getPairStiffness();
	fEDEnergy      = fEAMProperty->getElecDensEnergy();
	fEDForce       = fEAMProperty->getElecDensForce();
	fEDStiffness   = fEAMProperty->getElecDensStiffness();
	fEmbedEnergy   = fEAMProperty->getEmbedEnergy();
	fEmbedForce    = fEAMProperty->getEmbedForce();
	fEmbedStiffness = fEAMProperty->getEmbedStiffness();
	
	/* set lattice parameter */
	fLatticeParameter = fEAMProperty->GetLatticeParameter();
}

/*
* Compute unit strain energy density:
*
*     unit strain energy = energy/atom
*/
double EAM_particle::ComputeUnitEnergy(void)
{
	double rho = 0.0;
	double energy = 0.0;
	
	const int* pcount = fCounts.Pointer();
	for (int i = 0; i < fNumBonds; i++)
	{
		int ci = *pcount++;

		double  ri = fBonds[i];
		double a = fPairEnergy(ri, NULL, NULL);
		double phi = a*a/ri;
		double rho = fEDEnergy(ri, NULL, NULL);

		rho    += ci*rho;
		energy += ci*0.5*phi;
	}
	
	energy += fEmbedEnergy(rho, NULL, NULL);
	return energy;
}

/*
* Compute unit 2nd PK stress:
*
*     unit 2nd PK stress = SIJ*(volume per cell/atoms per cell)
*/
void EAM_particle::ComputeUnitStress(dSymMatrixT& stress)
{
	/* total atomic density */
	double rho = TotalElectronDensity();
	double dFdrho = fEmbedForce(rho, NULL, NULL);	// checks out

	/* assemble stress */
	stress = 0.0;
	
	for (int i = 0; i < fNumBonds; i++)	// fNumBonds, ri, ci check out
	{
		double ri = fBonds[i];
		int    ci = fCounts[i];	
		double a = fPairEnergy(ri, NULL, NULL);
		double Potential = a*a/ri;
		double DPotential = 2.0*fPairForce(ri, NULL, NULL)*a/ri-Potential/ri;
		double DDensity = fEDForce(ri, NULL, NULL);	
		double coeff = (1.0/ri)*ci*(0.5*DPotential + dFdrho*DDensity);
		fLattice.BondComponentTensor2(i,fBondTensor2);
		stress.AddScaled(coeff,fBondTensor2);
	}
}
	   	    	
/*
* Compute unit material tangent moduli:
*
*   unit material tangent moduli = CIJKL*(volume per cell/atoms per cell)
*/
void EAM_particle::ComputeUnitModuli(dMatrixT& moduli)
{
	/* compute total electron density */
	double rho = TotalElectronDensity();

	/* initialize */
	moduli = 0.0;

	/* mixed bond contribution */
	FormMixedBondContribution(rho, moduli);
	
	/* single bond contribution */
	FormSingleBondContribution(rho, moduli);
}

/* lattice parameter */
double EAM_particle::LatticeParameter(void) const
{
	return fLatticeParameter;
}

/*
* Compute the total electron density.
*/
double EAM_particle::TotalElectronDensity(void)
{
	double rho = 0.0;
	const int* pcount = fCounts.Pointer();

	for (int i = 0; i < fNumBonds; i++)
	{
		double ri = fBonds[i];
		double pedensity = fEDEnergy(ri, NULL, NULL);
		rho += (*pcount++)*(pedensity);
	}
	return rho;
}

/* return the embedding force for a given electron density */
double EAM_particle::ReturnEmbeddingForce(double rho)
{
	return fEmbedForce(rho, NULL, NULL);
}

/**********************************************************************
* Private
**********************************************************************/

/*
* Form matrix of mixed pair potential and embedding
* energy derivatives.  NOTE: computes the UPPER triangle
* ONLY.
*/
void EAM_particle::FormMixedDerivatives(double rho)
{
	double dFdrho = fEmbedForce(rho, NULL, NULL);
	double d2Fdrho2 = fEmbedStiffness(rho, NULL, NULL);

	/* form upper triangle only */
	for (int j = 0; j < fNumBonds; j++)
	{
		double rj = fBonds[j];
		int    cj = fCounts[j];

		double a = fPairEnergy(rj, NULL, NULL);
		double b = fPairForce(rj, NULL, NULL);
		double Potential = a*a/rj;
		double DPotential = 2.0*b*a/rj-Potential/rj;
		double DDPj = 2.0*fPairStiffness(rj, NULL, NULL)*a/rj;
		DDPj+=2.0*b*b/rj;
		DDPj-=(2.0*DPotential/rj);
		double DDDj = fEDStiffness(rj, NULL, NULL);
		double DDj = fEDForce(rj, NULL, NULL);
	
		for (int i = 0; i <= j; i++)
		{
			double Amn = 0.0;
			double ri = fBonds[i];
			int    ci = fCounts[i];
					
			/* diagonal terms */
			if (i == j)
			{
				/* pair potential */
				Amn += 0.5*cj*DDPj;	
				
				/* embedding energy */
				Amn += cj*dFdrho*DDDj;
			}
			double DDensity = fEDForce(ri, NULL, NULL);
		
			/* mixed embedding energy term */
			Amn += ci*cj*d2Fdrho2*DDensity*DDj;
		
			fAmn(i,j) = Amn/(ri*rj);
		}
	}
	
	fAmn.CopySymmetric();
}	

/*
* Moduli tensor contributions.
*/
void EAM_particle::FormSingleBondContribution(double rho, dMatrixT& moduli)
{
	/* Embedding energy derivative */
	double dFdrho = fEmbedForce(rho, NULL, NULL);
	
	for (int i = 0; i < fNumBonds; i++)
	{
		double ri = fBonds[i];
		double a = fPairEnergy(ri, NULL, NULL);
		double Potential = a*a/ri;
		double DPotential = 2.0*fPairForce(ri, NULL, NULL)*a/ri-Potential/ri;
		double DDensity = fEDForce(ri, NULL, NULL);
		double coeff = -fCounts[i]*(0.5*DPotential + dFdrho*DDensity)/(ri*ri*ri);
		fLattice.BondComponentTensor4(i,fBondTensor4);		
		moduli.AddScaled(coeff,fBondTensor4);
	}
}

void EAM_particle::FormMixedBondContribution(double rho, dMatrixT& moduli)
{
	/* batch fetch */
	fLattice.BatchBondComponentTensor2(fTensor2Table);

	/* mixed bond derivatives */
	FormMixedDerivatives(rho);

/*	
	for (int n = 0; n < fNumBonds; n++)
	{
		for (int m = 0; m < fNumBonds; m++)
		{
			double Amn = fAmn(m,n);
			
			for (int j = 0; j < fModuliDim; j++)
			{
				double tempj = Amn*fTensor2Table(n,j);
			
				for (int i = 0; i <= j; i++)
					moduli(i,j) += tempj*fTensor2Table(m,i);
			}
		}
	}
*/
	
	/* reversing the loops */
	for (int j = 0; j < fModuliDim; j++)
		for (int i = 0; i <= j; i++)
		{
			register double cij = 0.0;
			double* pm = fTensor2Table(0) + i;
			double* pn = fTensor2Table(0) + j;
		
			for (int n = 0; n < fNumBonds; n++)
			{
				double* pAmn = fAmn(n);
				double* pm2  = pm;
				double  nj   = *pn;
			
				for (int m = 0; m < fNumBonds; m++)
				{
					cij += nj*(*pAmn)*(*pm2);

					pAmn++;
					pm2 += fModuliDim;
				}
				
				pn += fModuliDim;
			}
			
			moduli(i,j) += cij;
		}

	moduli.CopySymmetric();
}
