/* $Id: EAM_particle.cpp,v 1.1.2.2 2004-02-25 17:24:46 paklein Exp $ */
/* created: hspark(02/25/2004) */
#include "EAM_particle.h"
#include <iostream.h> //TEMP
#include "CBLatticeT.h"
#include "C1FunctionT.h"

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
	fEDFunction(NULL)
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

	fEDFunction    = fEAMProperty->getElecDensEnergy();
	
	fEmbedEnergy   = fEAMProperty->getEmbedEnergy();
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

		double   r = fBonds[i];
		double phi = fPairEnergy(r, NULL, NULL);
		double rho = fDensity(r, NULL, NULL);

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
	double dFdrho = fEmbeddingEnergy->DFunction(rho);	

	/* assemble stress */
	stress = 0.0;
	
	dArrayT& DPotential = fPairPotential->MapDFunction(fBonds, fBond1);
	dArrayT& DDensity   = fElectronDensity->MapDFunction(fBonds, fBond2);
	for (int i = 0; i < fNumBonds; i++)
	{
		double ri = fBonds[i];
		int    ci = fCounts[i];		
		double coeff = (1.0/ri)*ci*(0.5*DPotential[i] + dFdrho*DDensity[i]);
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
	double dFdrho   = fEmbeddingEnergy->DFunction(rho);
	double d2Fdrho2 = fEmbeddingEnergy->DDFunction(rho);

	/* batched calls */
	dArrayT& DDensity    = fElectronDensity->MapDFunction(fBonds, fBond1);
	dArrayT& DDDensity   = fElectronDensity->MapDDFunction(fBonds, fBond2);
	dArrayT& DDPotential = fPairPotential->MapDDFunction(fBonds, fBond3);

	/* form upper triangle only */
	for (int j = 0; j < fNumBonds; j++)
	{
		double rj = fBonds[j];
		int    cj = fCounts[j];

		double DDPj = DDPotential[j];
		double DDDj = DDDensity[j];
		double DDj  = DDensity[j];
	
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
		
			/* mixed embedding energy term */
			Amn += ci*cj*d2Fdrho2*DDensity[i]*DDj;
		
			fAmn(i,j) = Amn/(ri*rj);
		}
	}
	
	fAmn.CopySymmetric();
}	

/*
* Compute the total electron density.
*/
double EAM_particle::TotalElectronDensity(void)
{
	/* compute total atomic density */
	dArrayT& ElectronDensity = fElectronDensity->MapFunction(fBonds, fBond1);

	double rho = 0.0;
	const int* pcount = fCounts.Pointer();
	double* pedensity = ElectronDensity.Pointer();

	for (int i = 0; i < fNumBonds; i++)
		rho += (*pcount++)*(*pedensity++);

	return rho;
}

/*
* Moduli tensor contributions.
*/
void EAM_particle::FormSingleBondContribution(double rho, dMatrixT& moduli)
{
	/* batch fetch */
	dArrayT& DPotential = fPairPotential->MapDFunction(fBonds, fBond1);
	dArrayT& DDensity   = fElectronDensity->MapDFunction(fBonds, fBond2);
	
	/* Embedding energy derivative */
	double dFdrho = fEmbeddingEnergy->DFunction(rho);	
	
	for (int i = 0; i < fNumBonds; i++)
	{
		double ri = fBonds[i];
	
		double coeff = -fCounts[i]*(0.5*DPotential[i] +
		                              dFdrho*DDensity[i])/(ri*ri*ri);
	
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
