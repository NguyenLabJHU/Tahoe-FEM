/* $Id: EAM_particle.cpp,v 1.2.8.4 2004-07-13 16:42:35 paklein Exp $ */
/* created: hspark(02/25/2004) */
#include "EAM_particle.h"
#include <iostream.h> //TEMP
#include "CBLatticeT.h"

/* EAM property types */
#include "ParadynEAMT.h"

using namespace Tahoe;

/* constructor */
EAM_particle::EAM_particle(CBLatticeT& lattice, const StringT& param_file): 
	fLattice(lattice),
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

/* Destructor */
EAM_particle::~EAM_particle(void) {
	delete fEAMProperty;
}

/* set "glue" functions and dimension work space */
void EAM_particle::Initialize(int nsd, int numbonds)
{
	/* dimension work space */
	int nstrs = dSymMatrixT::NumValues(nsd);
	fBondTensor4.Dimension(nstrs);
	fAmn.Dimension(numbonds);
	fBondTensor2.Dimension(nstrs);
	fTensor2Table.Dimension(numbonds, nstrs);
	fBond1.Dimension(numbonds);
	fBond2.Dimension(numbonds);
	fBond3.Dimension(numbonds);
}

/*
* Compute unit strain energy density:
*
*     unit strain energy = energy/atom
*/
double EAM_particle::ComputeUnitEnergy(void)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	double rho = 0.0;
	double energy = 0.0;
	
	const int* pcount = counts.Pointer();
	int nb = r.Length();
	for (int i = 0; i < nb; i++)
	{
		int ci = *pcount++;

		double  ri = r[i];
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
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	/* total atomic density */
	double rho = TotalElectronDensity();
	double dFdrho = fEmbedForce(rho, NULL, NULL);	// checks out

	/* assemble stress */
	stress = 0.0;
	int nb = r.Length();	
	for (int i = 0; i < nb; i++)	// fNumBonds, ri, ci check out
	{
		double ri = r[i];
		int    ci = counts[i];	
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

/*
* Compute the total electron density.
*/
double EAM_particle::TotalElectronDensity(void)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	double rho = 0.0;
	const int* pcount = counts.Pointer();
	int nb = r.Length();
	for (int i = 0; i < nb; i++)
	{
		double ri = r[i];
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
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	double dFdrho = fEmbedForce(rho, NULL, NULL);
	double d2Fdrho2 = fEmbedStiffness(rho, NULL, NULL);

	/* form upper triangle only */
	int nb = r.Length();
	for (int j = 0; j < nb; j++)
	{
		double rj = r[j];
		int    cj = counts[j];

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
			double ri = r[i];
			int    ci = counts[i];
					
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
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	/* Embedding energy derivative */
	double dFdrho = fEmbedForce(rho, NULL, NULL);
	int nb = r.Length();	
	for (int i = 0; i < nb; i++)
	{
		double ri = r[i];
		double a = fPairEnergy(ri, NULL, NULL);
		double Potential = a*a/ri;
		double DPotential = 2.0*fPairForce(ri, NULL, NULL)*a/ri-Potential/ri;
		double DDensity = fEDForce(ri, NULL, NULL);
		double coeff = -counts[i]*(0.5*DPotential + dFdrho*DDensity)/(ri*ri*ri);
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
	int dim = moduli.Rows();
	int nb = fLattice.NumberOfBonds();
	for (int j = 0; j < dim; j++)
		for (int i = 0; i <= j; i++)
		{
			register double cij = 0.0;
			double* pm = fTensor2Table(0) + i;
			double* pn = fTensor2Table(0) + j;
		
			for (int n = 0; n < nb; n++)
			{
				double* pAmn = fAmn(n);
				double* pm2  = pm;
				double  nj   = *pn;
			
				for (int m = 0; m < nb; m++)
				{
					cij += nj*(*pAmn)*(*pm2);

					pAmn++;
					pm2 += dim;
				}
				
				pn += dim;
			}
			
			moduli(i,j) += cij;
		}

	moduli.CopySymmetric();
}
