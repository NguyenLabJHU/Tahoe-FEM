/* $Id: EAM.cpp,v 1.8 2006-07-03 01:38:28 hspark Exp $ */
/* created: paklein (12/02/1996) */
#include "EAM.h"
#include "CBLatticeT.h"
#include "C1FunctionT.h"

using namespace Tahoe;

/* constructor */
EAM::EAM(CBLatticeT& lattice):
	fPairPotential(NULL),
	fEmbeddingEnergy(NULL),
	fElectronDensity(NULL),
	fLattice(lattice)
{

}

/* Destructor */
EAM::~EAM(void)
{
	delete fPairPotential;
	delete fEmbeddingEnergy;
	delete fElectronDensity;
}

/* Set "glue" functions and dimension work space */
void EAM::Initialize(int nsd, int numbonds)
{
	/* glue functions */
	SetPairPotential();
	SetEmbeddingEnergy();
	SetElectronDensity(); 		

	/* dimension work space */
	int nstrs = dSymMatrixT::NumValues(nsd);
	fBondTensor4.Dimension(nstrs),
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
double EAM::ComputeUnitEnergy(void)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	/* compute total atomic density */	
	dArrayT& ElectronDensity = fElectronDensity->MapFunction(r, fBond1);
	dArrayT& PairPotential   = fPairPotential->MapFunction(r, fBond2);
	
	double rho = 0.0;
	double energy = 0.0;
	
	const int* pcount = counts.Pointer();
	double* prho = ElectronDensity.Pointer();
	double* pphi = PairPotential.Pointer();
	int nb = r.Length();
	for (int i = 0; i < nb; i++)
	{
		int    ci = *pcount++;
		
		rho    += ci*(*prho++);
		energy += ci*0.5*(*pphi++);
	}
	
	energy += fEmbeddingEnergy->Function(rho);
	return energy;
}

/*
* Compute unit 2nd PK stress:
*
*     unit 2nd PK stress = SIJ*(volume per cell/atoms per cell)
*/
void EAM::ComputeUnitStress(dSymMatrixT& stress)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	/* total atomic density */
	dArrayT rho = ElectronDensityAtNeighbors();
	double rho0 = rho[0];
	double dFdrho_i = fEmbeddingEnergy->DFunction(rho[0]);

	/* assemble stress */
	stress = 0.0;
	dArrayT& DPotential = fPairPotential->MapDFunction(r, fBond1);
	dArrayT& DDensity   = fElectronDensity->MapDFunction(r, fBond2);
	int nb = r.Length();

	for (int i = 0; i < nb; i++)
	{
		double dFdrho_j = fEmbeddingEnergy->DFunction(rho[i+1]);
		double ri = r[i];
		int    ci = counts[i];		
		double coeff = (1.0/ri)*ci*0.5*(DPotential[i] + dFdrho_j*DDensity[i] + dFdrho_i*DDensity[i]);
		//double coeff = (1.0/ri)*ci*(0.5*DPotential[i] + dFdrho_j*DDensity[i]); // old expression
		fLattice.BondComponentTensor2(i,fBondTensor2);
		stress.AddScaled(coeff,fBondTensor2);
	}
}
	   	    	
/*
* Compute unit material tangent moduli:
*
*   unit material tangent moduli = CIJKL*(volume per cell/atoms per cell)
*/
void EAM::ComputeUnitModuli(dMatrixT& moduli)
{
	/* compute total electron density */
	dArrayT rho = ElectronDensityAtNeighbors();

	/* initialize */
	moduli = 0.0;

	/* mixed bond contribution */
	FormMixedBondContribution(rho, moduli);
	
	/* single bond contribution */
	FormSingleBondContribution(rho, moduli);
}

/* compute the total electron density at atom_i and its neighbors atom_j */
dArrayT EAM::ElectronDensityAtNeighbors(void)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();
	dArrayT rhoall(counts.Length()+1);	// rho_i is the first value, rho_j are others

	/* compute total atomic density */
	dArrayT& ElectronDensity = fElectronDensity->MapFunction(r, fBond1);

	double rho = 0.0;
	const int* pcount = counts.Pointer();
	double* pedensity = ElectronDensity.Pointer();
	int nb = r.Length();

	for (int i = 0; i < nb; i++)
		rho += (*pcount++)*(*pedensity++);

	rhoall = rho;
	return rhoall;
}

/* compute the total electron density */
double EAM::TotalElectronDensity(void)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	/* compute total atomic density */
	dArrayT& ElectronDensity = fElectronDensity->MapFunction(r, fBond1);

	double rho = 0.0;
	const int* pcount = counts.Pointer();
	double* pedensity = ElectronDensity.Pointer();
	int nb = r.Length();

	for (int i = 0; i < nb; i++)
		rho += (*pcount++)*(*pedensity++);

	return rho;
}

/**********************************************************************
 * Private
 **********************************************************************/

/*
* Form matrix of mixed pair potential and embedding
* energy derivatives.  NOTE: computes the UPPER triangle
* ONLY.
*/
void EAM::FormMixedDerivatives(dArrayT rho)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	//double dFdrho   = fEmbeddingEnergy->DFunction(rho);
	//double d2Fdrho2 = fEmbeddingEnergy->DDFunction(rho);

	/* batched calls */
	dArrayT& DDensity    = fElectronDensity->MapDFunction(r, fBond1);
	dArrayT& DDDensity   = fElectronDensity->MapDDFunction(r, fBond2);
	dArrayT& DDPotential = fPairPotential->MapDDFunction(r, fBond3);

	/* form upper triangle only */
	int nb = r.Length();
	for (int j = 0; j < nb; j++)
	{
		double rj = r[j];
		int    cj = counts[j];

		double DDPj = DDPotential[j];
		double DDDj = DDDensity[j];
		double DDj  = DDensity[j];
	
		/* Embedding energy derivative */
		double dFdrho = fEmbeddingEnergy->DFunction(rho[j]);	
		double d2Fdrho2 = fEmbeddingEnergy->DDFunction(rho[j]);
		
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
		
			/* mixed embedding energy term */
			Amn += ci*cj*d2Fdrho2*DDensity[i]*DDj;
		
			fAmn(i,j) = Amn/(ri*rj);
		}
	}
	
	fAmn.CopySymmetric();
}	

/* moduli tensor contributions */
void EAM::FormSingleBondContribution(dArrayT rho, dMatrixT& moduli)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	/* batch fetch */
	dArrayT& DPotential = fPairPotential->MapDFunction(r, fBond1);
	dArrayT& DDensity   = fElectronDensity->MapDFunction(r, fBond2);

	/* Embedding energy derivative */
	//double dFdrho = fEmbeddingEnergy->DFunction(rho[i]);	
	
	int nb = r.Length();
	for (int i = 0; i < nb; i++)
	{
		/* Embedding energy derivative */
		double dFdrho = fEmbeddingEnergy->DFunction(rho[i]);	

		double ri = r[i];
	
		double coeff = -counts[i]*(0.5*DPotential[i] + dFdrho*DDensity[i])/(ri*ri*ri);
	
		fLattice.BondComponentTensor4(i,fBondTensor4);		
		moduli.AddScaled(coeff,fBondTensor4);
	}
}

void EAM::FormMixedBondContribution(dArrayT rho, dMatrixT& moduli)
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
