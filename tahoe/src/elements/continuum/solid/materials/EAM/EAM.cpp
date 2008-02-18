/* $Id: EAM.cpp,v 1.16 2008-02-18 22:52:59 hspark Exp $ */
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
	fRepRho(3),
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
	
	/* Dimension surface cluster work space */
	const dArrayT& rb = fLattice.DeformedBulk();
	const dArrayT& rs1 = fLattice.DeformedSurf1();
	const dArrayT& rs2 = fLattice.DeformedSurf2();
	fBond4.Dimension(rb.Length());
	fBond5.Dimension(rs1.Length());
	fBond6.Dimension(rs2.Length());
	fBond7.Dimension(rb.Length());
	fBond8.Dimension(rs1.Length());
	fBond9.Dimension(rs2.Length());
	fBondTensor2b.Dimension(nstrs);
	//fIntType.Dimension(6,2);
	fIntType.Dimension(9,2);
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
* Compute unit surface strain energy density:
*
*     unit strain energy = energy/atom
*/
double EAM::ComputeUnitSurfaceEnergy(void)
{
	/* Get deformed lengths for bulk, surface1 and surface2 clusters */
	const dArrayT& rb = fLattice.DeformedBulk();
	const dArrayT& rs1 = fLattice.DeformedSurf1();
	const dArrayT& rs2 = fLattice.DeformedSurf2();
	
	/* Get bond counts for each cluster type */
	const iArrayT& countsb = fLattice.BulkCounts();
	const iArrayT& countss1 = fLattice.Surf1Counts();
	const iArrayT& countss2 = fLattice.Surf2Counts();
	
	/* values of electron density, pair potential */
	double rhob = 0.0;
	double rhos1 = 0.0;
	double rhos2 = 0.0;
	double energyb = 0.0;
	double energys1 = 0.0;
	double energys2 = 0.0;
	double totalgamma = 0.0;

	/* compute electron density and pair potential at each neighbor for each atom type */
	dArrayT& ElecDensBulk = fElectronDensity->MapFunction(rb, fBond4);
	dArrayT& ElecDensSurf1 = fElectronDensity->MapFunction(rs1, fBond5);
	dArrayT& ElecDensSurf2 = fElectronDensity->MapFunction(rs2, fBond6);
	dArrayT& PairBulk = fPairPotential->MapFunction(rb, fBond7);
	dArrayT& PairSurf1 = fPairPotential->MapFunction(rs1, fBond8);
	dArrayT& PairSurf2 = fPairPotential->MapFunction(rs2, fBond9);

	const int* pcountb = countsb.Pointer();
	const int* pcounts1 = countss1.Pointer();
	const int* pcounts2 = countss2.Pointer();
	double* pedensb = ElecDensBulk.Pointer();
	double* pedenss1 = ElecDensSurf1.Pointer();
	double* pedenss2 = ElecDensSurf2.Pointer();
	double* pphib = PairBulk.Pointer();
	double* pphis1 = PairSurf1.Pointer();
	double* pphis2 = PairSurf2.Pointer();

	for (int i = 0; i < rb.Length(); i++)
	{
		int    cib = *pcountb++;
		rhob    += cib*(*pedensb++);
		energyb += cib*0.5*(*pphib++);
	}
	energyb += fEmbeddingEnergy->Function(rhob);
	energyb *= 1.5;	
	//1.5 due to splitting of layer 4 bulk energy to balance # atoms subtracted

	for (int i = 0; i < rs1.Length(); i++)
	{
		int    cis1 = *pcounts1++;
		rhos1    += cis1*(*pedenss1++);
		energys1 += cis1*0.5*(*pphis1++);
	}
	energys1 += fEmbeddingEnergy->Function(rhos1);
		
	for (int i = 0; i < rs2.Length(); i++)
	{
		int    cis2 = *pcounts2++;
		rhos2    += cis2*(*pedenss2++);
		energys2 += cis2*0.5*(*pphis2++);
	}
	energys2 += fEmbeddingEnergy->Function(rhos2);

	totalgamma += energyb;
	totalgamma += energys1;
	totalgamma += energys2;

	return totalgamma;	// overprediction of atoms/area as compared to reality
}


/*
* Compute unit 2nd PK stress - used by bulk CB solver
*
*     unit 2nd PK stress = SIJ*(volume per cell/atoms per cell)
*/
void EAM::ComputeUnitStress(dSymMatrixT& stress)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	/* total atomic density */
	double rho = TotalElectronDensity();
	double dFdrho = fEmbeddingEnergy->DFunction(rho);

	/* assemble stress */
	stress = 0.0;
	dArrayT& DPotential = fPairPotential->MapDFunction(r, fBond1);
	dArrayT& DDensity   = fElectronDensity->MapDFunction(r, fBond2);
	int nb = r.Length();
	
	for (int i = 0; i < nb; i++)
	{
		double ri = r[i];
		int    ci = counts[i];		
		double coeff = (1.0/ri)*ci*(0.5*DPotential[i] + dFdrho*DDensity[i]); // old expression
		//double coeff = (1.0/ri)*ci*0.5*(DPotential[i] + dFdrho*DDensity[i] + dFdrho*DDensity[i]);
		fLattice.BondComponentTensor2(i,fBondTensor2);
		stress.AddScaled(coeff,fBondTensor2);
	}
	double asdf = ComputeUnitEnergy();
}

/*
* Compute unit 2nd PK surface stress - used by surface CB solver
*
*     unit 2nd PK surface stress = SIJ*(area per cell/atoms per cell)
*/
void EAM::ComputeUnitSurfaceStress(dSymMatrixT& stress)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();
	const iArrayT& atom_type = fLattice.AtomTypes();	// added by HSP
	const dArray2DT& bonds = fLattice.Bonds();	// can access initial bonds
	
	/* calculate representative electron densities */
	ComputeElectronDensity();
	
	/* Create interaction table using calculated electron densities */
	// ADD THREE NEW INTERACTION TYPES HERE
	fIntType(0,0) = (fRepRho[1]);	// surface1
	fIntType(0,1) = (fRepRho[1]);	// surface1
	fIntType(1,0) = (fRepRho[1]);	// surface1
	fIntType(1,1) = (fRepRho[2]);	// surface2
	fIntType(2,0) = (fRepRho[1]);	// surface1
	fIntType(2,1) = (fRepRho[0]);	// bulk
	fIntType(3,0) = (fRepRho[2]);	// surface2
	fIntType(3,1) = (fRepRho[1]);	// surface1
	fIntType(4,0) = (fRepRho[2]);	// surface2
	fIntType(4,1) = (fRepRho[2]);	// surface2
	fIntType(5,0) = (fRepRho[2]);	// surface2
	fIntType(5,1) = (fRepRho[0]);	// bulk
	fIntType(6,0) = (fRepRho[0]);	// bulk	- new additions for S3 and S4 layers
	fIntType(6,1) = (fRepRho[1]);	// surface1
	fIntType(7,0) = (fRepRho[0]);	// bulk
	fIntType(7,1) = (fRepRho[2]);	// surface2
	fIntType(8,0) = (fRepRho[0]);	// bulk
	fIntType(8,1) = (fRepRho[0]);	// bulk

	/* assemble stress */
	stress = 0.0;
	dArrayT& DPotential = fPairPotential->MapDFunction(r, fBond1);
	dArrayT& DDensity   = fElectronDensity->MapDFunction(r, fBond2);
	int nb = r.Length();
	dArrayT blah(2);

	for (int i = 0; i < nb; i++)
	{
		// define dFdrho_i and dFdrho_j here depending on rho_i and rho_j
		int type = atom_type[i];	// what is the interaction type?
		fIntType.RowCopy(type, blah);	// get correct electron densities for i-j interactions
		double dFdrho_i = fEmbeddingEnergy->DFunction(blah[0]);
		double dFdrho_j = fEmbeddingEnergy->DFunction(blah[1]);
		double ri = r[i];
		int    ci = counts[i];		
		double coeff = (1.0/ri)*ci*0.5*(DPotential[i] + dFdrho_j*DDensity[i] + dFdrho_i*DDensity[i]);
		fLattice.BondComponentTensor2(i,fBondTensor2b);
		stress.AddScaled(coeff,fBondTensor2b);
	}
	double asdf = ComputeUnitSurfaceEnergy();
}
	   	    	
/*
* Compute unit material tangent moduli:
*
*   unit material tangent moduli = CIJKL*(volume per cell/atoms per cell)
*/
void EAM::ComputeUnitModuli(dMatrixT& moduli)
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

/* calculate rho for bulk, surface1 and surface2 atoms */
void EAM::ComputeElectronDensity(void)
{
	/* Get deformed lengths for bulk, surface1 and surface2 clusters */
	const dArrayT& rb = fLattice.DeformedBulk();
	const dArrayT& rs1 = fLattice.DeformedSurf1();
	const dArrayT& rs2 = fLattice.DeformedSurf2();
	
	/* Get bond counts for each cluster type */
	const iArrayT& countsb = fLattice.BulkCounts();
	const iArrayT& countss1 = fLattice.Surf1Counts();
	const iArrayT& countss2 = fLattice.Surf2Counts();
	
	/* values of electron density */
	double rhob = 0.0;
	double rhos1 = 0.0;
	double rhos2 = 0.0;

	/* compute electron density at each neighbor for each atom type */
	dArrayT& ElecDensBulk = fElectronDensity->MapFunction(rb, fBond4);
	dArrayT& ElecDensSurf1 = fElectronDensity->MapFunction(rs1, fBond5);
	dArrayT& ElecDensSurf2 = fElectronDensity->MapFunction(rs2, fBond6);

	const int* pcountb = countsb.Pointer();
	const int* pcounts1 = countss1.Pointer();
	const int* pcounts2 = countss2.Pointer();
	double* pedensb = ElecDensBulk.Pointer();
	double* pedenss1 = ElecDensSurf1.Pointer();
	double* pedenss2 = ElecDensSurf2.Pointer();
	int nbb = rb.Length();
	int nbs1 = rs1.Length();
	int nbs2 = rs2.Length();

	/* bulk atoms */
	for (int i = 0; i < nbb; i++)
		rhob += (*pcountb++)*(*pedensb++);

	/* surface1 atoms */
	for (int i = 0; i < nbs1; i++)
		rhos1 += (*pcounts1++)*(*pedenss1++);	
	
	/* surface2 atoms */
	for (int i = 0; i < nbs2; i++)
		rhos2 += (*pcounts2++)*(*pedenss2++);
	
	double blah[3] = {rhob,rhos1,rhos2};
	fRepRho = blah;
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
void EAM::FormMixedDerivatives(double rho)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	double dFdrho   = fEmbeddingEnergy->DFunction(rho);
	double d2Fdrho2 = fEmbeddingEnergy->DDFunction(rho);

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
		//double dFdrho = fEmbeddingEnergy->DFunction(rho[j]);	
		//double d2Fdrho2 = fEmbeddingEnergy->DDFunction(rho[j]);
		
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
void EAM::FormSingleBondContribution(double rho, dMatrixT& moduli)
{
	const dArrayT& r = fLattice.DeformedLengths();
	const iArrayT& counts = fLattice.BondCounts();

	/* batch fetch */
	dArrayT& DPotential = fPairPotential->MapDFunction(r, fBond1);
	dArrayT& DDensity   = fElectronDensity->MapDFunction(r, fBond2);

	/* Embedding energy derivative */
	double dFdrho = fEmbeddingEnergy->DFunction(rho);	
	
	int nb = r.Length();
	for (int i = 0; i < nb; i++)
	{
		/* Embedding energy derivative */
		//double dFdrho = fEmbeddingEnergy->DFunction(rho[i]);	

		double ri = r[i];
	
		double coeff = -counts[i]*(0.5*DPotential[i] + dFdrho*DDensity[i])/(ri*ri*ri);
	
		fLattice.BondComponentTensor4(i,fBondTensor4);		
		moduli.AddScaled(coeff,fBondTensor4);
	}
}

void EAM::FormMixedBondContribution(double rho, dMatrixT& moduli)
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
