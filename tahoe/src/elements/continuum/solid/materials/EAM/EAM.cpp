/* $Id: EAM.cpp,v 1.17 2009-06-06 20:20:41 hspark Exp $ */
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
	fEdgeRho(9),
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
	fBondTensor2c.Dimension(nstrs);

	/* Interaction table for surface and edges */
	fIntType.Dimension(9,2);
	fIntType2.Dimension(44,2);
	
	/* Dimension edge cluster work space */
	/* Get deformed lengths for edge clusters */
	const dArrayT& rbe = fLattice.DeformedEdgeBulk();
	const dArrayT& re1 = fLattice.DeformedEdge1();
	const dArrayT& re2 = fLattice.DeformedEdge2();
	const dArrayT& re3 = fLattice.DeformedEdge3();
	const dArrayT& re4 = fLattice.DeformedEdge4();
	const dArrayT& re5 = fLattice.DeformedEdge5();
	const dArrayT& re6 = fLattice.DeformedEdge6();
	const dArrayT& re7 = fLattice.DeformedEdge7();
	const dArrayT& re8 = fLattice.DeformedEdge8();	
	fBond10.Dimension(rbe.Length());
	fBond11.Dimension(re1.Length());
	fBond12.Dimension(re2.Length());
	fBond13.Dimension(re3.Length());
	fBond14.Dimension(re4.Length());
	fBond15.Dimension(re5.Length());	
	fBond16.Dimension(re6.Length());
	fBond17.Dimension(re7.Length());
	fBond18.Dimension(re8.Length());	
	fBond19.Dimension(rbe.Length());
	fBond20.Dimension(re1.Length());
	fBond21.Dimension(re2.Length());
	fBond22.Dimension(re3.Length());
	fBond23.Dimension(re4.Length());
	fBond24.Dimension(re5.Length());	
	fBond25.Dimension(re6.Length());
	fBond26.Dimension(re7.Length());
	fBond27.Dimension(re8.Length());	
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
	energyb *= 2.0;	
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
* Compute unit surface edge energy density:
*
*     unit strain energy = energy/atom
*/
double EAM::ComputeUnitEdgeEnergy(void)
{
	/* Get deformed lengths for edge clusters */
	const dArrayT& rb = fLattice.DeformedEdgeBulk();
	const dArrayT& re1 = fLattice.DeformedEdge1();
	const dArrayT& re2 = fLattice.DeformedEdge2();
	const dArrayT& re3 = fLattice.DeformedEdge3();
	const dArrayT& re4 = fLattice.DeformedEdge4();
	const dArrayT& re5 = fLattice.DeformedEdge5();
	const dArrayT& re6 = fLattice.DeformedEdge6();
	const dArrayT& re7 = fLattice.DeformedEdge7();
	const dArrayT& re8 = fLattice.DeformedEdge8();	
	
	/* Get bond counts for each cluster type */
	const iArrayT& countsb = fLattice.EdgeBulkCounts();
	const iArrayT& countse1 = fLattice.Edge1Counts();
	const iArrayT& countse2 = fLattice.Edge2Counts();
	const iArrayT& countse3 = fLattice.Edge3Counts();
	const iArrayT& countse4 = fLattice.Edge4Counts();
	const iArrayT& countse5 = fLattice.Edge5Counts();
	const iArrayT& countse6 = fLattice.Edge6Counts();
	const iArrayT& countse7 = fLattice.Edge7Counts();
	const iArrayT& countse8 = fLattice.Edge8Counts();	
	
	/* values of electron density and pair potential */
	double rhob = 0.0;
	double rhoe1 = 0.0;
	double rhoe2 = 0.0;
	double rhoe3 = 0.0;
	double rhoe4 = 0.0;
	double rhoe5 = 0.0;
	double rhoe6 = 0.0;
	double rhoe7 = 0.0;
	double rhoe8 = 0.0;	
	double totalgamma = 0.0;
	double energyb = 0.0;
	double energye1 = 0.0;
	double energye2 = 0.0;
	double energye3 = 0.0;
	double energye4 = 0.0;
	double energye5 = 0.0;
	double energye6 = 0.0;
	double energye7 = 0.0;
	double energye8 = 0.0;	

	/* compute electron density at each neighbor for each atom type */
	dArrayT& ElecDensBulk = fElectronDensity->MapFunction(rb, fBond10);
	dArrayT& ElecDensEdge1 = fElectronDensity->MapFunction(re1, fBond11);
	dArrayT& ElecDensEdge2 = fElectronDensity->MapFunction(re2, fBond12);
	dArrayT& ElecDensEdge3 = fElectronDensity->MapFunction(re3, fBond13);
	dArrayT& ElecDensEdge4 = fElectronDensity->MapFunction(re4, fBond14);
	dArrayT& ElecDensEdge5 = fElectronDensity->MapFunction(re5, fBond15);
	dArrayT& ElecDensEdge6 = fElectronDensity->MapFunction(re6, fBond16);
	dArrayT& ElecDensEdge7 = fElectronDensity->MapFunction(re7, fBond17);
	dArrayT& ElecDensEdge8 = fElectronDensity->MapFunction(re8, fBond18);	
	dArrayT& PairBulk = fPairPotential->MapFunction(rb,fBond19);
	dArrayT& PairEdge1 = fPairPotential->MapFunction(re1,fBond20);
	dArrayT& PairEdge2 = fPairPotential->MapFunction(re2,fBond21);
	dArrayT& PairEdge3 = fPairPotential->MapFunction(re3,fBond22);
	dArrayT& PairEdge4 = fPairPotential->MapFunction(re4,fBond23);
	dArrayT& PairEdge5 = fPairPotential->MapFunction(re5,fBond24);
	dArrayT& PairEdge6 = fPairPotential->MapFunction(re6,fBond25);
	dArrayT& PairEdge7 = fPairPotential->MapFunction(re7,fBond26);
	dArrayT& PairEdge8 = fPairPotential->MapFunction(re8,fBond27);
	
	const int* pcountb = countsb.Pointer();
	const int* pcounte1 = countse1.Pointer();
	const int* pcounte2 = countse2.Pointer();
	const int* pcounte3 = countse3.Pointer();
	const int* pcounte4 = countse4.Pointer();
	const int* pcounte5 = countse5.Pointer();
	const int* pcounte6 = countse6.Pointer();
	const int* pcounte7 = countse7.Pointer();
	const int* pcounte8 = countse8.Pointer();	
	double* pedensb = ElecDensBulk.Pointer();
	double* pedense1 = ElecDensEdge1.Pointer();
	double* pedense2 = ElecDensEdge2.Pointer();
	double* pedense3 = ElecDensEdge3.Pointer();
	double* pedense4 = ElecDensEdge4.Pointer();
	double* pedense5 = ElecDensEdge5.Pointer();
	double* pedense6 = ElecDensEdge6.Pointer();
	double* pedense7 = ElecDensEdge7.Pointer();
	double* pedense8 = ElecDensEdge8.Pointer();	
	double* pphib = PairBulk.Pointer();
	double* pphie1 = PairEdge1.Pointer();
	double* pphie2 = PairEdge2.Pointer();
	double* pphie3 = PairEdge3.Pointer();
	double* pphie4 = PairEdge4.Pointer();
	double* pphie5 = PairEdge5.Pointer();
	double* pphie6 = PairEdge6.Pointer();
	double* pphie7 = PairEdge7.Pointer();
	double* pphie8 = PairEdge8.Pointer();	

	for (int i = 0; i < rb.Length(); i++)
	{
		int    cib = *pcountb++;
		rhob    += cib*(*pedensb++);
		energyb += cib*0.5*(*pphib++);
	}
	energyb += fEmbeddingEnergy->Function(rhob);

	for (int i = 0; i < re1.Length(); i++)
	{
		int    cie1 = *pcounte1++;
		rhoe1    += cie1*(*pedense1++);
		energye1 += cie1*0.5*(*pphie1++);
	}
	energye1 += fEmbeddingEnergy->Function(rhoe1);
		
	for (int i = 0; i < re2.Length(); i++)
	{
		int    cie2 = *pcounte2++;
		rhoe2    += cie2*(*pedense2++);
		energye2 += cie2*0.5*(*pphie2++);
	}
	energye2 += fEmbeddingEnergy->Function(rhoe2);
	for (int i = 0; i < re3.Length(); i++)
	{
		int    cie3 = *pcounte3++;
		rhoe3    += cie3*(*pedense3++);
		energye3 += cie3*0.5*(*pphie3++);
	}
	energye3 += fEmbeddingEnergy->Function(rhoe3);
		
	for (int i = 0; i < re4.Length(); i++)
	{
		int    cie4 = *pcounte4++;
		rhoe4    += cie4*(*pedense4++);
		energye4 += cie4*0.5*(*pphie4++);
	}
	energye4 += fEmbeddingEnergy->Function(rhoe4);
	for (int i = 0; i < re5.Length(); i++)
	{
		int    cie5 = *pcounte5++;
		rhoe5    += cie5*(*pedense5++);
		energye5 += cie5*0.5*(*pphie5++);
	}
	energye5 += fEmbeddingEnergy->Function(rhoe5);
		
	for (int i = 0; i < re6.Length(); i++)
	{
		int    cie6 = *pcounte6++;
		rhoe6    += cie6*(*pedense6++);
		energye6 += cie6*0.5*(*pphie6++);
	}
	energye6 += fEmbeddingEnergy->Function(rhoe6);
	for (int i = 0; i < re7.Length(); i++)
	{
		int    cie7 = *pcounte7++;
		rhoe7    += cie7*(*pedense7++);
		energye7 += cie7*0.5*(*pphie7++);
	}
	energye7 += fEmbeddingEnergy->Function(rhoe7);
		
	for (int i = 0; i < re8.Length(); i++)
	{
		int    cie8 = *pcounte8++;
		rhoe8    += cie8*(*pedense8++);
		energye8 += cie8*0.5*(*pphie8++);
	}
	energye8 += fEmbeddingEnergy->Function(rhoe8);	
	
	totalgamma += energyb;
	totalgamma += energye1;
	totalgamma += energye2;
	totalgamma += energye3;
	totalgamma += energye4;
	totalgamma += energye5;
	totalgamma += energye6;
	totalgamma += energye7;
	totalgamma += energye8;	
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
	const dArrayT& r = fLattice.DeformedSurfaceLengths();
	const iArrayT& counts = fLattice.BondCounts();
	const iArrayT& atom_type = fLattice.AtomTypes();	// added by HSP
	const dArray2DT& bonds = fLattice.Bonds();	// can access initial bonds
	
	/* calculate representative electron densities */
	ComputeElectronDensity();
	
	/* Create interaction table using calculated electron densities */
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
}

/*
* Compute unit 2nd PK edge stress - used by surface CB solver
*
*     unit 2nd PK edge stress = SIJ*(area per cell/atoms per cell)
*/
void EAM::ComputeUnitEdgeStress(dSymMatrixT& stress)
{
	const dArrayT& r = fLattice.DeformedEdgeLengths();
	const iArrayT& counts = fLattice.BondCounts();
	const iArrayT& edge_type = fLattice.EdgeTypes();	// added by HSP
	const dArray2DT& bonds = fLattice.Bonds();	// can access initial bonds

	/* calculate representative electron densities */
	ComputeEdgeElectronDensity();
	
	/* Create interaction table using calculated electron densities */
	// Will have 44 interactions for edges (fIntType2)
	// edge0 = bulk
	// edge1 = edge
	// edge2 = 1 atom away on surface 1 (-X normal)
	// edge3 = 2 atoms away on surface 1 (-X normal)
	// edge4 = 1 atom away on surface 1 (+Z normal)
	// edge5 = 2 atoms away on surface 1 (+Z normal)
	// edge6 = quasi surface 2
	// edge7 = real surface 2 (-X normal)
	// edge8 = real surface 2 (+Z normal)
	fIntType2(0,0) = (fEdgeRho[1]); // edge
	fIntType2(0,1) = (fEdgeRho[1]); // edge
	fIntType2(1,0) = (fEdgeRho[1]); // edge
	fIntType2(1,1) = (fEdgeRho[4]); // 1 atom away (+Z)
	fIntType2(2,0) = (fEdgeRho[1]); // edge
	fIntType2(2,1) = (fEdgeRho[6]); // quasi surface 2
	fIntType2(3,0) = (fEdgeRho[1]); // edge
	fIntType2(3,1) = (fEdgeRho[2]); // 1 atom away (-X)
	fIntType2(4,0) = (fEdgeRho[1]); // edge
	fIntType2(4,1) = (fEdgeRho[5]); // 2 atoms away (+Z)
	fIntType2(5,0) = (fEdgeRho[1]); // edge
	fIntType2(5,1) = (fEdgeRho[3]); // 2 atoms away (-X)
	fIntType2(6,0) = (fEdgeRho[1]); // edge
	fIntType2(6,1) = (fEdgeRho[8]); // real surface 2 (+Z)
	fIntType2(7,0) = (fEdgeRho[1]); // edge
	fIntType2(7,1) = (fEdgeRho[7]); // real surface 2 (-X)

	fIntType2(8,0) = (fEdgeRho[2]); // 1 atom away (-X)
	fIntType2(8,1) = (fEdgeRho[3]); // 2 atoms away (-X)
	fIntType2(9,0) = (fEdgeRho[2]); // 1 atom away (-X)
	fIntType2(9,1) = (fEdgeRho[2]); // 1 atom away (-X)
	fIntType2(10,0) = (fEdgeRho[2]); // 1 atom away (-X)
	fIntType2(10,1) = (fEdgeRho[8]); // real surface 2 (+Z)
	fIntType2(11,0) = (fEdgeRho[2]); // 1 atom away (-X)
	fIntType2(11,1) = (fEdgeRho[7]); // real surface 2 (-X)
	fIntType2(12,0) = (fEdgeRho[2]); // 1 atom away (-X)
	fIntType2(12,1) = (fEdgeRho[4]); // 1 atom away (+Z)
	fIntType2(13,0) = (fEdgeRho[2]); // 1 atom away (-X)
	fIntType2(13,1) = (fEdgeRho[6]); // quasi surface 2
	fIntType2(14,0) = (fEdgeRho[2]); // 1 atom away (-X)
	fIntType2(14,1) = (fEdgeRho[0]); // bulk
	fIntType2(15,0) = (fEdgeRho[2]); // 1 atom away (-X)
	fIntType2(15,1) = (fEdgeRho[5]); // 2 atoms away (+Z)

	fIntType2(16,0) = (fEdgeRho[3]); // 2 atoms away (-X)
	fIntType2(16,1) = (fEdgeRho[3]); // 2 atoms away (-X)
	fIntType2(17,0) = (fEdgeRho[3]); // 2 atoms away (-X)
	fIntType2(17,1) = (fEdgeRho[0]); // bulk
	fIntType2(18,0) = (fEdgeRho[3]); // 2 atoms away (-X)
	fIntType2(18,1) = (fEdgeRho[7]); // real surface 2 (-X)
	fIntType2(19,0) = (fEdgeRho[3]); // 2 atoms away (-X)
	fIntType2(19,1) = (fEdgeRho[6]); // quasi surface 2
	fIntType2(20,0) = (fEdgeRho[3]); // 2 atoms away (-X)
	fIntType2(20,1) = (fEdgeRho[2]); // 1 atom away (-X)
	fIntType2(21,0) = (fEdgeRho[3]); // 2 atoms away (-X)
	fIntType2(21,1) = (fEdgeRho[4]); // 1 atom away (+Z)
	fIntType2(22,0) = (fEdgeRho[3]); // 2 atoms away (-X)
	fIntType2(22,1) = (fEdgeRho[8]); // real surface 2 (+Z)

	fIntType2(23,0) = (fEdgeRho[4]); // 1 atom away (+Z)
	fIntType2(23,1) = (fEdgeRho[7]); // real surface 2 (-X)
	fIntType2(24,0) = (fEdgeRho[4]); // 1 atom away (+Z)
	fIntType2(24,1) = (fEdgeRho[4]); // 1 atom away (+Z)
	fIntType2(25,0) = (fEdgeRho[4]); // 1 atom away (+Z)
	fIntType2(25,1) = (fEdgeRho[5]); // 2 atoms away (+Z)
	fIntType2(26,0) = (fEdgeRho[4]); // 1 atom away (+Z)
	fIntType2(26,1) = (fEdgeRho[8]); // real surface 2 (+Z)
	fIntType2(27,0) = (fEdgeRho[4]); // 1 atom away (+Z)
	fIntType2(27,1) = (fEdgeRho[6]); // quasi surface 2
	fIntType2(28,0) = (fEdgeRho[4]); // 1 atom away (+Z)
	fIntType2(28,1) = (fEdgeRho[3]); // 2 atoms away (-X)
	fIntType2(29,0) = (fEdgeRho[4]); // 1 atom away (+Z)
	fIntType2(29,1) = (fEdgeRho[0]); // bulk

	fIntType2(30,0) = (fEdgeRho[5]); // 2 atoms away (+Z)
	fIntType2(30,1) = (fEdgeRho[0]); // bulk
	fIntType2(31,0) = (fEdgeRho[5]); // 2 atoms away (+Z)
	fIntType2(31,1) = (fEdgeRho[5]); // 2 atoms away (+Z)
	fIntType2(32,0) = (fEdgeRho[5]); // 2 atoms away (+Z)
	fIntType2(32,1) = (fEdgeRho[6]); // quasi surface 2
	fIntType2(33,0) = (fEdgeRho[5]); // 2 atoms away (+Z)
	fIntType2(33,1) = (fEdgeRho[8]); // real surface 2 (+Z)
	fIntType2(34,0) = (fEdgeRho[5]); // 2 atoms away (+Z)
	fIntType2(34,1) = (fEdgeRho[7]); // real surface 2 (-X)

	fIntType2(35,0) = (fEdgeRho[6]); // quasi surface 2
	fIntType2(35,1) = (fEdgeRho[7]); // real surface 2 (-X)
	fIntType2(36,0) = (fEdgeRho[6]); // quasi surface 2
	fIntType2(36,1) = (fEdgeRho[6]); // quasi surface 2
	fIntType2(37,0) = (fEdgeRho[6]); // quasi surface 2
	fIntType2(37,1) = (fEdgeRho[8]); // real surface 2 (+Z)
	fIntType2(38,0) = (fEdgeRho[6]); // quasi surface 2
	fIntType2(38,1) = (fEdgeRho[0]); // bulk 

	fIntType2(39,0) = (fEdgeRho[7]); // real surface 2 (-X)
	fIntType2(39,1) = (fEdgeRho[7]); // real surface 2 (-X)
	fIntType2(40,0) = (fEdgeRho[7]); // real surface 2 (-X)
	fIntType2(40,1) = (fEdgeRho[0]); // bulk
	fIntType2(41,0) = (fEdgeRho[7]); // real surface 2 (-X)
	fIntType2(41,1) = (fEdgeRho[8]); // real surface 2 (+Z)
	fIntType2(42,0) = (fEdgeRho[8]); // real surface 2 (+Z)
	fIntType2(42,1) = (fEdgeRho[0]); // bulk 
	fIntType2(43,0) = (fEdgeRho[8]); // real surface 2 (+Z)
	fIntType2(43,1) = (fEdgeRho[8]); // real surface 2 (+Z)	

	/* assemble stress */
	stress = 0.0;
	dArrayT& DPotential = fPairPotential->MapDFunction(r, fBond1);
	dArrayT& DDensity   = fElectronDensity->MapDFunction(r, fBond2);
	int nb = r.Length();
	dArrayT blah(2);

	for (int i = 0; i < nb; i++)
	{
		// define dFdrho_i and dFdrho_j here depending on rho_i and rho_j
		int type = edge_type[i];	// what is the interaction type?
		fIntType2.RowCopy(type, blah);	// get correct electron densities for i-j interactions
		double dFdrho_i = fEmbeddingEnergy->DFunction(blah[0]);
		double dFdrho_j = fEmbeddingEnergy->DFunction(blah[1]);
		double ri = r[i];
		int    ci = counts[i];		
		double coeff = (1.0/ri)*ci*0.5*(DPotential[i] + dFdrho_j*DDensity[i] + dFdrho_i*DDensity[i]);
		fLattice.BondComponentTensor2(i,fBondTensor2c);
		stress.AddScaled(coeff,fBondTensor2c);
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

/* calculate rho for edge CB */
void EAM::ComputeEdgeElectronDensity(void)
{
	/* Get deformed lengths for edge clusters */
	const dArrayT& rb = fLattice.DeformedEdgeBulk();
	const dArrayT& re1 = fLattice.DeformedEdge1();
	const dArrayT& re2 = fLattice.DeformedEdge2();
	const dArrayT& re3 = fLattice.DeformedEdge3();
	const dArrayT& re4 = fLattice.DeformedEdge4();
	const dArrayT& re5 = fLattice.DeformedEdge5();
	const dArrayT& re6 = fLattice.DeformedEdge6();
	const dArrayT& re7 = fLattice.DeformedEdge7();
	const dArrayT& re8 = fLattice.DeformedEdge8();	
	
	/* Get bond counts for each cluster type */
	const iArrayT& countsb = fLattice.EdgeBulkCounts();
	const iArrayT& countse1 = fLattice.Edge1Counts();
	const iArrayT& countse2 = fLattice.Edge2Counts();
	const iArrayT& countse3 = fLattice.Edge3Counts();
	const iArrayT& countse4 = fLattice.Edge4Counts();
	const iArrayT& countse5 = fLattice.Edge5Counts();
	const iArrayT& countse6 = fLattice.Edge6Counts();
	const iArrayT& countse7 = fLattice.Edge7Counts();
	const iArrayT& countse8 = fLattice.Edge8Counts();	
	
	/* values of electron density */
	double rhob = 0.0;
	double rhoe1 = 0.0;
	double rhoe2 = 0.0;
	double rhoe3 = 0.0;
	double rhoe4 = 0.0;
	double rhoe5 = 0.0;
	double rhoe6 = 0.0;
	double rhoe7 = 0.0;
	double rhoe8 = 0.0;	

	/* compute electron density at each neighbor for each atom type */
	dArrayT& ElecDensBulk = fElectronDensity->MapFunction(rb, fBond10);
	dArrayT& ElecDensEdge1 = fElectronDensity->MapFunction(re1, fBond11);
	dArrayT& ElecDensEdge2 = fElectronDensity->MapFunction(re2, fBond12);
	dArrayT& ElecDensEdge3 = fElectronDensity->MapFunction(re3, fBond13);
	dArrayT& ElecDensEdge4 = fElectronDensity->MapFunction(re4, fBond14);
	dArrayT& ElecDensEdge5 = fElectronDensity->MapFunction(re5, fBond15);
	dArrayT& ElecDensEdge6 = fElectronDensity->MapFunction(re6, fBond16);
	dArrayT& ElecDensEdge7 = fElectronDensity->MapFunction(re7, fBond17);
	dArrayT& ElecDensEdge8 = fElectronDensity->MapFunction(re8, fBond18);	

	const int* pcountb = countsb.Pointer();
	const int* pcounte1 = countse1.Pointer();
	const int* pcounte2 = countse2.Pointer();
	const int* pcounte3 = countse3.Pointer();
	const int* pcounte4 = countse4.Pointer();
	const int* pcounte5 = countse5.Pointer();
	const int* pcounte6 = countse6.Pointer();
	const int* pcounte7 = countse7.Pointer();
	const int* pcounte8 = countse8.Pointer();	
	double* pedensb = ElecDensBulk.Pointer();
	double* pedense1 = ElecDensEdge1.Pointer();
	double* pedense2 = ElecDensEdge2.Pointer();
	double* pedense3 = ElecDensEdge3.Pointer();
	double* pedense4 = ElecDensEdge4.Pointer();
	double* pedense5 = ElecDensEdge5.Pointer();
	double* pedense6 = ElecDensEdge6.Pointer();
	double* pedense7 = ElecDensEdge7.Pointer();
	double* pedense8 = ElecDensEdge8.Pointer();	

	/* bulk atoms */
	for (int i = 0; i < rb.Length(); i++)
		rhob += (*pcountb++)*(*pedensb++);

	/* edge1 atoms */
	for (int i = 0; i < re1.Length(); i++)
		rhoe1 += (*pcounte1++)*(*pedense1++);	
	
	/* edge2 atoms */
	for (int i = 0; i < re2.Length(); i++)
		rhoe2 += (*pcounte2++)*(*pedense2++);
	
	/* edge3 atoms */
	for (int i = 0; i < re3.Length(); i++)
		rhoe3 += (*pcounte3++)*(*pedense3++);	
	
	/* edge4 atoms */
	for (int i = 0; i < re4.Length(); i++)
		rhoe4 += (*pcounte4++)*(*pedense4++);
		
	/* edge5 atoms */
	for (int i = 0; i < re5.Length(); i++)
		rhoe5 += (*pcounte5++)*(*pedense5++);	
	
	/* edge6 atoms */
	for (int i = 0; i < re6.Length(); i++)
		rhoe6 += (*pcounte6++)*(*pedense6++);
		
	/* edge7 atoms */
	for (int i = 0; i < re7.Length(); i++)
		rhoe7 += (*pcounte7++)*(*pedense7++);	
	
	/* edge8 atoms */
	for (int i = 0; i < re8.Length(); i++)
		rhoe8 += (*pcounte8++)*(*pedense8++);		
	
	double blah[9] = {rhob,rhoe1,rhoe2,rhoe3,rhoe4,rhoe5,rhoe6,rhoe7,rhoe8};
	fEdgeRho = blah;
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
