/* $Id: BondLengthsT.cpp,v 1.1.1.1 2001-01-29 08:20:26 paklein Exp $ */
/* created: paklein (05/20/1997)                                          */
/* Class to compute/manage all bond lengths and derivatives               */
/* for the 2 unit cell, diamond cubic, modified Cauchy-Born,              */
/* constitutive equations.                                                */

#include "BondLengthsT.h"
#include <math.h>
#include "dArrayT.h"

/* internal parameters */
const int	kNumAtoms  = 8;
const int   kNumBonds  = 8;
const int	kNSD       = 3;
const int	kStressDim = dSymMatrixT::NumValues(kNSD);

/* constructor */
BondLengthsT::BondLengthsT(const dMatrixT& Q):
	fR(kNumBonds,kNSD), fRmod(kNumBonds,kNSD), fl(kNumBonds),
	fdl_dXsi(kNumBonds,kNSD), fd2l_dXsidXsi(kNumBonds, kNSD*kNSD),
	fdl_dC(kNumBonds, kNSD*kNSD), fd2l_dCdC(kNumBonds, kStressDim*kStressDim),
	fd2l_dCdXsi(kNumBonds, kStressDim*kNSD),
	fTempVec(kNSD),
	fTempMat1(kNSD),
	fTempRank4(kStressDim),
	fTempSymMat1(kNSD),
	fTempMixed(kStressDim,kNSD)
{
	Initialize(Q);
}

/* Destructor */
BondLengthsT::~BondLengthsT(void) { }

/* set free dof - triggers recomputation */
void BondLengthsT::SetdXsi(const dMatrixT& CIJ, const dArrayT& Xsi)
{
	/* initialize fRmod */
	fRmod = fR;

	dArrayT	 shRmod, shl;
	dMatrixT shMat;

	for (int i = 0; i < fl.Length(); i++)
	{
		shRmod.Set(kNSD, fRmod(i));
		shRmod += Xsi;
		
		CIJ.MultTx(shRmod,fTempVec);				
		
		/* deformed lengths */
		fl[i] = sqrt(dArrayT::Dot(shRmod,fTempVec));
		
		/* dl/dXsi */
		shl.Set(kNSD, fdl_dXsi(i));
		shl.SetToScaled(1.0/fl[i], fTempVec);
		
		/* d2l/dXsi2 */
		shMat.Set(kNSD, kNSD, fd2l_dXsidXsi(i));
		shMat.SetToScaled(1.0/fl[i], CIJ);
		
		fTempMat1.Outer(fTempVec,fTempVec);
		shMat.AddScaled(-pow(fl[i],-3.0), fTempMat1);		
	}	
}

void BondLengthsT::SetdC(const dMatrixT& CIJ)
{
	/* Derivatives wrt. C */
	dR(fRmod, fl, CIJ, fdl_dC, fd2l_dCdC);
}

/* set free dof - triggers recomputation */
void BondLengthsT::SetAll(const dMatrixT& CIJ)
{
	/* Derivatives wrt. C  (and wrt. Xsi) */
	SetdC(CIJ);
	
	/* shallow workspace */
	dArrayT	 shRmod;
	dMatrixT shMat;
	
	/* mixed derivative */
	for (int i = 0; i < fl.Length(); i++)
	{
		shRmod.Set(kNSD, fRmod(i));		
		CIJ.MultTx(shRmod,fTempVec);
		
		fTempMat1.Outer(shRmod,shRmod);
		fTempSymMat1.FromMatrix(fTempMat1);	
						
		/* d2l/dCdXsi */
		shMat.Set(kStressDim, kNSD, fd2l_dCdXsi(i));
		
		fTempMixed.Outer(fTempSymMat1,fTempVec);
		shMat.SetToScaled(-0.5*pow(fl[i],-3.0),fTempMixed);
		
		Symmetrized(fTempMixed, shRmod);
		shMat.AddScaled(0.5/fl[i],fTempMixed);
	}	
}

/* Lengths */
const dArrayT& BondLengthsT::Lengths(void) const { return fl; }

/* Derivatives wrt. Xsi */
const dArray2DT& BondLengthsT::dl_dXsi(void) const { return fdl_dXsi; }
const dArray2DT& BondLengthsT::d2l_dXsidXsi(void) const { return fd2l_dXsidXsi; }

/* Derivatives wrt. C */
const dArray2DT& BondLengthsT::dl_hat_dC(void) const { return fdl_dC; }
const dArray2DT& BondLengthsT::d2l_hat_dCdC(void) const { return fd2l_dCdC; }
const dArray2DT& BondLengthsT::d2l_hat_dCdXsi(void) const { return fd2l_dCdXsi; }

/**********************************************************************
* Protected
**********************************************************************/

/* Compute dR/dC and d2R/dCdC */
void BondLengthsT::dR(const dArray2DT& R, const dArrayT& l, const dMatrixT& C,
	dArray2DT& dC, dArray2DT& dCdC)
{
	/* shallow workspace */
	dArrayT	 shR;
	dMatrixT shMat;

	for (int i = 0; i < l.Length(); i++)
	{
		shR.Set(kNSD, R(i));		
		C.MultTx(shR,fTempVec);				
						
		/* dl/dC */
		shMat.Set(kNSD, kNSD, dC(i));
		fTempMat1.Outer(shR,shR);

		shMat.SetToScaled(0.5/l[i],fTempMat1);

		/* d2l/dC2 */
		shMat.Set(kStressDim, kStressDim, dCdC(i));
		fTempSymMat1.FromMatrix(fTempMat1);
		fTempRank4.Outer(fTempSymMat1,fTempSymMat1);	
		
		shMat.SetToScaled(-0.25*pow(l[i],-3.0),fTempRank4);
	}	
}	

/**********************************************************************
* Private
**********************************************************************/

/* form symmetrized contribution to d2l_dCdXsi */
void BondLengthsT::Symmetrized(dMatrixT& mat, const dArrayT& Rmod)
{
	/* strictly 3D */
	if (Rmod.Length() != 3) throw eGeneralFail;

	int dex[6][2] =
		{{0,0},
		 {1,1},
		 {2,2},
		 {1,2},
		 {0,2},
		 {0,1}};

	for (int M = 0; M < 3; M++)
		for (int Rhat = 0; Rhat < 6; Rhat++)
		{
			int R = dex[Rhat][0];
			int S = dex[Rhat][1];
		
			mat(Rhat,M) = Rmod[R]*((M == S) ? 1 : 0) +
			              Rmod[S]*((M == R) ? 1 : 0);
		}

}

/* called by constructor */
void BondLengthsT::Initialize(const dMatrixT& Q)
{
	double sqrt3 = sqrt(3.0);
	double kLJ   = pow(2.0,1.0/6.0);

	/* undeformed positions of atoms in DC unit cell */
	double Xdata[kNumAtoms][kNSD] =
	   {{          0.0,           0.0,           0.0},
		{    kLJ/sqrt3,     kLJ/sqrt3,     kLJ/sqrt3},
		{2.0*kLJ/sqrt3, 2.0*kLJ/sqrt3,           0.0},
{   -kLJ/sqrt3,    -kLJ/sqrt3,     kLJ/sqrt3},
{2.0*kLJ/sqrt3,           0.0, 2.0*kLJ/sqrt3},
{   -kLJ/sqrt3,     kLJ/sqrt3,    -kLJ/sqrt3},
{          0.0, 2.0*kLJ/sqrt3, 2.0*kLJ/sqrt3},
{    kLJ/sqrt3,    -kLJ/sqrt3,    -kLJ/sqrt3}};

/* atoms in each pair - 2nd atom is free in each case */
	int pairdata[kNumBonds][2] =
		{{0,1},
		 {0,3},
		 {0,5},
		 {0,7},
		 {2,1},
		 {4,1},
		 {6,1},
		 {0,1}};

	/* compute undeformed bond vectors */
	dArrayT Xp, Xf, Ri;
	for (int i = 0; i < kNumBonds; i++)
	{
		Xp.Set( kNSD, Xdata[ pairdata[i][0] ] );
		Xf.Set( kNSD, Xdata[ pairdata[i][1] ] );
	
		fR.RowAlias(i,Ri);
		
		Ri  = Xf;
		Ri -= Xp;
	}
	
	/* transform bonds */
	if ( Q.Rows() > 0 )
	{
		dArrayT shbond, dpbond(kNSD);
	
		for (int i = 0; i < kNumBonds; i++)
		{
			/* get bond vector */
			fR.RowAlias(i, shbond);
			
			/* temp */
			dpbond = shbond;
		
			/* transform */
			Q.MultTx(dpbond, shbond);		
		}
	}
}
