/* $Id: LengthsAndAnglesT.cpp,v 1.1.1.1 2001-01-29 08:20:26 paklein Exp $ */
/* created: paklein (05/26/1997)                                          */
/* Class to compute/manage all bond angles and derivatives                */
/* for the 2 unit cell, diamond cubic, modified Cauchy-Born,              */
/* constitutive equations.                                                */

#include "LengthsAndAnglesT.h"
#include <math.h>
#include "dArrayT.h"
#include "iArray2DT.h"

/* internal parameters */
const int	kNSD       = 3;
const int	kNumDOF	   = 3;
const int	kStressDim = dSymMatrixT::NumValues(kNSD);

/* Constructor */
LengthsAndAnglesT::LengthsAndAnglesT(const dMatrixT& Q, const iArray2DT& pairs):
	BondLengthsT(Q),
	fPairs(pairs),
	fR_3(fPairs.MajorDim(), kNSD),
	fl_3(fPairs.MajorDim()),
	fdl_3_dC(fPairs.MajorDim(), kNSD*kNSD),
	fd2l_3_dCdC(fPairs.MajorDim(), kStressDim*kStressDim),

//Cos(theta) = alpha/beta
	fCos12(fPairs.MajorDim()),

	fdCos_dXsi(fPairs.MajorDim(), kNSD),
	fd2Cos_dXsidXsi(fPairs.MajorDim(), kNSD*kNSD),
	fdCos_dC(fPairs.MajorDim(), kNSD*kNSD),
	fd2Cos_dCdC(fPairs.MajorDim(), kStressDim*kStressDim),
	fd2Cos_dCdXsi(fPairs.MajorDim(), kStressDim*kNumDOF),
	
//alpha
	fda_dXsi(kNSD),
	fd2a_dXsidXsi(kNSD),
	fdah_dC(kNSD),
	fd2ah_dCdC(kStressDim),
	fd2ah_dCdXsi(kStressDim,kNSD),

//beta
	fdb_dXsi(kNSD),
	fd2b_dXsidXsi(kNSD),
	fdbh_dC(kNSD),
	fd2bh_dCdC(kStressDim),
	fd2bh_dCdXsi(kStressDim,kNSD),
	
//workspace	
	fTempMat2(kNSD),
	fTempSymMat2(kNSD)
{
	Initialize();
}

/* Accessors */
const dArrayT& LengthsAndAnglesT::Cosines(void) const { return fCos12; }

/* Derivatives wrt Xsi */
const dArray2DT& LengthsAndAnglesT::dCos_dXsi(void) const { return fdCos_dXsi; }
const dArray2DT& LengthsAndAnglesT::d2Cos_dXsidXsi(void) const { return fd2Cos_dXsidXsi; }

/* Derivatives wrt. C */
const dArray2DT& LengthsAndAnglesT::dCos_hat_dC(void) const { return fdCos_dC; }
const dArray2DT& LengthsAndAnglesT::d2Cos_hat_dCdC(void) const { return fd2Cos_dCdC; }
const dArray2DT& LengthsAndAnglesT::d2Cos_hat_dCdXsi(void) const { return fd2Cos_dCdXsi; }

/* Set free dof - triggers recomputation */
void LengthsAndAnglesT::SetdXsi(const dMatrixT& CIJ, const dArrayT& Xsi)
{
	/* inherited */
	BondLengthsT::SetdXsi(CIJ, Xsi);

	dArrayT	 shR_3;
	dArrayT	 dl1dXsi, dl2dXsi;
	dMatrixT d2l1dXsidXsi, d2l2dXsidXsi;
	dArrayT	 dCosdXsi;
	dMatrixT d2CosdXsidXsi;

	for (int i = 0; i < fPairs.MajorDim(); i++)
	{
		int n1 = fPairs(i,0);
		int n2 = fPairs(i,1);
		
		shR_3.Set(kNSD,fR_3(i));

		CIJ.MultTx(shR_3,fTempVec);				
		
		/* deformed length of 3rd edge */
		double l3 = fl_3[i] = sqrt(dArrayT::Dot(shR_3,fTempVec));
		double l1 = fl[n1];		
		double l2 = fl[n2];
		
		/* angle factors */
		double alpha = l1*l1 + l2*l2 - l3*l3;
		double beta  = 2.0*l1*l2;
		
		fCos12[i] = alpha/beta;
		
		/* d/dXsi */
		dl1dXsi.Set(kNumDOF,fdl_dXsi(n1));
		dl2dXsi.Set(kNumDOF,fdl_dXsi(n2));
		
			/* alpha and beta */		
		fda_dXsi.SetToCombination(2.0*l1, dl1dXsi, 2.0*l2, dl2dXsi);
		fdb_dXsi.SetToCombination(2.0*l1, dl2dXsi, 2.0*l2, dl1dXsi);

			/* cos */
		dCosdXsi.Set(kNumDOF,fdCos_dXsi(i));
		dCosdXsi.SetToCombination(1.0/beta          , fda_dXsi,
		                          -alpha/(beta*beta), fdb_dXsi);

		/* d2/dXsidXsi */
		d2l1dXsidXsi.Set(kNumDOF,kNumDOF,fd2l_dXsidXsi(n1));
		d2l2dXsidXsi.Set(kNumDOF,kNumDOF,fd2l_dXsidXsi(n2));

			/* alpha */			
		fd2a_dXsidXsi.SetToCombination(2.0*l1,d2l1dXsidXsi,
		                               2.0*l2,d2l2dXsidXsi);
		
		fTempMat1.Outer(dl1dXsi,dl1dXsi);
		fTempMat2.Outer(dl2dXsi,dl2dXsi);

		fd2a_dXsidXsi.AddCombination(2.0, fTempMat1,
		                             2.0, fTempMat2);
			/* beta */
		fd2b_dXsidXsi.SetToCombination(2.0*l1,d2l2dXsidXsi,
		                               2.0*l2,d2l1dXsidXsi);

		fTempMat1.Outer(dl1dXsi,dl2dXsi);
		fTempMat2.Outer(dl2dXsi,dl1dXsi);

		fd2b_dXsidXsi.AddCombination(2.0, fTempMat1,
		                             2.0, fTempMat2);
		
			/* cos */
		d2CosdXsidXsi.Set(kNSD,kNSD,fd2Cos_dXsidXsi(i));

		d2CosdXsidXsi.Outer(fdb_dXsi,fdb_dXsi);
		d2CosdXsidXsi *= 2.0*alpha*pow(beta,-3.0);
		
		d2CosdXsidXsi.AddCombination(1.0/beta          , fd2a_dXsidXsi,
		                             -alpha/(beta*beta), fd2b_dXsidXsi);
		
		fTempMat1.Outer(fdb_dXsi,fda_dXsi);
		fTempMat2.Outer(fda_dXsi,fdb_dXsi);

		d2CosdXsidXsi.AddCombination(-pow(beta,-2.0), fTempMat1,
		                             -pow(beta,-2.0), fTempMat2);		
	}	
}

void LengthsAndAnglesT::SetdC(const dMatrixT& CIJ)
{
	/* inherited */
	BondLengthsT::SetdC(CIJ);

	/* Derivatives of l_3 wrt. C */
	dR(fR_3, fl_3, CIJ, fdl_3_dC, fd2l_3_dCdC);
	
	dMatrixT dl1dC, dl2dC, dl3dC;
	dMatrixT dCoshdC;

	for (int i = 0; i < fPairs.MajorDim(); i++)
	{
		int n1 = fPairs(i,0);
		int n2 = fPairs(i,1);
		
		/* deformed lengths */
		double l3 = fl_3[i];
		double l1 = fl[n1];		
		double l2 = fl[n2];
		
		/* angle factors */
		double alpha = l1*l1 + l2*l2 - l3*l3;
		double beta  = 2.0*l1*l2;
				
		/* d/dC */
		dl1dC.Set(kNSD,kNSD,fdl_dC(n1));
		dl2dC.Set(kNSD,kNSD,fdl_dC(n2));
		dl3dC.Set(kNSD,kNSD,fdl_3_dC(i));
		
			/* alpha and beta */
		fdah_dC.SetToScaled(-2.0*l3,dl3dC);			
		fdah_dC.AddCombination(2.0*l1, dl1dC, 2.0*l2, dl2dC);
		
		fdbh_dC.SetToCombination(2.0*l1, dl2dC, 2.0*l2, dl1dC);

			/* cos */
		dCoshdC.Set(kNSD,kNSD,fdCos_dC(i));
		dCoshdC.SetToCombination(1.0/beta          , fdah_dC,
		                         -alpha/(beta*beta), fdbh_dC);
	}
}

/* Set free dof - triggers recomputation */
void LengthsAndAnglesT::SetAll(const dMatrixT& CIJ)
{
	/* inherited */
	BondLengthsT::SetAll(CIJ);

	/* Derivatives of l_3 wrt. C */
	dR(fR_3, fl_3, CIJ, fdl_3_dC, fd2l_3_dCdC);
	
	dMatrixT dl1dC, dl2dC, dl3dC;
	dMatrixT d2l1dCdC, d2l2dCdC, d2l3dCdC;
	dMatrixT dCoshdC;
	dMatrixT d2CosdCdC, d2CosdCdXsi;
	dArrayT  dl1dXsi, dl2dXsi;
	dMatrixT d2l1dCdXsi, d2l2dCdXsi;

	for (int i = 0; i < fPairs.MajorDim(); i++)
	{
		int n1 = fPairs(i,0);
		int n2 = fPairs(i,1);
		
		/* deformed lengths */
		double l1 = fl[n1];		
		double l2 = fl[n2];
		double l3 = fl_3[i];
		
		/* angle factors */
		double alpha = l1*l1 + l2*l2 - l3*l3;
		double beta  = 2.0*l1*l2;
				
		/* d/dC */
		dl1dC.Set(kNSD,kNSD,fdl_dC(n1));
		dl2dC.Set(kNSD,kNSD,fdl_dC(n2));
		dl3dC.Set(kNSD,kNSD,fdl_3_dC(i));
		
			/* alpha and beta */
		fdah_dC.SetToScaled(-2.0*l3,dl3dC);			
		fdah_dC.AddCombination(2.0*l1, dl1dC, 2.0*l2, dl2dC);
		
		fdbh_dC.SetToCombination(2.0*l1, dl2dC, 2.0*l2, dl1dC);

			/* cos */
		dCoshdC.Set(kNSD,kNSD,fdCos_dC(i));
		dCoshdC.SetToCombination(1.0/beta          , fdah_dC,
		                         -alpha/(beta*beta), fdbh_dC);
		
		/* d2/dCdC */
		d2l1dCdC.Set(kStressDim, kStressDim, fd2l_dCdC(n1));
		d2l2dCdC.Set(kStressDim, kStressDim, fd2l_dCdC(n2));
		d2l3dCdC.Set(kStressDim, kStressDim, fd2l_3_dCdC(i));
		
			/* alpha */
		fd2ah_dCdC.SetToScaled(-2.0*l3, d2l3dCdC);
		fd2ah_dCdC.AddCombination(2.0*l1, d2l1dCdC, 2.0*l2, d2l2dCdC);

		fTempSymMat2.FromMatrix(dl3dC);
		fTempRank4.Outer(fTempSymMat2,fTempSymMat2);
		fd2ah_dCdC.AddScaled(-2.0,fTempRank4);		
			
		fTempSymMat1.FromMatrix(dl1dC);	
		fTempRank4.Outer(fTempSymMat1,fTempSymMat1);
		fd2ah_dCdC.AddScaled(2.0,fTempRank4);		

		fTempSymMat2.FromMatrix(dl2dC);
		fTempRank4.Outer(fTempSymMat2,fTempSymMat2);
		fd2ah_dCdC.AddScaled(2.0,fTempRank4);		

			/* beta */
		fd2bh_dCdC.SetToCombination(2.0*l2, d2l1dCdC, 2.0*l1, d2l2dCdC);

		fTempRank4.Outer(fTempSymMat2,fTempSymMat1);
		fd2bh_dCdC.AddScaled(2.0,fTempRank4);		

		fTempRank4.Outer(fTempSymMat1,fTempSymMat2);
		fd2bh_dCdC.AddScaled(2.0,fTempRank4);
		
			/* cos */
		d2CosdCdC.Set(kStressDim,kStressDim,fd2Cos_dCdC(i));	
			
		fTempSymMat1.FromMatrix(fdah_dC);
		fTempSymMat2.FromMatrix(fdbh_dC);

		d2CosdCdC.Outer(fTempSymMat2,fTempSymMat2);
		d2CosdCdC *= 2.0*alpha*pow(beta,-3.0);

		d2CosdCdC.AddCombination(1.0/beta          , fd2ah_dCdC,
		                         -alpha/(beta*beta), fd2bh_dCdC);
		
		fTempRank4.Outer(fTempSymMat1,fTempSymMat2);
		d2CosdCdC.AddScaled(-pow(beta,-2.0), fTempRank4);
		
		fTempRank4.Outer(fTempSymMat2,fTempSymMat1);
		d2CosdCdC.AddScaled(-pow(beta,-2.0), fTempRank4);
		
		/* d2/dCdXsi */
		fTempSymMat1.FromMatrix(dl1dC);
		fTempSymMat2.FromMatrix(dl2dC);

		dl1dXsi.Set(kNumDOF, fdl_dXsi(n1));
		dl2dXsi.Set(kNumDOF, fdl_dXsi(n2));

		d2l1dCdXsi.Set(kStressDim, kNumDOF, fd2l_dCdXsi(n1));
		d2l2dCdXsi.Set(kStressDim, kNumDOF, fd2l_dCdXsi(n2));

			/* alpha */
		fd2ah_dCdXsi.SetToCombination(2.0*l1,d2l1dCdXsi,2.0*l2,d2l2dCdXsi);	
	
		fTempMixed.Outer(fTempSymMat1,dl1dXsi);
		fd2ah_dCdXsi.AddScaled(2.0,fTempMixed);

		fTempMixed.Outer(fTempSymMat2,dl2dXsi);
		fd2ah_dCdXsi.AddScaled(2.0,fTempMixed);

			/* beta */					
		fd2bh_dCdXsi.SetToCombination(2.0*l2,d2l1dCdXsi,2.0*l1,d2l2dCdXsi);	
	
		fTempMixed.Outer(fTempSymMat2,dl1dXsi);
		fd2bh_dCdXsi.AddScaled(2.0,fTempMixed);

		fTempMixed.Outer(fTempSymMat1,dl2dXsi);
		fd2bh_dCdXsi.AddScaled(2.0,fTempMixed);
		
			/* cos */
		d2CosdCdXsi.Set(kStressDim, kNumDOF, fd2Cos_dCdXsi(i));

				/* need dbdXsi and dadXsi again */
		fda_dXsi.SetToCombination(2.0*l1, dl1dXsi, 2.0*l2, dl2dXsi);
		fdb_dXsi.SetToCombination(2.0*l1, dl2dXsi, 2.0*l2, dl1dXsi);

		d2CosdCdXsi.Outer(fTempSymMat2,fdb_dXsi);
		d2CosdCdXsi *= 2.0*alpha*pow(beta,-3.0);

		d2CosdCdXsi.AddCombination(1.0/beta          , fd2ah_dCdXsi,
		                           -alpha/(beta*beta), fd2bh_dCdXsi);
		
		fTempMixed.Outer(fTempSymMat1,fdb_dXsi);
		d2CosdCdXsi.AddScaled(-pow(beta,-2.0),fTempMixed);

		fTempMixed.Outer(fTempSymMat2,fda_dXsi);
		d2CosdCdXsi.AddScaled(-pow(beta,-2.0),fTempMixed);
	}	
}

/**********************************************************************
* Private
**********************************************************************/

/* called by constructor */
void LengthsAndAnglesT::Initialize(void)
{
	dArrayT R1, R2, R3;

	/* compute undeformed 3rd edge */
	for (int i = 0; i < fPairs.MajorDim(); i++)
	{
		/* fetch pointers */
		fR.RowAlias(fPairs(i,0), R1);
		fR.RowAlias(fPairs(i,1), R2);
		fR_3.RowAlias(i,R3);
		
		/* compute vector */
		R3.SetToCombination(1.0, R2,-1.0, R1);			
	}
}
