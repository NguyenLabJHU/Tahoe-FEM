/* $Id: SmallStrainT.cpp,v 1.6.2.2 2002-09-22 23:08:57 paklein Exp $ */

#include "SmallStrainT.h"
#include "ShapeFunctionT.h"
#include "SSStructMatT.h"
#include "MaterialListT.h"

/* constructor */

using namespace Tahoe;

SmallStrainT::SmallStrainT(const ElementSupportT& support, const FieldT& field):
	ElasticT(support, field),
	fNeedsOffset(-1),
	fGradU(NumSD())
{

}

/* called immediately after constructor */
void SmallStrainT::Initialize(void)
{
	/* inherited */
	ElasticT::Initialize();

	/* what's needed */
	bool need_strain = false;
	bool need_strain_last = false;
	for (int i = 0; i < fMaterialNeeds.Length(); i++)
	{
		const ArrayT<bool>& needs = fMaterialNeeds[i];
		need_strain = need_strain || needs[fNeedsOffset + kstrain];
		need_strain_last = need_strain_last || needs[fNeedsOffset + kstrain_last];
	}

	/* allocate deformation gradient list */
	if (need_strain)
	{
		fStrain_List.Allocate(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fStrain_List[i].Allocate(NumSD());
	}
	
	/* allocate "last" deformation gradient list */
	if (need_strain_last)
	{
		fStrain_last_List.Allocate(NumIP());
		for (int i = 0; i < NumIP(); i++)
			fStrain_last_List[i].Allocate(NumSD());
	}
}

/***********************************************************************
* Protected
***********************************************************************/

/* construct list of materials from the input stream */
void SmallStrainT::ReadMaterialData(ifstreamT& in)
{
	/* inherited */
	ElasticT::ReadMaterialData(in);

	/* offset to class needs flags */
	fNeedsOffset = fMaterialNeeds[0].Length();
	
	/* set material needs */
	for (int i = 0; i < fMaterialNeeds.Length(); i++)
	{
		/* needs array */
		ArrayT<bool>& needs = fMaterialNeeds[i];

		/* resize array */
		needs.Resize(needs.Length() + 2, true);

		/* casts are safe since class contructs materials list */
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
		SSStructMatT* mat = (SSStructMatT*) pcont_mat;

		/* collect needs */
		needs[fNeedsOffset + kstrain     ] = mat->Need_Strain();
		needs[fNeedsOffset + kstrain_last] = mat->Need_Strain_last();
		
		/* consistency */
		needs[kNeedDisp] = needs[kNeedDisp] || needs[fNeedsOffset + kstrain];
		needs[KNeedLastDisp] = needs[KNeedLastDisp] || needs[fNeedsOffset + kstrain_last];
	}
}

/* initialize local field arrays. Allocate B-bar workspace if needed. */
void SmallStrainT::SetLocalArrays(void)
{
	/* inherited */
	ElasticT::SetLocalArrays();

	/* using B-bar */
	if (fStrainDispOpt == kMeanDilBbar) {
		fLocDispTranspose.Dimension(fLocDisp.Length());
		fMeanDilatation.Dimension(NumSD(), NumElementNodes());
	}
}

/* calculate the internal force contribution ("-k*d") */
void SmallStrainT::FormKd(double constK)
{
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();
	
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		/* strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanDilatation, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);

		/* B^T * Cauchy stress */
		fB.MultTx(fCurrMaterial->s_ij(), fNEEvec);
		
		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);
		
		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();
	}	
}

/* form the element stiffness matrix */
void SmallStrainT::FormStiffness(double constK)
{
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integrate element stiffness */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		double scale = constK*(*Det++)*(*Weight++);
	
		/* strain displacement matrix */
		if (fStrainDispOpt == kMeanDilBbar)
			Set_B_bar(fShapes->Derivatives_U(), fMeanDilatation, fB);
		else
			Set_B(fShapes->Derivatives_U(), fB);

		/* get D matrix */
		fD.SetToScaled(scale, fCurrMaterial->c_ijkl());
							
		/* multiply b(transpose) * db, taking account of symmetry, */
		/* and accumulate in elstif */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
}

/* compute the measures of strain/deformation over the element */
void SmallStrainT::SetGlobalShape(void)
{
	/* inherited */
	ElasticT::SetGlobalShape();

	/* material information */
	int material_number = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[material_number];
	
	/* using B-bar */
	if (fStrainDispOpt == kMeanDilBbar)
	{
		/* compute mean dilatation */
		SetMeanDilatation(fMeanDilatation);

		/* loop over integration points */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* set B-bar */
			int ip = fShapes->CurrIP();
			Set_B_bar(fShapes->Derivatives_U(ip), fMeanDilatation, fB);
	
			/* deformation gradient */
			if (needs[fNeedsOffset + kstrain])
			{
				/* transpose displacement array */
				fLocDisp.ReturnTranspose(fLocDispTranspose);

				/* compute strain using B-bar */
				dSymMatrixT& strain = fStrain_List[ip];
				fB.Multx(fLocDispTranspose, strain);
				strain.ScaleOffDiagonal(0.5);
			}

			/* "last" deformation gradient */
			if (needs[fNeedsOffset + kstrain_last])
			{
				/* transpose displacement array */
				fLocLastDisp.ReturnTranspose(fLocDispTranspose);

				/* compute strain using B-bar */
				dSymMatrixT& strain = fStrain_last_List[ip];
				fB.Multx(fLocDispTranspose, strain);
				strain.ScaleOffDiagonal(0.5);
			}
		}		
	}
	else
	{
		/* loop over integration points */
		for (int i = 0; i < NumIP(); i++)
		{
			/* deformation gradient */
			if (needs[fNeedsOffset + kstrain])
			{
				/* displacement gradient */
				fShapes->GradU(fLocDisp, fGradU, i);

				/* symmetric part */
				 fStrain_List[i].Symmetrize(fGradU);
			}

			/* "last" deformation gradient */
			if (needs[fNeedsOffset + kstrain_last])
			{
				/* displacement gradient */
				fShapes->GradU(fLocLastDisp, fGradU, i);

				/* symmetric part */
				 fStrain_last_List[i].Symmetrize(fGradU);
			}
		}
	}
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* compute mean dilatation, Hughes (4.5.23) */
void SmallStrainT::SetMeanDilatation(dArray2DT& mean_dilatation) const
{
	int nip = NumIP();
	const double* det = fShapes->IPDets();
	const double*   w = fShapes->IPWeights();

	/* volume */
	double vol = 0.0;
	for (int i = 0; i < nip; i++)
		vol += w[i]*det[i];

	/* initialize */
	mean_dilatation = 0.0;			

	/* integrate */
	for (int i = 0; i < nip; i++)
		mean_dilatation.AddScaled(w[i]*det[i]/vol, fShapes->Derivatives_U(i));
}

/* set B-bar as given by Hughes (4.5.11-16) */
void SmallStrainT::Set_B_bar(const dArray2DT& DNa, const dArray2DT& mean_dilatation, 
	dMatrixT& B)
{
#if __option(extended_errorcheck)
	if (B.Rows() != dSymMatrixT::NumValues(DNa.MajorDim()) ||
	    B.Cols() != DNa.Length() ||
	    mean_dilatation.MinorDim() != DNa.MinorDim() ||
	    mean_dilatation.MajorDim() != DNa.MajorDim())
	    throw eSizeMismatch;
#endif

	int nnd = DNa.MinorDim();
	double* pB = B.Pointer();

	/* 1D */
	if (DNa.MajorDim() == 1)
	{
		cout << "\n SmallStrainT::Set_B_bar: not implemented yet for 1D B-bar" << endl;
		throw eGeneralFail;
	}
	/* 2D */
	else if (DNa.MajorDim() == 2)
	{
		double* pNax = DNa(0);
		double* pNay = DNa(1);
			
		double* pBmx = mean_dilatation(0);
		double* pBmy = mean_dilatation(1);
			
		for (int i = 0; i < nnd; i++)
		{
			double factx = ((*pBmx++) - (*pNax))/3.0;
			double facty = ((*pBmy++) - (*pNay))/3.0;
			
			/* Hughes (4.5.11-16) */
			*pB++ = *pNax + factx;
			*pB++ = factx;
			*pB++ = *pNay;
	
			*pB++ = facty;
			*pB++ = *pNay + facty;
			*pB++ = *pNax;
				
			pNax++; pNay++;
		}
	}
	/* 3D */
	else		
	{
		double* pNax = DNa(0);
		double* pNay = DNa(1);
		double* pNaz = DNa(2);

		double* pBmx = mean_dilatation(0);
		double* pBmy = mean_dilatation(1);
		double* pBmz = mean_dilatation(2);
			
		for (int i = 0; i < nnd; i++)
		{
			double factx = ((*pBmx++) - (*pNax))/3.0;
			double facty = ((*pBmy++) - (*pNay))/3.0;
			double factz = ((*pBmz++) - (*pNaz))/3.0;

			/* Hughes (4.5.11-16) */
			*pB++ = *pNax + factx;
			*pB++ = factx;
			*pB++ = factx;
			*pB++ = 0.0;
			*pB++ = *pNaz;
			*pB++ = *pNay;

			*pB++ = facty;
			*pB++ = *pNay + facty;
			*pB++ = facty;
			*pB++ = *pNaz;
			*pB++ = 0.0;
			*pB++ = *pNax;
	
			*pB++ = factz;
			*pB++ = factz;
			*pB++ = *pNaz + factz;
			*pB++ = *pNay;
			*pB++ = *pNax;
			*pB++ = 0.0;
				
			pNax++; pNay++; pNaz++;
		}
	}
}
