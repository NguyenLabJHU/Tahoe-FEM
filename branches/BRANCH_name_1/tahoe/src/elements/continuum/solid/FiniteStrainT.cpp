/* $Id: FiniteStrainT.cpp,v 1.8.2.1 2002-06-27 18:02:46 cjkimme Exp $ */

#include "FiniteStrainT.h"
#include "ShapeFunctionT.h"
#include "FDStructMatT.h"
#include "MaterialList2DT.h"
#include "MaterialList3DT.h"

/* constructor */

using namespace Tahoe;

FiniteStrainT::FiniteStrainT(const ElementSupportT& support, const FieldT& field):
	ElasticT(support, field),
	fNeedsOffset(-1)
{

}

/* called immediately after constructor */
void FiniteStrainT::Initialize(void)
{
	/* inherited */
	ElasticT::Initialize();

	/* what's needed */
	bool need_F = false;
	bool need_F_last = false;
	for (int i = 0; i < fMaterialList->Length(); i++)
	{
		need_F = need_F || Needs_F(i);		
		need_F_last = need_F_last || Needs_F_last(i);
	}	

	/* allocate deformation gradient list */
	if (need_F)
	{
		int nip = NumIP();
		int nsd = NumSD();
		fF_all.Allocate(nip*nsd*nsd);
		fF_List.Allocate(nip);
		for (int i = 0; i < nip; i++)
			fF_List[i].Set(nsd, nsd, fF_all.Pointer(i*nsd*nsd));
	}
	
	/* allocate "last" deformation gradient list */
	if (need_F_last)
	{
		int nip = NumIP();
		int nsd = NumSD();
		fF_last_all.Allocate(nip*nsd*nsd);
		fF_last_List.Allocate(nip);
		for (int i = 0; i < nip; i++)
			fF_last_List[i].Set(nsd, nsd, fF_last_all.Pointer(i*nsd*nsd));
	}
}

/* compute field gradients with respect to current coordinates */
void FiniteStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const
{
	/* field gradient */
	fShapes->GradU(u, grad_u);
}

/* compute field gradients with respect to current coordinates */
void FiniteStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, 
	int ip) const
{
	/* field gradient */
	fShapes->GradU(u, grad_u, ip);
}

/* compute field gradients with respect to reference coordinates */
void FiniteStrainT::ComputeGradient_reference(const LocalArrayT& u, 
	dMatrixT& grad_u) const
{
#pragma unused(u)
#pragma unused(grad_u)
	cout << "\n FiniteStrainT::ComputeGradient_reference: not implemented" << endl;
	throw;
}

void FiniteStrainT::ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u, 
	int ip) const
{
#pragma unused(u)
#pragma unused(grad_u)
#pragma unused(ip)
	cout << "\n FiniteStrainT::ComputeGradient_reference: not implemented" << endl;
	throw;
}

/***********************************************************************
* Protected
***********************************************************************/

/* construct materials manager and read data */
MaterialListT* FiniteStrainT::NewMaterialList(int size) const
{
	if (NumSD() == 2)
		return new MaterialList2DT(size, *this);
	else if (NumSD() == 3)
		return new MaterialList3DT(size, *this);
	else
		return NULL;			
}

/* construct list of materials from the input stream */
void FiniteStrainT::ReadMaterialData(ifstreamT& in)
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
		FDStructMatT* mat = (FDStructMatT*) pcont_mat;

		/* collect needs */
		needs[fNeedsOffset + kF     ] = mat->Need_F();
		needs[fNeedsOffset + kF_last] = mat->Need_F_last();
		
		/* consistency */
		needs[kNeedDisp] = needs[kNeedDisp] || needs[fNeedsOffset + kF];
		needs[KNeedLastDisp] = needs[KNeedLastDisp] || needs[fNeedsOffset + kF_last];
	}
}

/* form shape functions and derivatives */
void FiniteStrainT::SetGlobalShape(void)
{
	/* inherited */
	ElasticT::SetGlobalShape();
	
	/* what needs to get computed */
	int material_number = CurrentElement().MaterialNumber();
	bool needs_F = Needs_F(material_number);
	bool needs_F_last = Needs_F_last(material_number);
	
	/* loop over integration points */
	for (int i = 0; i < NumIP(); i++)
	{
		/* deformation gradient */
		if (needs_F)
		{
			dMatrixT& mat = fF_List[i];

			/* displacement gradient */
			fShapes->GradU(fLocDisp, mat, i);

			/* add identity */
			mat.PlusIdentity();
		}

		/* "last" deformation gradient */
		if (needs_F_last)
		{
			dMatrixT& mat = fF_last_List[i];

			/* displacement gradient */
			fShapes->GradU(fLocLastDisp, mat, i);

			/* add identity */
			mat.PlusIdentity();
		}
	}
}

/* calculate the damping force contribution ("-c*v") */
void FiniteStrainT::FormCv(double constC)
{
#pragma unused(constC)
//NOTE: this "linear" damping was needed only to support
//      Rayleigh damping for linear problems. Will be removed 
//      soon.

#if 0
//approximate?
	/* clear workspace */
	fLHS = 0.0;
	fStressStiff = 0.0;

	/* form tangent stiffness */
	FormStiffness(constC);
	fLHS.CopySymmetric();

	/* reorder */
	fLocVel.ReturnTranspose(fTemp2);
	
	/* C*v */
	fLHS.MultTx(fTemp2, fNEEvec);
	
	/* Accumulate */
	fRHS += fNEEvec;
#endif
}

/* write all current element information to the stream */
void FiniteStrainT::CurrElementInfo(ostream& out) const
{
	/* inherited */
	ElasticT::CurrElementInfo(out);
	
	/* write deformation gradients */
	out << "\n i.p. deformation gradients:\n";
	for (int i = 0; i < fF_List.Length(); i++)
		out << " ip: " << i+1 << '\n'
		    << fF_List[i] << '\n';
	out << '\n';
}
