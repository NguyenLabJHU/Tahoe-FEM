/* $Id: FiniteStrainT.cpp,v 1.19 2003-12-28 08:23:20 paklein Exp $ */
#include "FiniteStrainT.h"

#include "ShapeFunctionT.h"
#include "FSSolidMatT.h"
#include "FSMatSupportT.h"

/* materials lists */
#include "SolidMatList1DT.h"
#include "SolidMatList2DT.h"
#include "SolidMatList3DT.h"

using namespace Tahoe;

/* constructor */
FiniteStrainT::FiniteStrainT(const ElementSupportT& support, const FieldT& field):
	SolidElementT(support, field),
	fNeedsOffset(-1),
	fCurrShapes(NULL),
	fFSMatSupport(NULL)
{
	/* disable any strain-displacement options */
	if (fStrainDispOpt != kStandardB)
	{
		cout << "\n FiniteStrainT::FiniteStrainT: no strain-displacement options\n" << endl;
		fStrainDispOpt = kStandardB;
	}
}

/* destructor */
FiniteStrainT::~FiniteStrainT(void)
{
	delete fFSMatSupport;
}

/* called immediately after constructor */
void FiniteStrainT::Initialize(void)
{
	/* inherited */
	SolidElementT::Initialize();

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
		fF_all.Dimension(nip*nsd*nsd);
		fF_List.Dimension(nip);
		for (int i = 0; i < nip; i++)
			fF_List[i].Set(nsd, nsd, fF_all.Pointer(i*nsd*nsd));
	}
	
	/* allocate "last" deformation gradient list */
	if (need_F_last)
	{
		int nip = NumIP();
		int nsd = NumSD();
		fF_last_all.Dimension(nip*nsd*nsd);
		fF_last_List.Dimension(nip);
		for (int i = 0; i < nip; i++)
			fF_last_List[i].Set(nsd, nsd, fF_last_all.Pointer(i*nsd*nsd));
	}
}

/* TEMPORARY */
void FiniteStrainT::InitialCondition(void)
{
	/* inherited */
	SolidElementT::InitialCondition();
	
	/* set the source for the iteration number */
	fFSMatSupport->SetIterationNumber(ElementSupport().IterationNumber(Group()));
}

/* compute field gradients with respect to current coordinates */
void FiniteStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const
{
	if (fCurrShapes)
	{
		/* field gradient */
		fCurrShapes->GradU(u, grad_u);
	}
	else
	{
		cout << "\n FiniteStrainT::ComputeGradient: shape functions wrt current coords not defined" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/* compute field gradients with respect to current coordinates */
void FiniteStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, 
	int ip) const
{
	if (fCurrShapes)
	{
		/* field gradient */
		fCurrShapes->GradU(u, grad_u, ip);
	}
	else
	{
		cout << "\n FiniteStrainT::ComputeGradient: shape functions wrt current coords not defined" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/***********************************************************************
* Protected
***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* FiniteStrainT::NewMaterialSupport(MaterialSupportT* p) const
{
	/* allocate */
	if (!p) p = new FSMatSupportT(NumSD(), NumDOF(), NumIP());

	/* inherited initializations */
	SolidElementT::NewMaterialSupport(p);
	
	/* set FiniteStrainT fields */
	FSMatSupportT* ps = TB_DYNAMIC_CAST(FSMatSupportT*, p);
	if (ps) {
		ps->SetDeformationGradient(&fF_List);
		ps->SetDeformationGradient_last(&fF_last_List);
	}

	return p;
}

/* construct materials manager and read data */
MaterialListT* FiniteStrainT::NewMaterialList(int nsd, int size)
{
	if (size > 0)
	{
		/* material support */
		if (!fFSMatSupport) {
			fFSMatSupport = TB_DYNAMIC_CAST(FSMatSupportT*, NewMaterialSupport());
			if (!fFSMatSupport) ExceptionT::GeneralFail("FiniteStrainT::NewMaterialList");
		}

		if (nsd == 1)
			return new SolidMatList1DT(size, *fFSMatSupport);
		else if (nsd == 2)
			return new SolidMatList2DT(size, *fFSMatSupport);
		else if (nsd == 3)
			return new SolidMatList3DT(size, *fFSMatSupport);
		else
			return NULL;
	}
	else
	{
		if (nsd == 1)
			return new SolidMatList1DT;
		else if (nsd == 2)
			return new SolidMatList2DT;
		else if (nsd == 3)
			return new SolidMatList3DT;
		else
			return NULL;
	}
}

/* construct list of materials from the input stream */
void FiniteStrainT::ReadMaterialData(ifstreamT& in)
{
	/* inherited */
	SolidElementT::ReadMaterialData(in);

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
		FSSolidMatT* mat = (FSSolidMatT*) pcont_mat;

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
	SolidElementT::SetGlobalShape();

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

/* write all current element information to the stream */
void FiniteStrainT::CurrElementInfo(ostream& out) const
{
	/* inherited */
	SolidElementT::CurrElementInfo(out);
	
	/* write deformation gradients */
	out << "\n i.p. deformation gradients:\n";
	for (int i = 0; i < fF_List.Length(); i++)
		out << " ip: " << i+1 << '\n'
		    << fF_List[i] << '\n';
	out << '\n';
}
