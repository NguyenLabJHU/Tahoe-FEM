/* $Id: FiniteStrainAxiT.cpp,v 1.4 2004-07-15 08:26:27 paklein Exp $ */
#include "FiniteStrainAxiT.h"

#include "ShapeFunctionT.h"
#include "FSSolidMatT.h"
#include "FSMatSupportT.h"
#include "ParameterContainerT.h"

using namespace Tahoe;

const double Pi = acos(-1.0);
const int kRadialDirection = 0; /* the x direction is radial */
const int kNSD = 2;

/* constructor */
FiniteStrainAxiT::FiniteStrainAxiT(const ElementSupportT& support):
	FiniteStrainT(support),
	fMat2D(kNSD),
	fLocCurrCoords(LocalArrayT::kCurrCoords)	
{
	SetName("large_strain_axi");
}

/* information about subordinate parameter lists */
void FiniteStrainAxiT::DefineSubs(SubListT& sub_list) const
{
	const char caller[] = "FiniteStrainAxiT::DefineSubs";

	/* inherited */
	FiniteStrainT::DefineSubs(sub_list);

	/* remove previous block definition */
	const char old_block_def[] = "large_strain_element_block";
	if (!sub_list.RemoveSub(old_block_def))
		ExceptionT::GeneralFail(caller, "did not find \"%s\"", old_block_def);

	/* element block/material specification */
	sub_list.AddSub("large_strain_axi_element_block", ParameterListT::OnePlus);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FiniteStrainAxiT::NewSub(const StringT& name) const
{
	if (name == "large_strain_axi_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("large_strain_material_3D");
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	}
	else /* inherited */
		return FiniteStrainT::NewSub(name);
}

/* accept parameter list */
void FiniteStrainAxiT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "FiniteStrainAxiT::TakeParameterList";

	/* inherited */
	FiniteStrainT::TakeParameterList(list);

	/* dimensions */
	int nip = NumIP();

	/* integration point radii over the current element */
	fRadius_X.Dimension(nip);
	fRadius_x.Dimension(nip);

	/* redimension with out-of-plane component */
	int nstrs = dSymMatrixT::NumValues(kNSD) + 1;
	fD.Dimension(nstrs);
	fB.Dimension(nstrs, NumSD()*NumElementNodes());

	/* allocate 3D deformation gradient list */
	if (fF_List.Length() > 0) {
		fF_all.Dimension(nip*3*3);
		for (int i = 0; i < nip; i++)
			fF_List[i].Set(3, 3, fF_all.Pointer(i*3*3));
	}
	
	/* allocate 3D "last" deformation gradient list */
	if (fF_last_List.Length() > 0) {
		fF_last_all.Dimension(nip*3*3);
		for (int i = 0; i < nip; i++)
			fF_last_List[i].Set(3, 3, fF_last_all.Pointer(i*3*3));
	}
}

/* extract the list of material parameters */
void FiniteStrainAxiT::CollectMaterialInfo(const ParameterListT& all_params, ParameterListT& mat_params) const
{
	const char caller[] = "FiniteStrainAxiT::CollectMaterialInfo";
	
	/* initialize */
	mat_params.Clear();
	mat_params.SetName("large_strain_material_3D");
	
	/* collected material parameters */
	int num_blocks = all_params.NumLists("large_strain_axi_element_block");
	for (int i = 0; i < num_blocks; i++) {

		/* block information */	
		const ParameterListT& block = all_params.GetList("large_strain_axi_element_block", i);
		
		/* collect material parameters */
		const ParameterListT& mat_list = block.GetList(mat_params.Name());
		const ArrayT<ParameterListT>& mat = mat_list.Lists();
		mat_params.AddList(mat[0]);
	}
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* initial local arrays */
void FiniteStrainAxiT::SetLocalArrays(void)
{
	/* inherited */
	FiniteStrainT::SetLocalArrays();

	/* allocate and set source */
	fLocCurrCoords.Dimension(NumElementNodes(), NumSD());
	ElementSupport().RegisterCoordinates(fLocCurrCoords);
}

/* construct a new material support and return a pointer */
MaterialSupportT* FiniteStrainAxiT::NewMaterialSupport(MaterialSupportT* p) const
{
	/* construct 3D support */
	if (!p) {
		p = new FSMatSupportT(NumDOF(), NumIP());
		p->SetNumSD(3);
	}

	/* inherited initializations */
	FiniteStrainT::NewMaterialSupport(p);

	return p;
}

/* form shape functions and derivatives */
void FiniteStrainAxiT::SetGlobalShape(void)
{
	/* skip call to FiniteStrainT::SetGlobalShape */
	SolidElementT::SetGlobalShape();

	/* what needs to get computed */
	int material_number = CurrentElement().MaterialNumber();
	bool needs_F = Needs_F(material_number);
	bool needs_F_last = Needs_F_last(material_number);
	
	/* get current element coordinates */
	SetLocalX(fLocCurrCoords);
	int nen = fLocCurrCoords.NumberOfNodes();
	int nun = fLocDisp.NumberOfNodes();
	
	/* loop over integration points */
	for (int i = 0; i < NumIP(); i++)
	{
		/* compute radii */
		const double* NaX = fShapes->IPShapeX(i);
		const double* NaU = fShapes->IPShapeU(i);
		const double* X_r = fLocInitCoords(kRadialDirection);
		const double* u_r = fLocDisp(kRadialDirection);
		const double* u_r_last = (needs_F_last) ? fLocLastDisp(kRadialDirection) : u_r; /* fLocLastDisp not used */
		double R = 0.0;
		double u = 0.0;
		double u_last = 0.0;
		if (nen == nun) 
		{
			for (int a = 0; a < nen; a++) {
				R += (*NaX)*(*X_r++);
				u += (*NaU)*(*u_r++);
				u_last += (*NaU)*(*u_r_last++);
				NaX++;
				NaU++;
			}
		}
		else /* separate loops for field and geometry */
		{
			for (int a = 0; a < nen; a++) {
				R += (*NaX)*(*X_r++);
				NaX++;
			}
			for (int a = 0; a < nun; a++) {
				u += (*NaU)*(*u_r++);
				u_last += (*NaU)*(*u_r_last++);
				NaU++;
			}		
		}
		double r = R + u;
		fRadius_X[i] = R;
		fRadius_x[i] = r;

		/* deformation gradient */
		if (needs_F)
		{
			/* 2D deformation gradient */
			fShapes->GradU(fLocDisp, fMat2D, i);
			fMat2D.PlusIdentity();

			/* make axisymmetric */
			dMatrixT& F3D = fF_List[i];
			F3D.Rank2ExpandFrom2D(fMat2D);
			F3D(2,2) = r/R;
		}

		/* "last" deformation gradient */
		if (needs_F_last)
		{
			/* 2D deformation gradient */
			fShapes->GradU(fLocLastDisp, fMat2D, i);
			fMat2D.PlusIdentity();

			/* make axisymmetric */
			dMatrixT& F3D = fF_last_List[i];
			F3D.Rank2ExpandFrom2D(fMat2D);
			F3D(2,2) = (R + u_last)/R;
		}
	}
}
