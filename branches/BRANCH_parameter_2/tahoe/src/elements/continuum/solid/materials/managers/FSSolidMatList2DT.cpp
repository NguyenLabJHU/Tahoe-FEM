/* $Id: FSSolidMatList2DT.cpp,v 1.1.2.5 2004-03-04 20:11:54 paklein Exp $ */
#include "FSSolidMatList2DT.h"
#include "FSMatSupportT.h"

#include "fstreamT.h"
#include "SolidMaterialsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#endif

#include "FDKStV2D.h"
#include "FDCubic2DT.h"
#include "SimoIso2D.h"
#include "QuadLog2D.h"
#include "QuadLogOgden2DT.h"

#ifdef CAUCHY_BORN_MATERIAL
#include "EAMFCC2D.h"
#include "LJTr2D.h"
#include "Hex2D.h"
#endif

#ifdef MODCBSW_MATERIAL
#include "ModCB2DT.h"
#include "SWDiamond110.h"
#include "SWDiamond100.h"
#endif

#ifdef VIB_MATERIAL
#include "IsoVIB2D.h"
#include "J2IsoVIB2DLinHardT.h"
#include "VIB2D.h"
#include "D2VIB2D_a.h"
#include "OgdenIsoVIB2D.h"
#endif

#ifdef PLASTICITY_MACRO_MATERIAL
#include "HyperEVP2D.h"
#include "BCJHypo2D.h"
#include "BCJHypoIsoDamageKE2D.h"
#include "BCJHypoIsoDamageYC2D.h"
#endif

#ifdef PLASTICITY_CRYSTAL_MATERIAL
#include "LocalCrystalPlast2D.h"
#include "GradCrystalPlast2D.h"
#include "LocalCrystalPlastFp2D.h"
#include "GradCrystalPlastFp2D.h"
#endif

#ifdef PLASTICITY_J2_MATERIAL
#include "J2Simo2D.h"
#include "J2QL2DLinHardT.h"
#endif

#ifdef VISCOELASTIC_MATERIALS_DEV
#include "RGVIB2D.h"
#include "RGSplit3D.h"
#include "FDSV_KStV2D.h"
#endif

#ifdef VISCOELASTICITY
#include "RGSplitT.h"
#endif

#ifdef ELASTIC_OGDEN_MATERIAL_DEV
#include "OgdenMaterialT.h"
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
#include "ABAQUS_BCJ.h"
#include "ABAQUS_BCJ_ISO.h"
#include "ABAQUS_VUMAT_BCJ.h"
#endif
#endif

#ifdef THERMO_VISCO_PLASTIC_MATERIAL
#include "tevp2D.h"
#include "povirk2D.h"
#endif

using namespace Tahoe;

/* constructor */
FSSolidMatList2DT::FSSolidMatList2DT(int length, const FSMatSupportT& support):
	SolidMatListT(length, support),
	fFSMatSupport(&support)
{
#pragma message("check spatial dimension of material support")
	SetName("large_strain_material_2D");
}

FSSolidMatList2DT::FSSolidMatList2DT(void):
	fFSMatSupport(NULL)
{
	SetName("large_strain_material_2D");
}

/* read material data from the input stream */
void FSSolidMatList2DT::ReadMaterialData(ifstreamT& in)
{
	const char caller[] = "FSSolidMatList2DT::ReadMaterialData";

	int i, matnum;
	SolidT::TypeT matcode;
	try {

	/* read material data */
	for (i = 0; i < fLength; i++)
 	{
		in >> matnum; matnum--;
		in >> matcode;
		/* checks */
		if (matnum < 0  || matnum >= fLength) ExceptionT::BadInputValue(caller);
		
		/* repeated material number */
		if (fArray[matnum] != NULL)
			ExceptionT::BadInputValue(caller, "repeated material number: %d", matnum + 1);
		
		/* add to the list of materials */
		switch (matcode)
		{
			case kFDKStV:
			{
				fArray[matnum] = new FDKStV2D(in, *fFSMatSupport);
				break;
			}
			case kFDCubic:
			{
				fArray[matnum] = new FDCubic2DT(in, *fFSMatSupport);
				break;
			}
			case kSimoIso:
			{
				fArray[matnum] = new SimoIso2D(in, *fFSMatSupport);
				break;
			}
			case kQuadLog:
			{
				fArray[matnum] = new QuadLog2D(in, *fFSMatSupport);
				break;
			}
			case kQuadLogOgden:
			{
				fArray[matnum] = new QuadLogOgden2DT(in, *fFSMatSupport);												
				break;
			}
			case kJ2Simo:
			{
#ifdef PLASTICITY_J2_MATERIAL
				fArray[matnum] = new J2Simo2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kJ2QL:
			{
#ifdef PLASTICITY_J2_MATERIAL
				fArray[matnum] = new J2QL2DLinHardT(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kMRSSKStV:
			{
#ifdef PLASTICITY_MR_MATERIAL_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new MRSSKStV2D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_MR_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLJTr2D:
			{
#ifdef CAUCHY_BORN_MATERIAL
				fArray[matnum] = new LJTr2D(in, *fFSMatSupport);
				break;
#else
				ExceptionT::BadInputValue(caller, "CAUCHY_BORN_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kHex2D:
			{
#ifdef CAUCHY_BORN_MATERIAL
				fArray[matnum] = new Hex2D(in, *fFSMatSupport);
				break;
#else
				ExceptionT::BadInputValue(caller, "CAUCHY_BORN_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kFCCEAM:
			{
#ifdef CAUCHY_BORN_MATERIAL
				int i_plane_code;
				in >> i_plane_code;
				EAMFCC2D::PlaneCodeT plane_code = (EAMFCC2D::PlaneCodeT) i_plane_code;
				if (plane_code != EAMFCC2D::kFCC001 &&
				    plane_code != EAMFCC2D::kFCC101 &&
				    plane_code != EAMFCC2D::kFCC111)
				{
					ExceptionT::BadInputValue(caller, "unrecognized plane code: %d", i_plane_code);
				}

				fArray[matnum] = new EAMFCC2D(in, *fFSMatSupport, plane_code);			
				break;
#else
				ExceptionT::BadInputValue(caller, "CAUCHY_BORN_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kmodCauchyBornDC:
			{
#ifdef MODCBSW_MATERIAL
				int i_plane_code;
				in >> i_plane_code;
				ModCB2DT::PlaneCodeT plane_code = (ModCB2DT::PlaneCodeT) i_plane_code;
				if (plane_code != ModCB2DT::kDC001 &&
				    plane_code != ModCB2DT::kDC101 &&
				    plane_code != ModCB2DT::kDC111)
				{
					ExceptionT::BadInputValue(caller, "unrecognized plane code: %d", i_plane_code);
				}

				fArray[matnum] = new ModCB2DT(in, *fFSMatSupport, true, plane_code);
				break;
#else
				ExceptionT::BadInputValue(caller, "MODCBSW_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kVIB:
			{
#ifdef VIB_MATERIAL
				fArray[matnum] = new VIB2D(in, *fFSMatSupport);
				fHasLocalizers = true;  				
				break;
#else
				ExceptionT::BadInputValue(caller, "VIB_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kIsoVIBSimo:
			{
#ifdef VIB_MATERIAL
				fArray[matnum] = new IsoVIB2D(in, *fFSMatSupport);
				fHasLocalizers = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VIB_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kIsoVIBOgden:
			{
#ifdef VIB_MATERIAL
				fArray[matnum] = new OgdenIsoVIB2D(in, *fFSMatSupport);
				fHasLocalizers = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VIB_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kIsoVIBSimoJ2:
			{
#ifdef VIB_MATERIAL
				fArray[matnum] = new J2IsoVIB2DLinHardT(in, *fFSMatSupport);
				fHasLocalizers = true;
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VIB_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kFCC:
			{
				ExceptionT::BadInputValue(caller, "material %d is 3D only", matcode);
			}
			case kFossumSSIso:
			{
#ifdef FOSSUM_MATERIAL_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new FossumSSIso2DT(in, *fSSMatSupport);
				fHasHistory = true;		       		      				
				break;
#else
				ExceptionT::BadInputValue(caller, "FOSSUM_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kRGSplitVE:
			{
#ifdef VISCOELASTICITY
				fArray[matnum] = new RGSplitT(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VISCOELASTICITY not enabled: %d", matcode);
#endif
			}
			case kThermoViscoPlastic:
			{
#ifdef THERMO_VISCO_PLASTIC_MATERIAL
				fArray[matnum] = new tevp2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "THERMO_VISCO_PLASTIC_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kPovirk2D:
			{
#ifdef THERMO_VISCO_PLASTIC_MATERIAL
				fArray[matnum] = new povirk2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "THERMO_VISCO_PLASTIC_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kHyperEVP:
			{
#ifdef PLASTICITY_MACRO_MATERIAL
				fArray[matnum] = new HyperEVP2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_MACRO_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kBCJHypo:
			{
#ifdef PLASTICITY_MACRO_MATERIAL
				fArray[matnum] = new BCJHypo2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_MACRO_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kBCJHypoIsoDmgKE:
			{
#ifdef PLASTICITY_MACRO_MATERIAL
				fArray[matnum] = new BCJHypoIsoDamageKE2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_MACRO_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kBCJHypoIsoDmgYC:
			{
#ifdef PLASTICITY_MACRO_MATERIAL
				fArray[matnum] = new BCJHypoIsoDamageYC2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_MACRO_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLocXtalPlast:
			{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
				fArray[matnum] = new LocalCrystalPlast2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLocXtalPlast_C:
			{
				ExceptionT::BadInputValue(caller, "model %d is not implemented in 2D", matcode);
			}
			case kGrdXtalPlast:
			{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
				fArray[matnum] = new GradCrystalPlast2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLocXtalPlastFp:
			{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
				fArray[matnum] = new LocalCrystalPlastFp2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kGrdXtalPlastFp:
			{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
				fArray[matnum] = new GradCrystalPlastFp2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kABAQUS_BCJ_ISO:
			{
#ifdef __F2C__
#if defined(ABAQUS_MATERIAL) && defined(ABAQUS_BCJ_MATERIAL_DEV)
	
				/* large strain interface */
				fArray[matnum] = new ABAQUS_BCJ_ISO(in, *fFSMatSupport);
					
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "ABAQUS_MATERIAL or ABAQUS_BCJ_MATERIAL_DEV not enabled: %d", matcode);
#endif
#else
				ExceptionT::BadInputValue(caller, "model requires f2c support: %d", kABAQUS_BCJ_ISO);
#endif /* __F2C__ */	
			}			
			case kABAQUS_VUMAT_BCJ:
			{
#ifdef __F2C__
#if defined(ABAQUS_MATERIAL) && defined(ABAQUS_BCJ_MATERIAL_DEV)
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new ABAQUS_VUMAT_BCJ(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "ABAQUS_MATERIAL or ABAQUS_BCJ_MATERIAL_DEV not enabled: %d", matcode);
#endif
#else
				ExceptionT::BadInputValue(caller, "model requires f2c support: %d", kABAQUS_BCJ);
#endif /* __F2C__ */
			}
			case kRGVIB:
			{
#ifdef VISCOELASTIC_MATERIALS_DEV
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new RGVIB2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VISCOELASTIC_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			case kRGSplit:
			{
#ifdef VISCOELASTIC_MATERIALS_DEV
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new RGSplit3D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VISCOELASTIC_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			case kFDSVKStV:
			{
#ifdef VISCOELASTIC_MATERIALS_DEV
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDSV_KStV2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VISCOELASTIC_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			case kSSSVKStV:
			{
#ifdef VISCOELASTIC_MATERIALS_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new SSSV_KStV2D(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VISCOELASTIC_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			case kOgdenMat:
			{
#ifdef ELASTIC_OGDEN_MATERIAL_DEV
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new OgdenMaterialT(in, *fFSMatSupport);
				break;
#else
				ExceptionT::BadInputValue(caller, "ELASTIC_OGDEN_MATERIAL_DEV not enabled: %d", matcode);
#endif
			}
			case kSSJ2LinHard:
			{
#ifdef J2PLASTICITY_MATERIALS_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new SSJ2LinHard2D(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "J2PLASITICITY_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			case kSSJ2LinHardplane:
			{
#ifdef J2PLASTICITY_MATERIALS_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new SSJ2LinHard3Dplane(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "J2PLASITICITY_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
//TEMP
#if 0
			case kSWDC100:
				fArray[matnum] = new SWDiamond100(in, *fFSMatSupport);
				break;

			case kSWDC110:
				fArray[matnum] = new SWDiamond110(in, *fFSMatSupport);
				break;


			case kD2VIB:
			{
#ifdef __NO_RTTI__
				cout << "\n FSSolidMatList2DT::ReadMaterialData: need RTTI" << endl;
				throw ExceptionT::kBadInputValue;
#endif
				const D2MeshFreeFSSolidT* D2ElementGroup;
				D2ElementGroup = TB_DYNAMIC_CAST(const D2MeshFreeFSSolidT*, &fElementGroup);
				if (!D2ElementGroup)
				{
					cout << "\n FSSolidMatList2DT::ReadMaterialData: material " << kD2VIB << " requires\n"
					     <<   "     a higher order gradient element" << endl;
					throw ExceptionT::kBadInputValue;
				}
			
				fArray[matnum] = new D2VIB2D_a(in, *D2ElementGroup);
				fHasLocalizers = true;  				
				break;
			}

#endif //TEMP

			default:
				ExceptionT::BadInputValue(caller, "unknown material code: %d", matcode);
		}

		/* safe cast since all structural */
		SolidMaterialT* pmat = (SolidMaterialT*) fArray[matnum];

		/* verify construction */
		if (!pmat) throw ExceptionT::kOutOfMemory;
		
		/* set thermal LTf pointer */
		int LTfnum = pmat->ThermalStrainSchedule();
		if (LTfnum > -1)
		{
			pmat->SetThermalSchedule(fSolidMatSupport->Schedule(LTfnum));
			
			/* set flag */
			fHasThermal = true;
		}
		
		/* perform initialization */
		pmat->Initialize();
	}  } /* end try */
	
	catch (ExceptionT::CodeT error)
	{
		ExceptionT::Throw(error, caller, "exception constructing material %d, index %d, code %d",
			i+1, matnum+1, matcode);
	}
}

/* return true if the list contains plane stress models */
bool FSSolidMatList2DT::HasPlaneStress(void) const
{
	/* check materials */
	for (int i = 0; i < Length(); i++)
	{
		/* get pointer to Material2DT */
		const ContinuumMaterialT* cont_mat = fArray[i];
		const SolidMaterialT* sol_mat = TB_DYNAMIC_CAST(const SolidMaterialT*, cont_mat);
		
		/* assume materials that don't have Material2DT are plane strain */
		if (sol_mat && sol_mat->Constraint() == SolidMaterialT::kPlaneStress) 
			return true;
	}
	return false;
}

/* information about subordinate parameter lists */
void FSSolidMatList2DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("fs_material_list_2D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FSSolidMatList2DT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
		SubListT& sub_sub_list) const
{
	if (sub == "fs_material_list_2D")
	{
		order = ParameterListT::Choice;
	
		sub_sub_list.AddSub("large_strain_cubic_2D");
		sub_sub_list.AddSub("large_strain_StVenant_2D");		
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSSolidMatList2DT::NewSub(const StringT& list_name) const
{
	/* try to construct material */
	FSSolidMatT* fs_solid_mat = NewFSSolidMat(list_name);
	if (fs_solid_mat)
		return fs_solid_mat;
	else /* inherited */
		return SolidMatListT::NewSub(list_name);
}

/* accept parameter list */
void FSSolidMatList2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SolidMatListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	AutoArrayT<FSSolidMatT*> materials;
	const ArrayT<ParameterListT>& subs = list.Lists();
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		FSSolidMatT* mat = NewFSSolidMat(sub.Name());
		if (mat) {
			materials.Append(mat);
			mat->TakeParameterList(sub);
		}
	}

	/* transfer */
	Dimension(materials.Length());
	for (int i = 0; i < materials.Length(); i++)
		fArray[i] = materials[i];
}

/* construct the specified material or NULL if the request cannot be completed */
FSSolidMatT* FSSolidMatList2DT::NewFSSolidMat(const StringT& name) const
{
	FSSolidMatT* mat = NULL;

	if (name == "large_strain_cubic_2D")
		mat = new FDCubic2DT;
	else if (name == "large_strain_StVenant_2D")
		mat = new FDKStV2D;

	/* set support */
	if (mat) mat->SetFSMatSupport(fFSMatSupport);

	return mat;

}
