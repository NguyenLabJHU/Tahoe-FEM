/* $Id: CSEAnisoT.cpp,v 1.3 2001-02-27 00:06:24 paklein Exp $ */
/* created: paklein (11/19/1997)                                          */
/* Cohesive surface elements with scalar traction potentials,             */
/* i.e., the traction potential is a function of the gap magnitude,       */
/* or effective gap magnitude only.                                       */

#include "CSEAnisoT.h"

#include <math.h>
#include <iostream.h>
#include <iomanip.h>

#include "fstreamT.h"
#include "Constants.h"
#include "SurfaceShapeT.h"
#include "FEManagerT.h"
#include "NodeManagerT.h"
#include "SurfacePotentialT.h"
#include "eControllerT.h"

/* potential functions */
#include "XuNeedleman2DT.h"
#include "XuNeedleman3DT.h"
#include "TvergHutch2DT.h"

//TEMP - catch NaN's
#include "math_utils.h"

/* constructor */
CSEAnisoT::CSEAnisoT(FEManagerT& fe_manager, bool rotate):
	CSEBaseT(fe_manager),
	fRotate(rotate),
	fCurrShapes(NULL),
	fQ(fNumSD),
	fdelta(fNumSD),
	fT(fNumSD),
	fddU(fNumSD)
{
	/* reset format for the element stiffness matrix */
	if (fRotate) fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
}

/* destructor */
CSEAnisoT::~CSEAnisoT(void)
{
	if (fRotate)
	{
		delete fCurrShapes;
		fCurrShapes = NULL;
	}
}

/* form of tangent matrix */
GlobalT::SystemTypeT CSEAnisoT::TangentType(void) const
{
	if (fRotate)
		/* tangent matrix is not symmetric */
		return GlobalT::kNonSymmetric;
	else
		return GlobalT::kSymmetric;
}

void CSEAnisoT::Initialize(void)
{
	/* inherited */
	CSEBaseT::Initialize();
	
	/* rotating local frame */
	if (fRotate)
	{
		/* shape functions wrt. current coordinates (linked parent domains) */
		fCurrShapes = new SurfaceShapeT(*fShapes, fLocCurrCoords);
		if (!fCurrShapes) throw eOutOfMemory;
		fCurrShapes->Initialize();
		
		/* allocate work space */
		fnsd_nee_1.Allocate(fNumSD, fNumElemEqnos);
		fnsd_nee_2.Allocate(fNumSD, fNumElemEqnos);
		fdQ.Allocate(fNumSD);
		for (int k = 0; k < fNumSD; k++)
			fdQ[k].Allocate(fNumSD, fNumElemEqnos);
	}
	else
		fCurrShapes = fShapes;

	/* streams */
	ifstreamT& in = fFEManager.Input();
	ostream&   out = fFEManager.Output();
		
	/* construct props */
	int numprops;
	in >> numprops;
	fSurfPots.Allocate(numprops);
	for (int i = 0; i < fSurfPots.Length(); i++)
	{
		int num, code;
		in >> num >> code;
		num--;

		/* check for repeated number */
		if (fSurfPots[num] != NULL) throw eBadInputValue;

		switch (code)
		{
			case SurfacePotentialT::kXuNeedleman:
			{			
				if (fNumDOF == 2)
					fSurfPots[num] = new XuNeedleman2DT(in);
				else
					fSurfPots[num] = new XuNeedleman3DT(in);
				break;
			}
			case SurfacePotentialT::kTvergaardHutchinson:
			{
				if (fNumDOF == 2)
					fSurfPots[num] = new TvergHutch2DT(in);
				else
				{
					cout << "\n CSEAnisoT::Initialize: potential not implemented for 3D: "
					     << code << endl; 				
					throw eBadInputValue;
				}
				break;
			}
			default:
				cout << "\n CSEAnisoT::Initialize: unknown potential code: " << code << endl;
				throw eBadInputValue;
		}
		if (!fSurfPots[num]) throw eOutOfMemory;
	}

	/* check compatibility of constitutive outputs */
	if (fSurfPots.Length() > 1 && fNodalOutputCodes[MaterialData])
		for (int k = 0; k < fSurfPots.Length(); k++)
		{
			const SurfacePotentialT* pot_k = fSurfPots[k];
			for (int i = k+1; i < fSurfPots.Length(); i++)
			{
				const SurfacePotentialT* pot_i = fSurfPots[i];
				if (!SurfacePotentialT::CompatibleOutput(*pot_k, *pot_i))
				{
					cout << "\n CSEAnisoT::Initialize: incompatible output between potentials\n"
					     <<   "     " << k+1 << " and " << i+1 << endl;
					throw eBadInputValue;
				}
			}
		}
			
	/* write */
	out << " Rotating local coordinate frame . . . . . . . . = " <<
	    ((fRotate) ? "ACTIVE" : "INACTIVE") << '\n';
	out << "\n Cohesive surface potentials:\n";
	out << " Number of potentials. . . . . . . . . . . . . . = ";
	out << fSurfPots.Length() << '\n';
	for (int j = 0; j < fSurfPots.Length(); j++)
	{
		out << "\n Potential number. . . . . . . . . . . . . . . . = " << j + 1 << '\n';
		out << " Potential name:\n";
		fSurfPots[j]->PrintName(out);
		fSurfPots[j]->Print(out);
	}
}

/***********************************************************************
* Protected
***********************************************************************/

void CSEAnisoT::LHSDriver(void)
{
	/* matrix format */
	dMatrixT::SymmetryFlagT format = (fRotate) ?
		dMatrixT::kWhole : dMatrixT::kUpperOnly;

	/* time-integration parameters */
	double constK = 0.0;
	int formK = fController->FormK(constK);
	if (!formK) return;

	/* node map of facet 1 */
	iArrayT facet1;
	(fShapes->NodesOnFacets()).RowAlias(0, facet1);

	Top();
	while ( NextElement() )
	{
		/* current element */
		const ElementCardT& element = CurrentElement();
	
		/* surface potential */
		SurfacePotentialT* surfpot = fSurfPots[element.MaterialNumber()];
	
		/* get ref geometry (1st facet only) */
		fNodes1.Collect(facet1, element.NodesX());
		fLocInitCoords1.SetLocal(fNodes1);

		/* get current geometry */
		SetLocalX(fLocCurrCoords); //EFFECTIVE_DVA

		/* initialize */
		fLHS = 0.0;

		/* loop over integration points */
		fShapes->TopIP();
		while (fShapes->NextIP())
		{
			/* integration weights */
			double w = fShapes->IPWeight();		

			/* coordinate transformations */
			double j0, j;
			if (fRotate)
			{
				j0 = fShapes->Jacobian();
				j  = fCurrShapes->Jacobian(fQ, fdQ);
			}
			else
				j0 = j = fShapes->Jacobian(fQ);

			/* check */
			if (j0 <= 0.0 || j <= 0.0) throw eBadJacobianDet;
		
			/* gap vector and gradient (facet1 to facet2) */
			const dArrayT&    delta = fShapes->InterpolateJumpU(fLocCurrCoords);
			const dMatrixT& d_delta = fShapes->Grad_d();

			/* stiffness in local frame */
			fQ.MultTx(delta, fdelta);
			const dMatrixT& K = surfpot->Stiffness(fdelta);
			
			/* rotation */
			if (fRotate)
			{
				/* traction in local frame */
				const dArrayT& T = surfpot->Traction(fdelta);

				/* 1st term */
				fT.SetToScaled(j0*w*constK, T);
				Q_ijk__u_j(fdQ, fT, fnsd_nee_1);
				fNEEmat.MultATB(d_delta, fnsd_nee_1);
				fLHS += fNEEmat;

				/* 2st term */
				fddU.SetToScaled(j0*w*constK, K);
				fnsd_nee_1.MultATB(fQ, d_delta);
				fnsd_nee_2.MultATB(fddU, fnsd_nee_1);
				u_i__Q_ijk(delta, fdQ, fnsd_nee_1);
				fNEEmat.MultATB(fnsd_nee_2, fnsd_nee_1);
				fLHS += fNEEmat;
			}

			/* 3rd term */
			fddU.MultQBQT(fQ, K);
			fddU *= j0*w*constK;
			fLHS.MultQTBQ(d_delta, fddU, format, dMatrixT::kAccumulate);	
		}
		
		/* assemble */
		AssembleLHS();
	}
}

void CSEAnisoT::RHSDriver(void)
{
	/* time-integration parameters */
	double constKd = 0.0;
	int formKd = fController->FormKd(constKd);
	if (!formKd) return;

	/* node map of facet 1 */
	iArrayT facet1;
	(fShapes->NodesOnFacets()).RowAlias(0, facet1);

	/* fracture surface area */
	fFractureArea = 0.0;

	Top();
	while (NextElement())
	{
		/* current element */
		ElementCardT& element = CurrentElement();
	
		/* get ref geometry (1st facet only) */
		fNodes1.Collect(facet1, element.NodesX());
		fLocInitCoords1.SetLocal(fNodes1);
	
		if (element.Flag() != kOFF)
		{
			/* surface potential */
			SurfacePotentialT* surfpot = fSurfPots[element.MaterialNumber()];
			
			/* get current geometry */
			SetLocalX(fLocCurrCoords); //EFFECTIVE_DVA
	
	  		/* initialize */
	  		fRHS = 0.0;
			
			/* loop over integration points */
			int all_failed = 1;
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* integration weights */
				double w = fShapes->IPWeight();		

				/* coordinate transformations */
				double j0, j;
				if (fRotate)
				{
					j0 = fShapes->Jacobian();
					j  = fCurrShapes->Jacobian(fQ);
				}
				else
					j0 = j = fShapes->Jacobian(fQ);
				
				/* check */
				if (j0 <= 0.0 || j <= 0.0)
				{
					cout << "\n CSEAnisoT::RHSDriver: jacobian error" << endl;
					throw eBadJacobianDet;
				}
	
				/* gap vector from facet1 to facet2 */
				const dArrayT& delta = fShapes->InterpolateJumpU(fLocCurrCoords);
	
				/* gap -> traction, in/out of local frame */
				fQ.MultTx(delta, fdelta);
				fQ.Multx(surfpot->Traction(fdelta), fT);
	
				/* expand */
				fShapes->Grad_d().MultTx(fT, fNEEvec);
	
				/* accumulate */
				fRHS.AddScaled(-j0*w*constKd, fNEEvec);
				
				/* check status */
				SurfacePotentialT::StatusT status = surfpot->Status(fdelta);
				if (status != SurfacePotentialT::Failed) all_failed = 0;
				
				/* fracture area */
				if (fOutputArea && status != SurfacePotentialT::Precritical)
					fFractureArea += j0*w;
			}

			//TEMP - catch NaN's
			for (int i = 0; i < fRHS.Length(); i++)
			{
				double x = fRHS[i];			
				if (is_NaN(x))
				{
					cout << "\n CSEAnisoT::RHSDriver: NaN detected" << endl;
					throw eBadJacobianDet;
				}
				else if (is_Inf(x))			
				{
					cout << "\n CSEAnisoT::RHSDriver: Inf detected" << endl;
					throw eBadJacobianDet;
				}
			}

			/* assemble */
			AssembleRHS();
			
			/* mark elements */
			if (all_failed)
			{
				int& flag = element.Flag();
				if (flag == kON) flag = kMarked;
			}
		}
		else if (fOutputArea)
		{
			/* integrate fracture area */
			fShapes->TopIP();
			while (fShapes->NextIP())
				fFractureArea += (fShapes->Jacobian())*(fShapes->IPWeight());
		}
	}	
}

/* nodal value calculations */
void CSEAnisoT::SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
	iArrayT& counts) const
{
	/* inherited */
	CSEBaseT::SetNodalOutputCodes(mode, flags, counts);

	/* resize for vectors not magnitudes */
	if (flags[NodalDispJump] == mode)
		counts[NodalDispJump] = fNumDOF;	
	if (flags[NodalTraction] == mode)
		counts[NodalTraction] = fNumDOF;
	if (flags[MaterialData] == mode)
		counts[MaterialData] = fSurfPots[0]->NumOutputVariables();
}

/* extrapolate the integration point stresses and strains and extrapolate */
void CSEAnisoT::ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	const iArrayT& e_codes, dArray2DT& e_values)
{
	/* number of output values */
	int n_out = n_codes.Sum();
	int e_out = e_codes.Sum();

	/* nothing to output */
	if (n_out == 0 && e_out == 0) return;

	/* reset averaging workspace */
	fNodes->ResetAverage(n_out);

	/* allocate element results space */
	e_values.Allocate(fNumElements, e_out);
	e_values = 0.0;

	/* work arrays */
	dArray2DT nodal_space(fNumElemNodes, n_out);
	dArray2DT nodal_all(fNumElemNodes, n_out);
	dArray2DT coords, disp;
	dArray2DT jump, T;
	dArray2DT matdat;	

	/* ip values */
	LocalArrayT loc_init_coords(LocalArrayT::kInitCoords, fNumElemNodes, fNumSD);
	LocalArrayT loc_disp(LocalArrayT::kDisp, fNumElemNodes, fNumDOF);
	fFEManager.RegisterLocal(loc_init_coords);
	fFEManager.RegisterLocal(loc_disp);
	dArrayT ipmat(n_codes[MaterialData]);
	
	/* set shallow copies */
	double* pall = nodal_space.Pointer();
	coords.Set(fNumElemNodes, n_codes[NodalCoord], pall) ; pall += coords.Length();
	disp.Set(fNumElemNodes, n_codes[NodalDisp], pall)    ; pall += disp.Length();
	jump.Set(fNumElemNodes, n_codes[NodalDispJump], pall); pall += jump.Length();
	T.Set(fNumElemNodes, n_codes[NodalTraction], pall)   ; pall += T.Length();
	matdat.Set(fNumElemNodes, n_codes[MaterialData], pall);

	/* element work arrays */
	dArrayT element_values(e_values.MinorDim());
	pall = element_values.Pointer();
	dArrayT centroid;
	if (e_codes[Centroid])
	{
		centroid.Set(fNumSD, pall); 
		pall += fNumSD;
	}
	double phi_tmp, area;
	double& phi = (e_codes[CohesiveEnergy]) ? *pall++ : phi_tmp;

	Top();
	while (NextElement())
	{
		/* current element */
		ElementCardT& element = CurrentElement();
	
		/* initialize */
	    nodal_space = 0.0;

		/* coordinates for whole element */
		if (n_codes[NodalCoord])
		{
			SetLocalX(loc_init_coords);
			loc_init_coords.ReturnTranspose(coords);
		}
		
		/* displacements for whole element */
		if (n_codes[NodalDisp])
		{
			SetLocalU(loc_disp);
			loc_disp.ReturnTranspose(disp);
		}

		/* compute output */
		if (element.Flag() == kON)
		{
	  		/* surface potential */
			SurfacePotentialT* surfpot = fSurfPots[element.MaterialNumber()];

			/* get current geometry */
			SetLocalX(fLocCurrCoords); //EFFECTIVE_DVA

			/* initialize element values */
			phi = area = 0;
			if (e_codes[Centroid]) centroid = 0.0;

			/* integrate */
			fShapes->TopIP();
			while (fShapes->NextIP())
			{
				/* element integration weight */
				double ip_w = fShapes->Jacobian()*fShapes->IPWeight();

				/* gap */
				const dArrayT& gap = fShapes->InterpolateJumpU(fLocCurrCoords);

				/* coordinate transformation */
				double j = fCurrShapes->Jacobian(fQ);
				fQ.MultTx(gap, fdelta);

				/* gap */				
				if (n_codes[NodalDispJump])
					fShapes->Extrapolate(fdelta, jump);				
				
				/* traction */
				if (n_codes[NodalTraction])
					fShapes->Extrapolate(surfpot->Traction(fdelta), T);
					
				/* material output data */
				if (n_codes[MaterialData])
				{
					surfpot->ComputeOutput(fdelta, ipmat);
					fShapes->Extrapolate(ipmat, matdat);
				}

				/* area-averaged centroid */
				if (e_codes[Centroid])
				{
					/* mass */
					area += ip_w;
				
					/* moment */
					centroid.AddScaled(ip_w, fShapes->IPCoords());
				}
				
				/* cohesive energy */
				if (e_codes[CohesiveEnergy])
				{
					/* surface potential */
					double potential = surfpot->Potential(fdelta);

					/* integrate */
					phi += potential*ip_w;
				}
			}
		}

		/* copy in the cols (in sequence of output) */
		int colcount = 0;
		nodal_all.BlockColumnCopyAt(disp  , colcount); colcount += disp.MinorDim();
		nodal_all.BlockColumnCopyAt(coords, colcount); colcount += coords.MinorDim();
		nodal_all.BlockColumnCopyAt(jump  , colcount); colcount += jump.MinorDim();
		nodal_all.BlockColumnCopyAt(T     , colcount); colcount += T.MinorDim();
		nodal_all.BlockColumnCopyAt(matdat, colcount);

		/* accumulate - extrapolation done from ip's to corners => X nodes */
		fNodes->AssembleAverage(element.NodesX(), nodal_all);
		
		/* element values */
		if (e_codes[Centroid]) centroid /= area;
		
		/* store results */
		e_values.SetRow(CurrElementNumber(), element_values);		
	}

	/* get nodally averaged values */
	fNodes->OutputUsedAverage(n_values);
}

void CSEAnisoT::GenerateOutputLabels(const iArrayT& n_codes, ArrayT<StringT>& n_labels,
	const iArrayT& e_codes, ArrayT<StringT>& e_labels) const
{
	/* inherited */
	CSEBaseT::GenerateOutputLabels(n_codes, n_labels, e_codes, e_labels);

	/* overwrite nodal labels */
	n_labels.Allocate(n_codes.Sum());
	int count = 0;
	if (n_codes[NodalDisp])
	{
		if (fNumDOF > 3) throw eGeneralFail;
		const char* dlabels[3] = {"D_X", "D_Y", "D_Z"};
		for (int i = 0; i < fNumDOF; i++)
			n_labels[count++] = dlabels[i];
	}

	if (n_codes[NodalCoord])
	{
		const char* xlabels[3] = {"x1", "x2", "x3"};
		for (int i = 0; i < fNumSD; i++)
			n_labels[count++] = xlabels[i];
	}

	if (n_codes[NodalDispJump])
	{
		const char* d_2D[2] = {"d_t", "d_n"};
		const char* d_3D[3] = {"d_t1", "d_t2", "d_n"};
		const char** dlabels;
		if (fNumDOF == 2)
			dlabels = d_2D;
		else if (fNumDOF == 3)
			dlabels = d_3D;
		else
			throw eGeneralFail;

		for (int i = 0; i < fNumDOF; i++)
			n_labels[count++] = dlabels[i];
	}

	if (n_codes[NodalTraction])
	{
		const char* t_2D[2] = {"T_t", "T_n"};
		const char* t_3D[3] = {"T_t1", "T_t2", "T_n"};
		const char** tlabels;
		if (fNumDOF == 2)
			tlabels = t_2D;
		else if (fNumDOF == 3)
			tlabels = t_3D;
		else
			throw eGeneralFail;

		for (int i = 0; i < fNumDOF; i++)
			n_labels[count++] = tlabels[i];
	}

	/* material output labels */
	if (n_codes[MaterialData])
	{
		ArrayT<StringT> matlabels;
		fSurfPots[0]->OutputLabels(matlabels);
		for (int i = 0; i < n_codes[MaterialData]; i++)
			n_labels[count++] = matlabels[i];
	}
}

/* write all current element information to the stream */
void CSEAnisoT::CurrElementInfo(ostream& out) const
{
	/* inherited */
	CSEBaseT::CurrElementInfo(out);
	
	/* current element configuration */
	out << " current integration point:" << fCurrShapes->CurrIP() << '\n';
	try
	{
		dMatrixT Qtemp(fQ);
		double j = fCurrShapes->Jacobian(Qtemp);
		out << " surface jacobian = " << j << '\n';
		out << " coordinate transformation:\n";
		out << fQ << '\n';
		out << " gap vector (global frame):\n";
		const dArrayT& delta = fShapes->InterpolateJumpU(fLocCurrCoords);
		out << delta << '\n';
		out << " gap vector (local frame):\n";
		dArrayT delta_temp(fdelta);
		fQ.MultTx(delta, delta_temp);
		out << delta_temp << '\n';	
	}
	
	catch (int error)
	{
		out << " CSEAnisoT::CurrElementInfo: error on surface jacobian\n";
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* operations with pseudo rank 3 (list in j) matrices */
void CSEAnisoT::u_i__Q_ijk(const dArrayT& u, const ArrayT<dMatrixT>& Q,
	dMatrixT& Qu)
{
	for (int i = 0; i < u.Length(); i++)
	{	
		Q[i].MultTx(u, fNEEvec);
		Qu.SetRow(i, fNEEvec);
	}
}

void CSEAnisoT::Q_ijk__u_j(const ArrayT<dMatrixT>& Q, const dArrayT& u,
	dMatrixT& Qu)
{
	if (Q.Length() == 2)
		Qu.SetToCombination(u[0], Q[0], u[1], Q[1]);
	else if (Q.Length() == 3)
		Qu.SetToCombination(u[0], Q[0], u[1], Q[1], u[2], Q[2]);
	else
		throw eGeneralFail;
}
