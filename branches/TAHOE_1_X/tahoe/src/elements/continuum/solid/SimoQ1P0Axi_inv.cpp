/* $Id: SimoQ1P0Axi_inv.cpp,v 1.2 2004-06-26 18:39:04 paklein Exp $ */
#include "SimoQ1P0Axi_inv.h"

#include "ShapeFunctionT.h"
#include "SolidMaterialT.h"
#include "SolidMatListT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"

const double Pi2 = 2.0*acos(-1.0);
const int kRadialDirection = 0; /* x <-> r */

using namespace Tahoe;

/* constructor */
SimoQ1P0Axi_inv::SimoQ1P0Axi_inv(const ElementSupportT& support, const FieldT& field):
	UpdatedLagrangianAxiT(support, field),
	fF_tmp(NumSD()),
	fOutputInit(false),
	fOutputCell(-1)
{

}

/* data initialization */
void SimoQ1P0Axi_inv::Initialize(void)
{
	const char caller[] = "SimoQ1P0Axi_inv::Initialize";

	/* inherited */
	UpdatedLagrangianAxiT::Initialize();

	/* check geometry code and number of element nodes -> Q1 */
	if (GeometryCode() == GeometryT::kQuadrilateral) {
		if (NumElementNodes() != 4) 
			ExceptionT::BadInputValue(caller, "expecting 4 node quad: %d", NumElementNodes());
	}
	else if (GeometryCode() == GeometryT::kHexahedron) {
		if (NumElementNodes() != 8) 
			ExceptionT::BadInputValue(caller, "expecting 8 node hex: %d", NumElementNodes());
	}
	else
		ExceptionT::BadInputValue(caller, "expecting hex or quad geometry: %d", GeometryCode());
	
	/* need to store last deformed element volume */
	fElementVolume.Dimension(NumElements());	
	fElementVolume = 0.0;
	fGamma.Dimension(NumElements());
	fGamma = 0.0;
	fGamma_last.Dimension(NumElements());
	fGamma_last = 0.0;
	
	/* element pressure */
	fPressure.Dimension(NumElements());
	fPressure = 0.0;
	
	/* determinant of the deformation gradient */
	fJacobian.Dimension(NumIP());
	fJacobian = 1.0;
	
	/* dimension work space */
	fMeanGradient.Dimension(NumSD(), NumElementNodes());
	fNEEmat.Dimension(fLHS);
	fdiff_b.Dimension(fGradNa);
	fb_bar.Dimension(fGradNa);
	fb_sig.Dimension(fGradNa);

	/* need to initialize previous Gamma */
	Top();
	while (NextElement())
	{
		/* inherited - computes gradients and standard 
		 * deformation gradients */
		UpdatedLagrangianAxiT::SetGlobalShape();

		/* compute mean of shape function gradients */
		double& V = fElementVolume[CurrElementNumber()]; /* reference volume */
		double& Gamma = fGamma_last[CurrElementNumber()];
		SetMeanGradient(fMeanGradient, V, Gamma);
	}

	/* check cell output */
	int index;
	if (ElementSupport().CommandLineOption("-track_group", index)) {
		const ArrayT<StringT>& argv = ElementSupport().Argv();
		int group = -99;
		group = atoi(argv[index+1]) - 1;
		if (group == ElementSupport().ElementGroupNumber(this))
			if (ElementSupport().CommandLineOption("-track_cell", index))
				fOutputCell = atoi(argv[index+1]) - 1;
	}
}

/* finalize current step - step is solved */
void SimoQ1P0Axi_inv::CloseStep(void)
{
	/* inherited */
	UpdatedLagrangianAxiT::CloseStep();
	
	/* store converged solution */
	fGamma_last = fGamma;
}
	
/* restore last converged state */
GlobalT::RelaxCodeT SimoQ1P0Axi_inv::ResetStep(void)
{
	/* inherited */
	GlobalT::RelaxCodeT relax = UpdatedLagrangianAxiT::ResetStep();
	
	/* store converged solution */
	fGamma = fGamma_last;

	return relax;
}

/* read restart information from stream */
void SimoQ1P0Axi_inv::ReadRestart(istream& in)
{
	/* inherited */
	UpdatedLagrangianAxiT::ReadRestart(in);
	
	/* read restart data */
	in >> fGamma;
	
	/* reset last state */
	fGamma_last = fGamma;
}

/* write restart information from stream */
void SimoQ1P0Axi_inv::WriteRestart(ostream& out) const
{
	/* inherited */
	UpdatedLagrangianAxiT::WriteRestart(out);
	
	/* read restart data */
	out << fGamma << '\n';
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form shape functions and derivatives */
void SimoQ1P0Axi_inv::SetGlobalShape(void)
{
	const char caller[] = "SimoQ1P0Axi_inv::SetGlobalShape";

	/* current element number */
	int elem = CurrElementNumber();

	/* inherited - computes gradients and standard 
	 * deformation gradients */
	UpdatedLagrangianAxiT::SetGlobalShape();

	/* compute mean of shape function gradients */
	double& V = fElementVolume[elem]; /* reference volume */
	double& Gamma = fGamma[elem];
	SetMeanGradient(fMeanGradient, V, Gamma);
	
	/* last Gamma */
	double& Gamma_last = fGamma_last[elem];

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
			/* "replace" dilatation */
			dMatrixT& F = fF_List[i];
			double J = F.Det();
			if (J < 0.0) ExceptionT::BadJacobianDet(caller, "J = %g", J);
			F *= pow(1.0/(Gamma*J), 1.0/3.0);
			
			/* store Jacobian */
			fJacobian[i] = J;
		}

		/* "last" deformation gradient */
		if (needs_F_last)
		{
			/* "replace" dilatation */
			dMatrixT& F = fF_last_List[i];
			double J = F.Det();
			F *= pow(1.0/(Gamma_last*J), 1.0/3.0);
		}
	}
}

/* form the element stiffness matrix */
void SimoQ1P0Axi_inv::FormStiffness(double constK)
{
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* current element info */
	int el = CurrElementNumber();
	double v = fElementVolume[el];
	double p_bar = fPressure[el];

	/* integration */
	const double* Det    = fCurrShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* initialize */
	fStressStiff = 0.0;
	fNEEvec = 0.0;

	int  nsd = NumSD();
	int ndof = NumDOF();
	int nen  = NumElementNodes();
	fCurrShapes->GradNa(fMeanGradient, fb_bar);	
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		int ip = fShapes->CurrIP();
		double r = fRadius_x[ip];
	
		/* scale factor */
		double scale = Pi2*r*constK*(*Det++)*(*Weight++);

		/* collect array of nodal shape functions */
		const double* Na_u = fCurrShapes->IPShapeU();
		fIPShape.Alias(nen, Na_u);
		double* u_r = fNEEvec.Pointer(kRadialDirection);
		for (int a = 0; a < nen; a++) {
			*u_r = *Na_u++;
			u_r += ndof;
		}
	
	/* S T R E S S   S T I F F N E S S */			
		/* compute Cauchy stress */
		const dSymMatrixT& cauchy = fCurrMaterial->s_ij();
		cauchy.ToMatrix(fStressMat);
		fMat2D.Rank2ReduceFrom3D(fStressMat);
		
		/* determinant of modified deformation gradient */
		double J_bar = DeformationGradient().Det();
		
		/* detF correction */
		double J_correction = J_bar/fJacobian[CurrIP()];
		double p = J_correction*cauchy.Trace()/3.0;

		/* get shape function gradients matrix */
		fCurrShapes->GradNa(fGradNa);
		fb_sig.MultAB(fMat2D, fGradNa);

		/* integration constants */		
		fMat2D *= scale*J_correction;
	
		/* using the stress symmetry */
		fStressStiff.MultQTBQ(fGradNa, fMat2D, format, dMatrixT::kAccumulate);

		/* contribution from out-of-plane stress */
		fLHS.Outer(fNEEvec, fNEEvec, scale*J_correction*fStressMat(2,2)/(r*r), dMatrixT::kAccumulate);

	/* M A T E R I A L   S T I F F N E S S */									
		/* strain displacement matrix */
		Set_B_bar_axi(fIPShape, fCurrShapes->Derivatives_U(), fMeanGradient, r, fB);

		/* get D matrix */
		fD.Rank4ReduceFrom3D(fCurrMaterial->c_ijkl());
		fD *= scale;
						
		/* accumulate */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);

		/* add axisymmetric contribution to b */
		Na_u = fCurrShapes->IPShapeU();
		double* b_r = fGradNa.Pointer(kRadialDirection);
		for (int a = 0; a < nen; a++) {
			*b_r += (*Na_u++)/r;
			b_r += nsd;
		}

		/* $div div$ term */	
		fNEEmat.Outer(fGradNa, fGradNa);
		fLHS.AddScaled(p_bar*scale, fNEEmat);

		fdiff_b.DiffOf(fGradNa, fb_bar);
		fNEEmat.Outer(fdiff_b, fdiff_b);
		fLHS.AddScaled(scale*2.0*p/3.0, fNEEmat);
		
		fNEEmat.Outer(fb_sig, fdiff_b);
		fNEEmat.Symmetrize();
		fLHS.AddScaled(-J_correction*scale*4.0/3.0, fNEEmat);

		bSp_bRq_to_KSqRp(fGradNa, fNEEmat);
		fLHS.AddScaled(scale*(p - p_bar), fNEEmat);
	}
						
	/* stress stiffness into fLHS */
	fLHS.Expand(fStressStiff, NumDOF(), dMatrixT::kAccumulate);
	
	/* $\bar{div}\bar{div}$ term */
	fNEEmat.Outer(fb_bar, fb_bar);
	fLHS.AddScaled(-p_bar*v, fNEEmat);
}

/* calculate the internal force contribution ("-k*d") */
void SimoQ1P0Axi_inv::FormKd(double constK)
{
	const double* Det    = fCurrShapes->IPDets();
	const double* Det0   = fShapes->IPDets();
	const double* Weight = fCurrShapes->IPWeights();

	/* collect incremental heat */
	bool need_heat = fElementHeat.Length() == fShapes->NumIP();

	/* current element number */
	int elem = CurrElementNumber();

	/* constant pressure */
	double& p_bar = fPressure[elem];
	p_bar = 0.0;

	bool hit_cell = false;
	int nen = NumElementNodes();
	fCurrShapes->TopIP();
	while ( fCurrShapes->NextIP() )
	{
		int ip = fShapes->CurrIP();
		double r = fRadius_x[ip];
		double R = fRadius_X[ip];
	
		/* collect array of nodal shape functions */
		fIPShape.Alias(nen, fShapes->IPShapeU());
	
		/* strain displacement matrix */
		Set_B_bar_axi(fIPShape, fCurrShapes->Derivatives_U(), fMeanGradient, r, fB);

		/* translate Cauchy stress to axisymmetric */
		const dSymMatrixT& cauchy = fCurrMaterial->s_ij();
		fStress2D_axi.ReduceFrom3D(cauchy);
		
		/* B^T * Cauchy stress */
		fB.MultTx(fStress2D_axi, fNEEvec);
		
		/* determinant of modified deformation gradient */
		double J_bar = DeformationGradient().Det();
		
		/* detF correction */
		double J_correction = J_bar/fJacobian[CurrIP()];
		
		/* integrate pressure */
		p_bar += Pi2*R*(*Weight)*(*Det0)*J_correction*cauchy.Trace()/3.0;
		
		/* accumulate */
		fRHS.AddScaled(constK*Pi2*r*(*Weight)*(*Det)*J_correction, fNEEvec);

		/* incremental heat generation */
		if (need_heat) 
			fElementHeat[fShapes->CurrIP()] += fCurrMaterial->IncrementalHeat();

		/* debugging output */
		int output_element = fOutputCell;
		if (CurrElementNumber() == output_element) {

			/* collect nodal velocities */
			if (CurrIP() == 0) 
				SetLocalU(fLocVel);

			/* step information */
			int step_number = ElementSupport().StepNumber();
			double time = ElementSupport().Time();
		
			/* acoustic wave speeds */
			dArrayT normal(3), speeds(3);
			normal[0] = 1.0;
			normal[1] = 0.0;
			normal[2] = 0.0;
			fCurrMaterial->WaveSpeeds(normal, speeds);

			/* neighborhood nodes */
			const iArrayT& nodes_u = CurrentElement().NodesU();

			/* get matrix of shape function gradients */
			fShapes->GradNa(fGradNa);

			/* include out-of-plane influence */
			const double* NaU = fShapes->IPShapeU();
			double r = fRadius_x[CurrIP()];
			for (int i = 0; i < nodes_u.Length(); i++)
				fGradNa(0,i) += (*NaU++)/r;
			
			/* file path */
			StringT path;
			path.FilePath(ElementSupport().Input().filename());
			
			/* write info for neighborhood nodes */
			for (int i = 0; i < nodes_u.Length(); i++) {

				/* file name */
				StringT node_file;
				node_file.Append("cell", output_element + 1);
				node_file.Append(".ip", CurrIP() + 1);
				node_file.Append(".nd", nodes_u[i] + 1);
				node_file.Append(".dat");
				node_file.Prepend(path);
				
				/* (re-)open stream */
				ofstreamT out;
				if (fOutputInit)
					out.open_append(node_file);
				else {
					out.open(node_file);

					/* Tecplot style data headers */				
					out << "VARIABLES = \"step\" \"time\" \"J\" \"J_bar\" \"Na_r\" \"Na_z\" \"v_r\" \"v_z\" \"c_d\" \"c_s1\" \"c_s2\"" << endl;
				}
					
				/* write output */
				int d_width = OutputWidth(out, &time);
				out << setw(kIntWidth) << step_number
				    << setw(d_width) << time
				    << setw(d_width) << fJacobian[CurrIP()] 
				    << setw(d_width) << J_bar
				    << setw(d_width) << fGradNa(0,i)
				    << setw(d_width) << fGradNa(1,i)
				    << setw(d_width) << fLocVel(i,0)
				    << setw(d_width) << fLocVel(i,1)
				    << speeds.no_wrap() << '\n';
				    
				/* close stream */
				out.close();
			}

			/* set flag */
			hit_cell = true;
		}

		/* next integration points */
		Det++; Det0++; Weight++;
	}

	/* append to results files */
	if (hit_cell) fOutputInit = true;	
	
	/* volume averaged */
	p_bar /= fElementVolume[CurrElementNumber()];
}

/* read materials data */
void SimoQ1P0Axi_inv::ReadMaterialData(ifstreamT& in)
{
	/* inherited */
	UpdatedLagrangianAxiT::ReadMaterialData(in);

	/* make sure 2D materials are plane strain */
	if (StructuralMaterialList().HasPlaneStress()) 
		ExceptionT::BadInputValue("SimoQ1P0Axi_inv::ReadMaterialData", "2D materials must be plane strain");
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* compute mean shape function gradient */
void SimoQ1P0Axi_inv::SetMeanGradient(dArray2DT& mean_gradient, double& V, double& Gamma) const
{
	/* assume same integration rule defined for current and references
	 * shape functions */
	int nip = NumIP();
	const double*   det = fCurrShapes->IPDets();
	const double* det_0 = fShapes->IPDets();
	const double*     w = fShapes->IPWeights();

	/* V and the inverse dilation */
	V = 0.0;
	Gamma = 0.0;
	for (int i = 0; i < nip; i++) {

		double R = fRadius_X[i];
		double r = fRadius_x[i];

		double dV = Pi2*R*w[i]*det_0[i];
		double dv = Pi2*r*w[i]*det[i];
	
		V += dV;
		double J = dv/dV;
		Gamma += dV/J;
	}
	Gamma /= V;

	/* initialize */
	mean_gradient = 0.0;			

	/* integrate */
	int nen = mean_gradient.MinorDim();	
	for (int i = 0; i < nip; i++) {
	
		double R = fRadius_X[i];
		double r = fRadius_x[i];		

		double dV = Pi2*R*w[i]*det_0[i];
		double dv = Pi2*r*w[i]*det[i];

		double J = dv/dV;
		double dV_by_JVGamma = dV/(V*Gamma*J);
		
		mean_gradient.AddScaled(dV_by_JVGamma, fCurrShapes->Derivatives_U(i));
		
		/* contribution from out-of-plane component */
		double* mean_r = mean_gradient(kRadialDirection);
		const double* pNaU = fCurrShapes->IPShapeU(i);
		for (int a = 0; a < nen; a++)
			*mean_r++ += dV_by_JVGamma*(*pNaU++)/r;		
	}
}

void SimoQ1P0Axi_inv::bSp_bRq_to_KSqRp(const dMatrixT& b, dMatrixT& K) const
{
#if __option(extended_errorcheck)
	/* dimension check */
	if (b.Length() != K.Rows() ||
	    K.Rows() != K.Cols()) ExceptionT::SizeMismatch("SimoQ1P0Axi_inv::bSp_bRq_to_KSqRp");
#endif

	int dim = K.Rows();
	int sub_dim = b.Rows();
	int S = 0;
	int p = 0;
	for (int i = 0; i < dim; i++)
	{
		int R = 0;
		int q = 0;
		for (int j = 0; j < dim; j++)
		{
			K(i,j) = b(q,S)*b(p,R);
		
			q++;
			if (q == sub_dim) {
				R++;
				q = 0;
			}
		}
		p++;
		if (p == sub_dim) {
			S++;
			p = 0;
		}
	}	
}
