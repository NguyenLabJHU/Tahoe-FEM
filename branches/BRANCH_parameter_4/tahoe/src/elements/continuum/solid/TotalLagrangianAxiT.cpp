/* $Id: TotalLagrangianAxiT.cpp,v 1.3 2004-06-26 18:35:31 paklein Exp $ */
#include "TotalLagrangianAxiT.h"

#include "ShapeFunctionT.h"
#include "SolidMaterialT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"

const double Pi2 = 2.0*acos(-1.0);
const int kRadialDirection = 0; /* x <-> r */

using namespace Tahoe;

/* constructor */
TotalLagrangianAxiT::TotalLagrangianAxiT(const ElementSupportT& support, const FieldT& field):
	FiniteStrainAxiT(support, field),
	fStressMat(3),
	fTempMat1(3),
	fTempMat2(3),
	fOutputInit(false),
	fOutputCell(-1)
{

}

/* data initialization */
void TotalLagrangianAxiT::Initialize(void)
{
	/* inherited */
	FiniteStrainAxiT::Initialize();

	/* dimension */
	fGradNa.Dimension(NumSD(), NumElementNodes());
	fStressStiff.Dimension(NumElementNodes());
	fTemp2.Dimension(NumElementNodes()*NumDOF());

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

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form the element stiffness matrix */
void TotalLagrangianAxiT::FormStiffness(double constK)
{		
	/* matrix format */
	dMatrixT::SymmetryFlagT format =
		(fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
		dMatrixT::kWhole :
		dMatrixT::kUpperOnly;

	/* integration */
	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	/* initialize */
	fStressStiff = 0.0;
	fNEEvec = 0.0;
	
	int ndof = NumSD();
	int nun = fLocDisp.NumberOfNodes();
	fShapes->TopIP();
	while ( fShapes->NextIP() )
	{
		int ip = fShapes->CurrIP();
		double R = fRadius_X[ip];
		double r = fRadius_x[ip];

		/* collect array of nodal shape functions */
		const double* Na_u = fShapes->IPShapeU();
		fIPShape.Alias(nun, Na_u);
		double* u_r = fNEEvec.Pointer(kRadialDirection);
		for (int a = 0; a < nun; a++) {
			*u_r = *Na_u++;
			u_r += ndof;
		}
			
	/* S T R E S S   S T I F F N E S S */
				
		/* Cauchy stress (and set deformation gradient) */
		(fCurrMaterial->s_ij()).ToMatrix(fStressMat);

		/* F */	
		fTempMat1 = DeformationGradient();
		double F_33 = fTempMat1(2,2);

		/* compute F^-1 in 2D */
		fMat2D.Rank2ReduceFrom3D(fTempMat1);
		double J = fMat2D.Det()*F_33;
		fMat2D.Inverse();
		fTempMat1.Rank2ExpandFrom2D(fMat2D);
		fTempMat1(2,2) = 1.0/F_33;

		/* chain rule shape function derivatives */
		fShapes->TransformDerivatives(fMat2D, fDNa_x);

		/* get shape function gradients matrix */
		fShapes->GradNa(fDNa_x, fGradNa);

		/* scale factor */
		double scale = Pi2*R*constK*(*Det++)*(*Weight++)*J;

		/* integration constants */		
		fStressMat *= scale;
	
		/* using the stress symmetry */
		fMat2D.Rank2ReduceFrom3D(fStressMat);
		fStressStiff.MultQTBQ(fGradNa, fMat2D, format,
			dMatrixT::kAccumulate);

		/* contribution from out-of-plane stress */
		fLHS.Outer(fNEEvec, fNEEvec, fStressMat(2,2)/(r*r), dMatrixT::kAccumulate);
			
	/* M A T E R I A L   S T I F F N E S S */

		/* strain displacement matrix */
		Set_B_axi(fIPShape, fDNa_x, r, fB);

		/* get D matrix */
		fD.Rank4ReduceFrom3D(fCurrMaterial->c_ijkl());
		fD *= scale;
		
		
						
		/* accumulate */
		fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);	
	}
						
	/* stress stiffness into fLHS */
	fLHS.Expand(fStressStiff, NumDOF(), dMatrixT::kAccumulate);
}

/* calculate the internal force contribution ("-k*d") */
void TotalLagrangianAxiT::FormKd(double constK)
{
	/* matrix alias to fTemp */
	dMatrixT fWP(NumSD(), fStressStiff.Rows(), fNEEvec.Pointer());

	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();

	bool hit_cell = false;
	int ndof = NumSD();
	int nun  = fLocDisp.NumberOfNodes();
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
		int ip = fShapes->CurrIP();
		double R = fRadius_X[ip];
	
		/* get Cauchy stress */
		(fCurrMaterial->s_ij()).ToMatrix(fTempMat1);

		/* F */
		fTempMat2 = DeformationGradient();
		double F_33 = fTempMat2(2,2);

		/* compute F^-1 in 2D */
		fMat2D.Rank2ReduceFrom3D(fTempMat2);
		double J = fMat2D.Det()*F_33;
		if (J <= 0.0) ExceptionT::BadJacobianDet("TotalLagrangianAxiT::FormKd");
		fMat2D.Inverse();
		fTempMat2.Rank2ExpandFrom2D(fMat2D);
		fTempMat2(2,2) = 1.0/F_33;

		/* compute PK1/J */
		fStressMat.MultABT(fTempMat1, fTempMat2);
		fMat2D.Rank2ReduceFrom3D(fStressMat);

		/* get matrix of shape function gradients */
		fShapes->GradNa(fGradNa);

		/* Wi,J PiJ */
		fWP.MultAB(fMat2D, fGradNa);

		/* accumulate */
		double scale = Pi2*R*J*constK*(*Weight++)*(*Det++);
		fRHS.AddScaled(scale, fNEEvec);
		
		/* contribution from out-of-plane component: x <-> r */
		scale *= fStressMat(2,2)/R;
		const double* NaU = fShapes->IPShapeU();
		double* pRHS = fRHS.Pointer(kRadialDirection);
		for (int a = 0; a < nun; a++) {
			*pRHS += scale*(*NaU++);
			pRHS += ndof;
		}

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

			/* include out-of-plane influence */
			const double* NaU = fShapes->IPShapeU();
			double R = fRadius_X[CurrIP()];
			for (int i = 0; i < nodes_u.Length(); i++)
				fGradNa(0,i) += (*NaU++)/R;

			/* transform shape function derivatives */
			dMatrixT gradNa(NumSD(), nodes_u.Length());
			const dMatrixT& F_3D = DeformationGradient();
			fMat2D.Rank2ReduceFrom3D(F_3D);
			fMat2D.Inverse();
			gradNa.MultATB(fMat2D, fGradNa);
			
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
					out << "VARIABLES = \"step\" \"time\" \"J\" \"Na_r\" \"Na_z\" \"v_r\" \"v_z\" \"c_d\" \"c_s1\" \"c_s2\"" << endl;
				}					
					
				/* write output */
				int d_width = OutputWidth(out, &time);
				out << setw(kIntWidth) << step_number
				    << setw(d_width) << time
				    << setw(d_width) << J
				    << setw(d_width) << gradNa(0,i)
				    << setw(d_width) << gradNa(1,i)
				    << setw(d_width) << fLocVel(i,0)
				    << setw(d_width) << fLocVel(i,1)
				    << speeds.no_wrap() << '\n';
				    
				/* close stream */
				out.close();
			}

			/* set flag */
			hit_cell = true;
		}
	}	
	
	/* append to results files */
	if (hit_cell) fOutputInit = true;
}
