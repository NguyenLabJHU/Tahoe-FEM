/* $Id: ViscousDragT.cpp,v 1.4 2004-07-15 08:26:14 paklein Exp $ */
#include "ViscousDragT.h"

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "ModelManagerT.h"
#include "ShapeFunctionT.h"
#include "OutputSetT.h"

using namespace Tahoe;

/* constructor */
ViscousDragT::ViscousDragT(const ElementSupportT& support):
	ElementBaseT(support),
	fViscosity(0.0)
{
	SetName("viscous_drag");
ExceptionT::GeneralFail("ViscousDragT::ViscousDragT", "out of date");
}

#if 0
/* class initialization */
void ViscousDragT::Initialize(void)
{
#pragma message("delete me")
	/* inherited */
	ElementBaseT::Initialize();

	/* read parameters */
	ifstreamT& in = ElementSupport().Input();
	in >> fViscosity
	   >> fID;
	   
	/* model manager */
	ModelManagerT& model = ElementSupport().ModelManager();
	const iArray2DT& connects = model.ElementGroup(fID);
	iArrayT nodes_used;
	nodes_used.Union(connects);
	fNodesUsed.Dimension(nodes_used.Length(), 1);
	fNodesUsed.SetColumn(0, nodes_used);
	nodes_used.Free();
	fNodalMass.Dimension(fNodesUsed.Length());
	fDragForce.Dimension(fNodesUsed.Length(), NumDOF());

	/* shape functions */
	GeometryT::CodeT geometry_code = model.ElementGroupGeometry(fID);
	LocalArrayT ref_coords(LocalArrayT::kInitCoords, connects.MinorDim(), NumSD());
	ElementSupport().RegisterCoordinates(ref_coords);
	ShapeFunctionT shape(geometry_code, 1, ref_coords);
	shape.Initialize();

	/* shape function information */
	const double* j0 = shape.IPDets();
	const double* w  = shape.IPWeights();
	const double* Na = shape.IPShapeU(0);

	/* calculate nodal mass */
	fGlobalToLocal.SetOutOfRange(InverseMapT::MinusOne);
	fGlobalToLocal.SetMap(fNodesUsed);
	int nel = connects.MajorDim();
	int nen = connects.MinorDim();
	fNodalMass = 0.0;
	iArrayT elem_nodes;
	for (int i = 0; i < nel; i++)
	{
		/* collect element info */
		connects.RowAlias(i, elem_nodes);
		ref_coords.SetLocal(elem_nodes);
		
		/* set shape functions */ 	
		shape.SetDerivatives();

		/* accumulate over element nodes - 1 ip*/
		for (int j = 0; j < nen; j++) {
			int index = fGlobalToLocal.Map(elem_nodes[j]);
			fNodalMass[index] += (*w)*(*j0)*Na[j]; 
		}
	}

	/* echo parameters */
	ofstreamT& out = ElementSupport().Output();
	out << " Viscosity per unit volume . . . . . . . . . . . = " << fViscosity << '\n';
	out << " Element block ID. . . . . . . . . . . . . . . . = " << fID << '\n';
	out << " Number of nodes used. . . . . . . . . . . . . . = " << fNodesUsed.Length() << '\n';
	if (ElementSupport().PrintInput())
		out << fNodesUsed.wrap(5) << '\n';
}
#endif

/* collecting element group equation numbers */
void ViscousDragT::Equations(AutoArrayT<const iArray2DT*>& eq_1, AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* inherited */
	ElementBaseT::Equations(eq_1, eq_2);

	/* collect equation numbers */
	fEqnos.Dimension(fNodesUsed.Length(), NumDOF());
	iArray2DT connects_temp(fNodesUsed.Length(), 1, fNodesUsed.Pointer());
	Field().SetLocalEqnos(connects_temp, fEqnos);
}

/* accumulate the residual force on the specified node */
void ViscousDragT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
	/* not my field */
	if (&field != &(Field())) return;

	/* local index */
	int index = fGlobalToLocal.Map(node);
	if (index < 0) return;

	double dt = ElementSupport().TimeStep();
	if (fabs(dt) < kSmall || fabs(fViscosity) < kSmall) return;

	/* displacements */
	const dArray2DT& disp = field(0,0);
	const dArray2DT& disp_last = field(-1,0);

	/* drag force */
	double visc = fNodalMass[index]*fViscosity/dt;
	int ndof = NumDOF();
	for (int i = 0; i < ndof; i++)
		force[i] += visc*(disp(index,i), disp_last(index,i));
}

void ViscousDragT::RegisterOutput(void)
{
	/* output labels */
	const char *d_labels[] = {"D_X", "D_Y", "D_Z"};
	const char *f_labels[] = {"F_X", "F_Y", "F_Z"};
	int ndof = NumDOF();
	if (ndof > 3) ExceptionT::GeneralFail("ViscousDragT::RegisterOutput", "ndof > 3: %d", ndof);
	ArrayT<StringT> labels(2*ndof);
	int index = 0;
	for (int i = 0; i < ndof; i++)
		labels[index++] = d_labels[i];
	for (int i = 0; i < ndof; i++)
		labels[index++] = f_labels[i];

	/* create output set */
	OutputSetT output(GeometryT::kPoint, fNodesUsed, labels);
	fOutputID = ElementSupport().RegisterOutput(output);
}

void ViscousDragT::WriteOutput(void) 
{
	/* write incremental viscous dissipation */
	ofstreamT& out = ElementSupport().Output();
	out << "\n ViscousDragT::WriteOutput: incremental dissipation = " << fIncrementalDissipation << endl;

	/* collect nodal values */
	int ndof = NumDOF();
	dArray2DT n_values(fNodesUsed.Length(), 2*ndof), e_values;
	const dArray2DT& disp = (Field())[0];

	/* write in by columns */
	int index = 0;
	int nnd = fNodesUsed.Length();
	for (int i = 0; i < n_values.MinorDim(); i++)
	{
		if (index < ndof) {
			for (int j = 0; j < nnd; j++)
				n_values(j, index) = disp(fNodesUsed[j], index);
		}
		else {
			for (int j = 0; j < nnd; j++)
				n_values(j, index) = fDragForce(j, index-ndof);
		}
		index++;
	}

	/* write output */
	ElementSupport().WriteOutput(fOutputID, n_values, e_values);
}

void ViscousDragT::SendOutput(int) {
	ExceptionT::GeneralFail("ViscousDragT::SendOutput", "not supported");
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* form group contribution to the stiffness matrix */
void ViscousDragT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
#pragma unused(sys_type)

	double dt = ElementSupport().TimeStep();
	if (fabs(dt) < kSmall || fabs(fViscosity) < kSmall) return;

	/* displacements */
	const FieldT& field = Field();
	const dArray2DT& disp = field(0,0);
	const dArray2DT& disp_last = field(-1,0);

	/* compute the nodal forces */
	iArrayT eqnos;
	ElementMatrixT stiffness(NumDOF(), ElementMatrixT::kDiagonal);
	for (int i = 0; i < fNodesUsed.Length(); i++)
	{
		stiffness.Identity(fNodalMass[i]*fViscosity/dt);

		/* assemble */
		fEqnos.RowAlias(i, eqnos);
		ElementSupport().AssembleLHS(field.Group(), stiffness, eqnos);		
	}	
}

/* form group contribution to the residual */
void ViscousDragT::RHSDriver(void)
{
	double dt = ElementSupport().TimeStep();
	if (fabs(dt) < kSmall || fabs(fViscosity) < kSmall) return;

	/* displacements */
	const FieldT& field = Field();
	const dArray2DT& disp = field(0,0);
	const dArray2DT& disp_last = field(-1,0);

	/* compute the nodal forces */
	fIncrementalDissipation = 0.0;
	fDragForce = 0.0;
	dArrayT inc_disp(NumDOF());
	for (int i = 0; i < fNodesUsed.Length(); i++)
	{
		/* incremental displacement */
		int node = fNodesUsed[i];
		inc_disp.DiffOf(disp(node), disp_last(node));
		
		/* viscosity */
		double viscosity = fNodalMass[i]*fViscosity/dt;
		fDragForce.AddToRowScaled(i, -viscosity, inc_disp);
		fIncrementalDissipation += viscosity*dArrayT::Dot(inc_disp, inc_disp);
	}
	
	/* assemble */
	ElementSupport().AssembleRHS(field.Group(), fDragForce, fEqnos);
}
