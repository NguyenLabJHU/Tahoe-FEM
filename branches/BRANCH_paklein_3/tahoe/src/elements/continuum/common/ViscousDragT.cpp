/* $Id: ViscousDragT.cpp,v 1.1.2.1 2003-09-10 13:35:27 paklein Exp $ */
#include "ViscousDragT.h"
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "ModelManagerT.h"
#include "ShapeFunctionT.h"

using namespace Tahoe;

/* constructor */
ViscousDragT::ViscousDragT(const ElementSupportT& support, const FieldT& field):
	ElementBaseT(support, field),
	fViscosity(0)
{

}

/* class initialization */
void ViscousDragT::Initialize(void)
{
	/* inherited */
	ElementBaseT::Initialize();

	/* read parameters */
	ifstreamT& in = ElementSupport().Input();
	in >> fViscosity
	   >> fID;
	   
	/* model manager */
	ModelManagerT& model = ElementSupport().Model();
	const iArray2DT& connects = model.ElementGroup(fID);
	fNodesUsed.Union(connects);
	fNodalMass.Dimension(fNodesUsed.Length());
	fDragForce.Dimension(fNodesUsed.Length(), NumDOF());

	/* shape functions */
	GeometryT::CodeT geometry_code = model.ElementGroupGeometry(fID);
	LocalArrayT ref_coords(LocalArrayT::kInitCoords, connects.MinorDim(), NumSD());
	ElementSupport().RegisterCoordinates(ref_coords);
	ShapeFunctionT shape(geometry_code, 1, ref_coords);
	shape.Initialize();

	/* calculate nodal mass */
	int nel = connects.MajorDim();
	int nen = connects.MinorDim();
	const double* j0 = shape.IPDets();
	const double* w  = shape.IPWeights();
	fNodalMass = 0.0;
	iArrayT elem_nodes;
	for (int i = 0; i < nel; i++)
	{
		/* collect element info */
		connects.RowAlias(i, elem_nodes);
		ref_coords.SetLocal(elem_nodes);
		
		/* set shape functions */ 	
		shape.SetDerivatives();

		/* accumulate over element nodes */
		for (int j = 0; j < nen; j++)
	}

	/* echo parameters */
	ofstreamT& out = ElementSupport().Output();
	out << " Viscosity per unit volume . . . . . . . . . . . = " << fViscosity << '\n';
	out << " Element block ID. . . . . . . . . . . . . . . . = " << fID << '\n';
	out << " Number of nodes used. . . . . . . . . . . . . . = " << fNodesUsed.Length() << '\n';
	if (ElementSupport().PrintInput())
		out << fNodesUsed.wrap(5) << '\n';
}

/* collecting element group equation numbers */
void Equations(AutoArrayT<const iArray2DT*>& eq_1, AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
	/* inherited */
	ElementBaseT::Equations(eq_1, eq_2);

	/* collect equation numbers */
	fEqnos.Dimension(fNodesUsed.Length(), NumDOF());
	Field().SetLocalEqnos(ElementSupport().Model().ElementGroup(fID), fEqnos);
}
