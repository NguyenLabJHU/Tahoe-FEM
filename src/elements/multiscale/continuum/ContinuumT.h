
#ifndef _CONTINUUMT_H_  
#define _CONTINUUMT_H_ 

#include "FEA.h"
#include "ModelT.h"

#include "ArrayT.h"
#include "nArrayT.h"
#include "nArray2DT.h"
#include "iArray2DT.h"
#include "StringT.h"

//--------------------------------------------------

class Element_SetT 
{
	public:	// NOTE: The notion of a global elmt number is unneccesary

		Element_SetT ( void ) { }
	 	~Element_SetT ( void ) { }

		StringT name;                 	// bottom flange 
		int 		element_type;		// Quad, Tri, Hex ...
		iArrayT constitutive_model;		// StV_Kich, Infin_Strain, BCJ ... (1 per field)
		iArrayT material_set;			// BCJ_Matl, Iso_Matl, .. (See enums in VMF_MaterialT) 
		iArrayT material_set_item;		// (e.x. the third instance in the defs of BCJ_Matl) 
		double geometric_parameter;		// Area of rod, thickness of quad, etc
		iArray2DT IEN;					// Nodal connectivity of each element (See Hughes)
		FEA_ShapeFunctionT Shapes;		// Basis functions for Finite Element Analysis
		iArray2DT Boundary_Surfaces; // See below: (of dimension n_ews x (1+n_es))
		// Boolean flags identifying sides of an element that comprise the body's boundary
		//MD_Particle_Basis  picBasis;	// Basis functions for Molecular Dynamics Analysis

		int n_el, n_en, n_ip, n_esf, n_sf; 	
		
		// Number of: elmnts, elmnt nodes, integ. pts, elmnts w/ surfaces, elmt surfaces (6 for hex) 
		// 						surfaces (total this group)
};

//--------------------------------------------------

class ContinuumT // Geometry, Material, Constitutive, and Interpolation Data
{
	public:

		enum FEA_CodeT 	{ kMonterey, kTahoe };
		enum E_TypeT 	{ kQuad, kHex, kTri, kTet, kRod, kPlate, kBeam };

		ContinuumT () { }
		ContinuumT (ModelT &Manager);
	 	~ContinuumT () { } 

		void Read_Tahoe_Node_Data 		( StringT &node_filename );
		void Read_Tahoe_Element_Data 	( StringT &elmt_set_filename,int e_set );
	  	void Read_Surface_Data 			( StringT &surface_file_name,int e_set );
		void Print 						(int input_file_type);

		dMatrixT Xo; 	// Node position at time t=0 				(node,sd)
		dMatrixT X; 	// Node position last time step			(node,sd)
		dMatrixT x; 	// Node position current time step	(node,sd)

		ArrayT < Element_SetT > Element_Set;

		int n_es, n_el, n_np, n_sd, n_fld; 
};

#endif 
