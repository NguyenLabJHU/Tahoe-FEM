/* $Id: ModelT.h,v 1.2 2003-05-05 00:58:08 paklein Exp $ */
#ifndef _MODELT_H_  
#define _MODELT_H_ 

#include "iArrayT.h"
#include "nArrayT.h"
#include "sArrayT.h"
#include "StringT.h"

using namespace Tahoe;

//--------------------------------------------------

class Field_T 
{
	public:
	
		StringT name;
		int 		n_dof;
		bool 		body_force;      // is a body force present (yes,no)
		bool 		nl_neumann;   	 // are neumann b.c. non-linear (i.e. functions of displacement) 
		StringT point_dirichlet_filename;		
		StringT point_neumann_filename;		
		StringT line_dirichlet_filename;		
		StringT line_neumann_filename;		
		StringT surface_dirichlet_filename;		
		StringT surface_neumann_filename;		

};

//--------------------------------------------------

class ModelT // Geometry, Material, Constitutive, and Interpolation Data
{
	public:

		ModelT() { }
		ModelT(StringT &model_filename);

		void Print (void);

		int input_file_type; 	// ex. { 0,1 } --> { monterey, tahoe }
		int constraint;  			// ex. plane strain
		int solver_type;			// ex. Newton
		int n_fld;						// number of unknown fields		
		int n_sd;							// number of spatial dimensions 
		int n_es;							// number of element sets
		int n_lc;							// number of load curves
		int n_ms;							// number of material sets 
		double initial_time;  // to 
		double delta_t;				// time step  
		int time_int_scheme; 	// Temporal integration scheme (Back Euler, Crank-Nic, ...)
		StringT model_name;		
		StringT node_filename;		
		StringT solver_filename;	
		sArrayT elmt_set_filenames; 				// element set connectivity file names
		sArrayT elmt_set_surface_filenames; // element set surface file names
		sArrayT load_curve_filenames; 			// load curve file names
		sArrayT material_set_filenames; 		// load curve file names

		nArrayT <Field_T> Field;					// Contains field specific data 
};

#endif 
