
#include "ifstreamT.h"
#include "ModelT.h"

//---------------------------------------------------------

ModelT::ModelT ( StringT &model_filename) 
{
	int i;
	model_name.Dimension (32); 
	solver_filename.Dimension (32);  					
	node_filename.Dimension (32);  					

	ifstreamT model (model_filename); 
	model.set_marker ( '#' );

	model >> model_name;
	model >> n_sd;
	model >> constraint;
	model >> initial_time;
	model >> delta_t;
	model >> time_int_scheme;
	model >> solver_type;
	model >> solver_filename;
	model >> node_filename;

	model >> input_file_type;
	model >> n_es;
	
  	elmt_set_filenames.Dimension(n_es);
  	elmt_set_surface_filenames.Dimension(n_es);
	for (i=0; i<n_es; i++)	{
		elmt_set_filenames[i].Dimension (32); 
		elmt_set_surface_filenames[i].Dimension (32); 
		model >> elmt_set_filenames[i];
		model >> elmt_set_surface_filenames[i];
	}

	model >> n_lc;
  	load_curve_filenames.Dimension(n_lc);
	for (i=0; i<n_lc; i++)	{
		load_curve_filenames[i].Dimension (32);
		model >> load_curve_filenames[i];
	}

	model >> n_ms;
  	material_set_filenames.Dimension(n_ms);
	for (i=0; i<n_ms; i++)	{
		material_set_filenames[i].Dimension (32);
		model >> material_set_filenames[i];
	}

	model >> n_fld;
 	Field.Dimension(n_fld);
	for (i=0; i<n_fld; i++)	{
		Field[i].name.Dimension (32);  	 											
		Field[i].point_dirichlet_filename.Dimension (32);  	 		
		Field[i].point_neumann_filename.Dimension (32);  		 		
		Field[i].line_dirichlet_filename.Dimension (32);  	 			
		Field[i].line_neumann_filename.Dimension (32);  		 			
		Field[i].surface_dirichlet_filename.Dimension (32);   		
		Field[i].surface_neumann_filename.Dimension (32);  	 		
		Field[i].surface_neumann_filename.Dimension (32);  	 		
					
		model >> Field[i].name;
		model >> Field[i].n_dof;
		model >> Field[i].body_force;
		model >> Field[i].nl_neumann;
		model >> Field[i].point_dirichlet_filename;
		model >> Field[i].point_neumann_filename;
		model >> Field[i].line_dirichlet_filename;
		model >> Field[i].line_neumann_filename;
		model >> Field[i].surface_dirichlet_filename;
		model >> Field[i].surface_neumann_filename;
	}

}

//---------------------------------------------------------

void ModelT::Print ( void )
{
	int i;

	cout << "\n\n\n\n Contents of ModelT follows:\n\n"; 

	cout << "model name: " 			<< model_name 			<<"\n";
	cout << "n_sd: " 						<< n_sd 						<<"\n";
	cout << "constraint: "			<< constraint 			<<"\n";
	cout << "initial_time: " 		<< initial_time			<<"\n";
	cout << "time_step: " 			<< delta_t 					<<"\n";
	cout << "time_int_scheme: " << time_int_scheme 	<<"\n";
	cout << "solver_type: " 		<< solver_type 			<<"\n";
	cout << "solver_filename: " << solver_filename 	<<"\n";
	cout << "node_filename: " 	<< node_filename 		<<"\n";
	cout << "input_file_type: " << input_file_type 	<<"\n";

	cout << "n_es: "						<< n_es <<"\n";
	for (i=0; i<n_es; i++) {	
		cout << "elmt_set_filenames "<<i<<": " 	<< elmt_set_filenames[i] <<"\n";
		cout << "elmt_set_surface_filenames "<<i<<": " 	<< elmt_set_surface_filenames[i] <<"\n";
	}

	cout << "n_lc: "						<< n_lc <<"\n";
	for (i=0; i<n_lc; i++)	
		cout << "load_curve_filename "<<i<<": " 	<< load_curve_filenames[i] <<"\n";

	cout << "n_ms: "						<< n_ms <<"\n";
	for (i=0; i<n_ms; i++)	
		cout << "material_set_filename "<<i<<": " << material_set_filenames[i] <<"\n";

	cout << "n_fld: " 					<< n_fld <<"\n";
	for (i=0; i<n_fld; i++)	{
		cout << "Field #" <<i<< ": name = "												<< Field[i].name 				<<"\n";
		cout << "Field #" <<i<< ": n_dof = "											<< Field[i].n_dof 			<<"\n";
		cout << "Field #" <<i<< ": body_force = "									<< Field[i].body_force 	<<"\n";
		cout << "Field #" <<i<< ": nl_neumann = "									<< Field[i].nl_neumann 	<<"\n";
		cout << "Field #" <<i<< ": point_dirichlet_filename = "		<< Field[i].point_dirichlet_filename 		<<"\n";
		cout << "Field #" <<i<< ": point_neumann_filename = "			<< Field[i].point_neumann_filename 			<<"\n";
		cout << "Field #" <<i<< ": line_dirichlet_filename = "		<< Field[i].line_dirichlet_filename 		<<"\n";
		cout << "Field #" <<i<< ": line_neumann_filename = "			<< Field[i].line_neumann_filename 			<<"\n";
		cout << "Field #" <<i<< ": surface_dirichlet_filename = "	<< Field[i].surface_dirichlet_filename 	<<"\n";
		cout << "Field #" <<i<< ": surface_neumann_filename = "		<< Field[i].surface_neumann_filename 		<<"\n";
	}
	cout <<"\n";
}
				
