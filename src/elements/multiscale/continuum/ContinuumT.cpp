
#include "ifstreamT.h"
#include "ContinuumT.h"
#include "StringT.h"
#include "iArray2DT.h"

//-----------------------------------------------------------------

ContinuumT::ContinuumT ( ModelT &Manager) 
{
	
	if (Manager.input_file_type == kMonterey) {

		//-- Read node data
					
		ifstreamT nodes (Manager.node_filename);
		nodes.set_marker  ( '#' );
		nodes >> n_np;
		Xo.Dimension (n_np, Manager.n_sd);
		nodes >> Xo;

		//-- Read element data
		
  	n_sd  = Manager.n_sd;	
  	n_es  = Manager.n_es;	
  	n_fld = Manager.n_fld;	
		Element_Set.Dimension(n_es);

		for (int E=0; E<n_es; E++) {

			ifstreamT el_set (Manager.elmt_set_filenames[E]); 
			el_set.set_marker  ( '#' );
			Element_Set[E].name.Dimension (32); 
			el_set >> Element_Set[E].name;
			el_set >> Element_Set[E].element_type; 
			el_set >> Element_Set[E].n_ip; 

			Element_Set[E].constitutive_model.Dimension (n_fld);
			Element_Set[E].material_set.Dimension 			(n_fld);
			Element_Set[E].material_set_item.Dimension	(n_fld);
			for (int f=0; f<n_fld; f++) {
				el_set >> Element_Set[E].constitutive_model[f]; 
				el_set >> Element_Set[E].material_set[f]; 
				el_set >> Element_Set[E].material_set_item[f]; 
				Element_Set[E].constitutive_model[f]--;
				Element_Set[E].material_set[f]--;		// Reduce all numbers to account for 0 
				Element_Set[E].material_set_item[f]--;
			}
			el_set >> Element_Set[E].geometric_parameter;  
			el_set >> Element_Set[E].n_el; 
			el_set >> Element_Set[E].n_en; 

			Element_Set[E].IEN.Dimension ( Element_Set[E].n_el, Element_Set[E].n_en );
			el_set >> Element_Set[E].IEN;

			Element_Set[E].IEN -= 1;  // Reduce all node numbers upon io to account for arrays starting at 0
		}

	} //-- End Monterey 


	if (Manager.input_file_type == kTahoe) { // <--- Single field, Single elmt group only

  	n_sd  = Manager.n_sd;	
		Read_Tahoe_Node_Data 		( Manager.node_filename );
		Read_Tahoe_Element_Data ( Manager.elmt_set_filenames[0],0 );
		Read_Surface_Data 			( Manager.elmt_set_surface_filenames[0],0 );

	} //-- End Tahoe


}	

//-----------------------------------------------------------------

void ContinuumT::Read_Tahoe_Node_Data ( StringT &node_filename )
{
		//-- Node Point DataA
		
		int i,j;
		double dwaste;
		int iwaste;

		ifstreamT nodes ( node_filename );
		nodes.set_marker  ( '#' );
		nodes >> n_np; 
		nodes >> iwaste; 
		Xo.Dimension (n_np, n_sd);

		for (i=0; i<n_np; i++) {
			nodes >> dwaste;
			for (j=0; j<n_sd; j++)
				nodes >> Xo(i,j);
		}

}

//-----------------------------------------------------------------

void ContinuumT::Read_Tahoe_Element_Data ( StringT &elmt_set_filename,int e_set )
{
		//-- Element Nodal Connectivity Data
		
		int i,j;
		double dwaste;

		n_es = 1;
		Element_Set.Dimension(n_es);
		ifstreamT el_set (elmt_set_filename); 
		el_set.set_marker  ( '#' );
		el_set >> Element_Set[e_set].n_el;										
		el_set >> Element_Set[e_set].n_en;				

		if ( Element_Set[e_set].n_en == 4)
						Element_Set[e_set].element_type = kQuad;

		if ( Element_Set[e_set].n_en == 8)
						Element_Set[e_set].element_type = kHex;

		Element_Set[e_set].IEN.Dimension ( Element_Set[e_set].n_el, Element_Set[e_set].n_en );

		for (i=0; i<Element_Set[e_set].n_el; i++) {
			el_set >> dwaste;
			for (j=0; j<Element_Set[e_set].n_en; j++)
				el_set >> Element_Set[e_set].IEN(i,j);
		}

		Element_Set[e_set].IEN -= 1;  // Reduce all node numbers upon io to account for arrays starting at 0

}

//-----------------------------------------------------------------

void ContinuumT::Read_Surface_Data ( StringT &surface_file_name,int e_set )
{

		//-- Element Exterior Surface Connectivity Data 
		
		if ( surface_file_name != "N/A" ) {
		
			int e,s, n_ews, n_esf;  // Number of: elmts w/ surfaces, element surfaces (6 for hex, 4 for tet.)
			ifstreamT surf_set ( surface_file_name ); 
			surf_set.set_marker  ( '#' );
			surf_set >> n_ews; 
			surf_set >> n_esf; 
			Element_Set[e_set].n_esf = n_esf;
			Element_Set[e_set].Boundary_Surfaces.Dimension (Element_Set[e_set].n_el,n_esf);
			Element_Set[e_set].Boundary_Surfaces = 0; 
			Element_Set[e_set].n_sf = 0; 

			iArray2DT Surf_Map (n_ews,n_esf+1);
			surf_set >> Surf_Map; // has dimension (n_ews x n_esf+1) must convert to (n_el x n_esf) next

			//-- Transfer surface data & keep tally of how many surfaces total this e_group
			
			for (e=0; e<n_ews; e++)
				for (s=0; s<n_esf; s++) 
					if ( Surf_Map(e,s+1) == 1 ) {
						Element_Set[e_set].n_sf++;
						Element_Set[e_set].Boundary_Surfaces ( Surf_Map(e,0)-1,s ) = Surf_Map ( e,s+1 ); 
					}

		}


}

//-----------------------------------------------------------------

void ContinuumT::Print(int input_file_type)
{

	if (input_file_type == kMonterey) {

		cout << "\n\n\n\n Contents of ContinuumT follows:\n\n";

		cout << "Xo: \n\n" << Xo 	<< "\n\n"; 
		//cout << "X:  \n\n" << X 	<< "\n\n"; 
		//cout << "x:  \n\n" << x 	<< "\n\n"; 
	
		for (int E=0; E<n_es; E++) {
			cout << "name: " 								<< Element_Set[E].name 									<< "\n";
			cout << "element_type: " 				<< Element_Set[E].element_type 					<< "\n";
			cout << "n_ip: " 								<< Element_Set[E].n_ip 									<< "\n";

			for (int f=0; f<n_fld; f++) {
				cout << "field: " 	<<f<< "\n";
				cout << "  constitutive_model: " 	<< Element_Set[E].constitutive_model[f]		<< "\n";
				cout << "  material set: " 				<< Element_Set[E].material_set[f] 				<< "\n";
				cout << "  material set item: " 	<< Element_Set[E].material_set_item[f] 		<< "\n";
			}
			cout << "geometric_parameter: " << Element_Set[E].geometric_parameter 	<< "\n";
			cout << "n_el: " 								<< Element_Set[E].n_el									<< "\n";
			cout << "n_en: " 								<< Element_Set[E].n_en									<< "\n";
			cout << "IEN Matrix: \n\n" 			<< Element_Set[E].IEN 									<< "\n\n";
			cout << "Bound Surf Matx: \n\n" << Element_Set[E].Boundary_Surfaces 		<< "\n\n";
		}

	} //-- End Monterey

	if (input_file_type == kTahoe) {

		cout << "n_np: " << n_np 		<< "\n\n"; 
		cout << "Xo: \n\n" << Xo 		<< "\n\n"; 
		cout << "n_el: " 						<< Element_Set[0].n_el	<< "\n";
		cout << "n_en: " 						<< Element_Set[0].n_en	<< "\n";
		cout << "n_sf: " 						<< Element_Set[0].n_sf	<< "\n";
		cout << "n_esf: " 					<< Element_Set[0].n_esf	<< "\n";
		cout << "IEN Matrix: \n\n" 	<< Element_Set[0].IEN 	<< "\n\n";
		cout << "\n\n" << Element_Set[0].Boundary_Surfaces 	<< "\n\n";

	} //-- End Tahoe

}
