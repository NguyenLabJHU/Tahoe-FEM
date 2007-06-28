#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iomanip>   // in order to use setw
#include <math.h>

using namespace std;

int main()
{

    ifstream fInput ("input.txt");
    ofstream fOutput;
    fOutput.open ("output.txt");
    string line;
    for (int i=0; i<3; i++)
      getline(fInput, line); 
    fOutput << line <<"\n";
    int n_sd;
    fInput >> n_sd;
    fOutput << "number of space dimension = " << n_sd <<"\n";
    getline(fInput, line); 
    int n_nodes,n_elements,n_element_blocks;
    fInput >> n_nodes >> n_elements >> n_element_blocks;
    fOutput << "number of nodes = " << n_nodes <<"\n";
    fOutput << "number of elements = " << n_elements <<"\n";
    fOutput << "number of element blocks = " << n_element_blocks <<"\n";
    getline(fInput, line); 
    int n_node_sets, n_side_sets;
    fInput >> n_node_sets >> n_side_sets;
    fOutput << "number of node sets = " << n_node_sets <<"\n";
    fOutput << "number of side sets = " << n_side_sets <<"\n \n";
    for (int i=0; i<6;i++)
	getline(fInput, line);
    fOutput << line <<"\n";
    double x[n_nodes],y[n_nodes],z[n_nodes];
    for (int i=0; i<n_nodes; i++)
    {
	fInput >> x[i] >> y[i] >> z[i];
	getline(fInput, line);
    }
    fOutput << setiosflags(ios::scientific) << setprecision(7);
    for (int i=0; i<n_nodes; i++)
	fOutput << setw(5) << i+1 << setw(20) << x[i] << setw(20) << y[i] << setw(20) << z[i] << "\n";

    fOutput << "\n";
    getline(fInput, line); 
    while (line !="! Element block    1")
      getline(fInput, line); 

	fOutput << line <<"\n";

  int*** Connectivity;
  Connectivity=new int**[n_element_blocks];
  int  n_elements_in_block[n_element_blocks];

    fOutput << "number of element blocks = " << n_element_blocks <<"\n";
    for (int z=0; z<n_element_blocks; z++)
    {
	fOutput << line <<"\n";
	int temp,elements_in_block;
	fInput >> temp >> elements_in_block;
	for (int i=0; i<3; i++)
	    getline(fInput, line);
	fOutput << line <<"\n";

	n_elements_in_block[z]=elements_in_block;

	Connectivity[z]=new int* [elements_in_block];

	for(int k=0;k<elements_in_block;k++)
	    Connectivity[z][k]=new int[27]; 
	for (int k=0; k<elements_in_block; k++)
	{
	    for(int i=0; i<8; i++) 
		fInput >> Connectivity[z][k][i];    
	    getline(fInput, line);
	    for(int i=8; i<16; i++) 
		fInput >> Connectivity[z][k][i];    
	    getline(fInput, line);
	    for(int i=16; i<24; i++) 
		fInput >> Connectivity[z][k][i];    
	    getline(fInput, line);
	    for(int i=24; i<27; i++) 
		fInput >> Connectivity[z][k][i];    
	    getline(fInput, line);
	}
	
	for (int k=0; k<elements_in_block; k++)
	{
	    fOutput << setw(5) << k+1;
	    for(int i=0; i<12; i++) 
		fOutput << setw(5) << Connectivity[z][k][i];  

	    fOutput << setw(5) << Connectivity[z][k][16] <<setw(5) << Connectivity[z][k][17] <<setw(5) << Connectivity[z][k][18] <<setw(5) << Connectivity[z][k][19] ; 
	    fOutput << setw(5) << Connectivity[z][k][12] <<setw(5) << Connectivity[z][k][13] <<setw(5) << Connectivity[z][k][14] <<setw(5) << Connectivity[z][k][15] ;   
	    fOutput << setw(5) << Connectivity[z][k][21] <<setw(5) << Connectivity[z][k][22] <<setw(5) << Connectivity[z][k][25] <<setw(5) << Connectivity[z][k][26] ; 
	    fOutput << setw(5) << Connectivity[z][k][23] <<setw(5) << Connectivity[z][k][24] <<setw(5) << Connectivity[z][k][20] ;   
 
	    fOutput << "\n";
	}
	getline(fInput, line); 
    }


    fOutput << "\n";
    fOutput << "number of node sets = " << n_node_sets <<"\n";
    for (int i=0; i<2; i++)
	    getline(fInput, line);
    int temp,nodes_in_set;
    int** Node_Sets;
    Node_Sets=new int*[n_node_sets];
    int n_nodes_in_set[n_node_sets];

    for (int z=0; z<n_node_sets; z++)
    {
	getline(fInput, line);
	fOutput << line <<"\n";
	fInput >> temp >> nodes_in_set;
	getline(fInput, line);
	Node_Sets[z]=new int [nodes_in_set];
	n_nodes_in_set[z] = nodes_in_set;
	for (int i=0; i<nodes_in_set; i++)
	{
	    fInput >> Node_Sets[z][i] >>temp ;
	    getline(fInput, line);
	}

	fOutput << "number of nodes in set = " << nodes_in_set <<"\n";
	for (int i=0; i<nodes_in_set; i++)
	    fOutput << setw(5) << Node_Sets[z][i];
	fOutput << "\n";

    }

    for (int i=0; i<2; i++)
    getline(fInput, line); 

/* creating .geom file */
    ofstream fGeom;
    fGeom.open ("output.geom");
    fGeom << "*version" <<"\n";
    fGeom << "1.0" <<"\n";
    fGeom << "*title" <<"\n";
    fGeom << ".geom file converted from cubit text output" <<"\n";
    fGeom << "*dimensions" <<"\n";
    fGeom << n_nodes <<setw(20)<<"# number of nodes" << "\n";
    fGeom << n_sd <<setw(35)<<"# number of spatial dimensions" << "\n \n";
    fGeom << n_element_blocks <<setw(29)<<"# number of element sets" << "\n";
    fGeom << "# [ID] [nel] [nen]" << "\n";
    for (int z=0; z<n_element_blocks; z++)
	fGeom << z+1 << setw(5) << n_elements_in_block[z] << setw(5) << 27 <<"\n";
    fGeom << "\n";
    fGeom << n_node_sets <<setw(26)<<"# number of node sets" << "\n";
    fGeom << "# [ID] [nnd]" << "\n";
    for (int z=0; z<n_node_sets; z++)
	fGeom << z+1 << setw(5) << n_nodes_in_set[z] <<"\n";
    fGeom << "\n";
    fGeom << 0 <<setw(26)<<"# number of side sets" << "\n"; 
    fGeom << "# end dimensions" << "\n"; 
    fGeom << "\n";
    fGeom << "*nodesets" << "\n"; 
    for (int z=0; z<n_node_sets; z++)
    {
    fGeom << "*set" << "\n"; 
    fGeom << n_nodes_in_set[z] << setw(19)<<"# number of nodes" << "\n";
    for (int i=0; i<n_nodes_in_set[z]; i++)
	fGeom << setw(5) << Node_Sets[z][i];
	fGeom << "\n";
    }
    fGeom << "# end node sets" << "\n";
    fGeom << "\n";

    fGeom << "*sidesets" << "\n";
    fGeom << "\n";

    fGeom << "*elements" << "\n";
    for (int z=0; z<n_element_blocks; z++)
    {
    fGeom << "*set" << "\n"; 
    fGeom << n_elements_in_block[z] << setw(22) << "# number of elements" <<"\n";
    fGeom << 27 << setw(27) << "# number of element nodes" <<"\n";
    for (int k=0; k<n_elements_in_block[z]; k++)
    {
	fGeom <<  k+1;
	for(int i=0; i<12; i++) 
	    fGeom << setw(7) << Connectivity[z][k][i];  

	fGeom << setw(7) << Connectivity[z][k][16] <<setw(7) << Connectivity[z][k][17] <<setw(7) << Connectivity[z][k][18] <<setw(7) << Connectivity[z][k][19] ; 
	fGeom << setw(7) << Connectivity[z][k][12] <<setw(7) << Connectivity[z][k][13] <<setw(7) << Connectivity[z][k][14] <<setw(7) << Connectivity[z][k][15] ;   
	fGeom << setw(7) << Connectivity[z][k][21] <<setw(7) << Connectivity[z][k][22] <<setw(7) << Connectivity[z][k][25] <<setw(7) << Connectivity[z][k][26] ; 
	fGeom << setw(7) << Connectivity[z][k][23] <<setw(7) << Connectivity[z][k][24] <<setw(7) << Connectivity[z][k][20] ;   
 
	fGeom << "\n";
    }

    }
    fGeom << "# end elements" << "\n";
    fGeom << "\n";
    fGeom << "*nodes" << "\n";
    fGeom << n_nodes << setw(19) << "# number of nodes" <<"\n";
    fGeom << n_sd << setw(34) << "# number of spatial dimensions" <<"\n";
    fGeom << setiosflags(ios::scientific) << setprecision(7);
    for (int i=0; i<n_nodes; i++)
	fGeom << setw(5) << i+1 << setw(20) << x[i] << setw(20) << y[i] << setw(20) << z[i] << "\n";


    fOutput.close();
    fInput.close();
    fGeom.close();
    delete [] Connectivity;
    delete [] Node_Sets;

    return 0;
}
