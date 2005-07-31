/* $Id: main.cpp,v 1.1.1.1 2002-11-15 00:51:07 paklein Exp $ */
#include <iostream.h>
#include "fstreamT.h"

/* files needed for MLS */
#include "ParentDomainT.h"
#include "D2MeshFreeSupport2DT.h"

/* I/O */
#include "OutputBaseT.h"
#include "OutputSetT.h"

/* Cauchy-Born constitutive model */
//#include "FDMatSupportT.h"
//#include "VIB2D.h"
//#include "LJTr2D.h"

//TEMP
//void WaveSpeeds(const dArrayT& normal, dArrayT& speeds)
#include "dSymMatrixT.h"
#include "MaterialSupportT.h"

enum Widths { iwidth = 10, dwidth = 12, dprecision = 5};

using namespace Tahoe;

int main(void)
{
	MaterialSupportT support(2,2,1);


	/* construct material support */
	int nsd  = 2;
	int ndof = 2;
	int nip  = 1;
//	FDMatSupportT material_support(nsd, ndof, nip);
	
	ArrayT<dMatrixT> F_List(nip);
	F_List[0].Dimension(nsd);
//	material_support.SetDeformationGradient(&F_List);

	/* construct constitutive model */
    StringT materials_params_file;
	cout << " material model parameters: ";
	cin >> materials_params_file;
	ifstreamT material_in('#', materials_params_file);
	if (!material_in.is_open())
	  {
		cout << " error opening file: " << materials_params_file << endl;
		throw;
	  }

//	VIB2D material(material_in, material_support);
//	material.Initialize();

    /* window function parameters */
    StringT window_params_file;
	cout << " file of window function parameters: ";
	cin >> window_params_file;
	ifstreamT window_in('#', window_params_file);
	if (!window_in.is_open())
	  {
		cout << " error opening file: " << window_params_file << endl;
		throw;
	  }

	StringT geo_file;
	cout << " geo_file (file or \"quit\"): ";
	cin >> geo_file;
	
	while (geo_file != "quit")
	{
		/* skip comments */
		ifstreamT geo(geo_file);
		if (!geo.is_open())
		{
			cout << " file not found" << endl;	
			throw;
		}
		
		for (int i = 0; i < 5; i++)
			geo.clear_line();
			
		/* read coordinates */	
		int nnd;
		geo >> nnd;	
		cout << " number of nodes: " << nnd << endl;
	
		int nsd = 2;
		//cout << " number of spatial dimensions: ";
		//cin >> nsd;	
		
		dArray2DT coords(nnd, nsd);
		if (nsd == 2) /* translate 3D coordinate file into 2D by dropping z */
		{
			dArray2DT coords_3D(nnd, 3);
			coords_3D.ReadNumbered(geo);
			
			dArrayT column(nnd);
			coords_3D.ColumnCopy(0, column);
			coords.SetColumn(0, column);
			coords_3D.ColumnCopy(1, column);
			coords.SetColumn(1, column);
		}
		else
			coords.ReadNumbered(geo);

		/* set MLS */
		iArray2DT junk2D(0,4);
		iArrayT   junk1D;
		ParentDomainT dummy_domain(GeometryT::kQuadrilateral, 1, 4);
		D2MeshFreeSupport2DT D2MLS(dummy_domain, coords, junk2D, junk1D, window_in);

		/* cutting facet */
		dArray2DT facet_coords(1,4);
		facet_coords(0,0) = -110.0;
		facet_coords(0,1) = -8.835;
		facet_coords(0,2) = -72.0;
		facet_coords(0,3) = -8.835;
		D2MLS.SetCuttingFacets(facet_coords, 2);

		/* set support size */
		double d_max;
		cout << " support size: ";
		cin >> d_max;
		if (d_max < 0)
			throw;
		else
		{
			iArrayT range(nnd);
			range.SetValueToPosition();
			dArray2DT Dmax(nnd, 1);
			Dmax = d_max;
			
			D2MLS.SetSupportParameters(range, Dmax);
		}
		
		/* write info */
		D2MLS.WriteStatistics(cout);
		cout << '\n';

		StringT disp_file;
		cout << " disp_file (file or \"quit\"): ";
		cin >> disp_file;
		while (disp_file != "quit")
		{
			/* skip comments */
			ifstreamT disp(disp_file);
			if (!disp.is_open())
			{
				cout << " file not found" << endl;	
				throw;
			}
			disp.clear_line();
		
			/* read displacements */
			dArray2DT displacements(nnd, nsd);
			if (nsd == 2)
			{
				dArray2DT d_3D(nnd, 3);
				disp >> d_3D;
				
				dArrayT column(nnd);
				d_3D.ColumnCopy(0, column);
				displacements.SetColumn(0, column);
				d_3D.ColumnCopy(1, column);
				displacements.SetColumn(1, column);
			}
			else
				disp >> displacements;
				
			/* streams */
			StringT root;
			root.Root(disp_file);

			/* set I/O */
			IOBaseT::FileTypeT file_type = IOBaseT::kEnSight;
			OutputBaseT* output = IOBaseT::NewOutput("atom_field", "0.0", "testing", disp_file, file_type, cout);
			output->SetCoordinates(coords, NULL);

			/* define "connectivities" */
			iArrayT range(coords.MajorDim());
			range.SetValueToPosition();
			iArray2DT connectivities(range.Length(), 1, range.Pointer());

			/* define output set */
			ArrayT<StringT> n_labels(20);
			n_labels[0] = "D_X";
			n_labels[1] = "D_Y";

			n_labels[2] = "DXDX";
			n_labels[3] = "DXDY";
			n_labels[4] = "DYDX";
			n_labels[5] = "DYDY";
			n_labels[6] = "strain_norm";
			n_labels[7] = "grad_norm";

			n_labels[8]  = "DXDXX";
			n_labels[9]  = "DXDYY";
			n_labels[10] = "DXDXY";
			n_labels[11] = "DYDXX";
			n_labels[12] = "DYDYY";
			n_labels[13] = "DYDXY";
			n_labels[14] = "grad_grad_norm";

			n_labels[15] = "S11";
			n_labels[16] = "S22";
			n_labels[17] = "S12";
			n_labels[18] = "CD";
			n_labels[19] = "CS";
			
			OutputSetT output_set(GeometryT::kPoint, connectivities, n_labels);
			int io_ID = output->AddElementSet(output_set);

			/* array of output values */
			dArray2DT e_values, n_values(connectivities.MajorDim(), n_labels.Length());
			dArrayT n_value;

			/* loop over nodes */
			int update = nnd/10;
			int count = 0, wrap = 0;
			dArrayT x(nsd);
			dArray2DT loc_disp(0, nsd);
			nVariArray2DT<double> loc_disp_man(0, loc_disp, nsd);
			dArrayT vec;
			VariArrayT<double> vec_man(0, vec);
			dArrayT normal(nsd), speeds(nsd);
			for (int i = 0; i < nnd; i++)
			{
				/* alias row of output table */
				n_values.RowAlias(i, n_value);
				n_value = 0.0;
			
				/* progress */
				if (++count == update)
				{
					cout << " node: " << i+1 << endl;
					count = 0;
				}
			
				/* set field */
				coords.RowAlias(i, x);
				if (D2MLS.SetFieldAt(x))
				{
					/* list of nodal neighbors */
					const ArrayT<int>& nodes = D2MLS.NeighborsAt();
			
					/* collect local displacements */
					loc_disp_man.SetMajorDimension(nodes.Length(), false);
					loc_disp.RowCollect(nodes, displacements);
					vec_man.SetLength(nodes.Length(), false);
					
					/* displacements */
					n_value[0] = displacements(i,0);
					n_value[1] = displacements(i,1);
				
					/* compute 1st gradients */
					const dArray2DT& Dphi = D2MLS.DFieldAt();

					loc_disp.ColumnCopy(0, vec);
					double& u1_1 = n_value[2] = Dphi.DotRow(0, vec.Pointer());
					double& u1_2 = n_value[3] = Dphi.DotRow(1, vec.Pointer());
				
					loc_disp.ColumnCopy(1, vec);
					double& u2_1 = n_value[4] = Dphi.DotRow(0, vec.Pointer());
					double& u2_2 = n_value[5] = Dphi.DotRow(1, vec.Pointer());

					/* strain norm */
					n_value[6] = sqrt(u1_1*u1_1 + u2_2*u2_2 + 0.5*(u1_2*u1_2 + u2_1*u2_1 + 2.0*(u1_2*u2_1)));
					n_value[7] = sqrt(u1_1*u1_1 + u1_2*u1_2 + u2_1*u2_1 + u2_2*u2_2);

					/* compute 2nd gradients */
					const dArray2DT& DDphi = D2MLS.DDFieldAt();

					loc_disp.ColumnCopy(0, vec);
					double& u1_11 = n_value[8]  = DDphi.DotRow(0, vec.Pointer());
					double& u1_22 = n_value[9]  = DDphi.DotRow(1, vec.Pointer());
					double& u1_12 = n_value[10] = DDphi.DotRow(2, vec.Pointer());

					loc_disp.ColumnCopy(1, vec);
					double& u2_11 = n_value[11] = DDphi.DotRow(0, vec.Pointer());
					double& u2_22 = n_value[12] = DDphi.DotRow(1, vec.Pointer());
					double& u2_12 = n_value[13] = DDphi.DotRow(2, vec.Pointer());
	
					/* strain gradient norm */
					n_value[14] = sqrt(u1_11*u1_11 + 2.0*u1_12*u1_12 + u1_22*u1_22 +
									   u2_11*u2_11 + 2.0*u2_12*u2_12 + u2_22*u2_22);

					/* evaluate continuum materials response */
					dMatrixT& F = F_List[0];
					F.Identity();
					F(0,0) += u1_1;
					F(0,1) += u1_2;
					F(1,0) += u2_1;
					F(1,1) += u2_2;
					
					/* cauchy stress */
					//const dSymMatrixT& cauchy = material.s_ij();
					dSymMatrixT cauchy(2);
					n_value[15] = cauchy(0,0);
					n_value[16] = cauchy(1,1);
					n_value[17] = cauchy(0,1);
					
					/* wave speeds */
					normal[0] = 1.0;
					normal[1] = 0.0;
					//material.WaveSpeeds(normal, speeds);
					n_value[18] = speeds[0];
					n_value[19] = speeds[1];
	
#if 0
					if (grad_grad_norm == 0.0) {
						l1_out << setw(dwidth) << 1000.0;
						l2_out << setw(dwidth) << 1000.0;
					} else {
						l1_out << setw(dwidth) << strain_norm/grad_grad_norm;
						l2_out << setw(dwidth) << grad_norm/grad_grad_norm;
					}
#endif
				}
				else /* could not fit MLS field */
				{
					/* only displacements */
					n_value[0] = displacements(i,0);
					n_value[1] = displacements(i,1);
				}
				
			}
			
			/* write output */
			output->WriteOutput(1.0, io_ID, n_values, e_values);
			
			/* clean up */
			delete output;

			/* next file */
			cout << " disp_file (file or \"quit\"): ";
			cin >> disp_file;
		}

		/* next */
		cout << " geo_file (file or \"quit\"): ";
		cin >> geo_file;
	}

	cout << "\n end." << endl;
}
