#include <iostream.h>

#include "fstreamT.h"
#include "StringT.h"
#include "dArray2DT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "ParentDomainT.h"
#include "D2MeshFreeSupport2DT.h"
#include "nVariArray2DT.h"
#include "VariArrayT.h"

static void OpenStreams(const StringT& root, int nsd, ArrayT<ofstream>& Du_out, 
	ArrayT<ofstream>& DDu_out);

enum Widths { iwidth = 10, dwidth = 12, dprecision = 5};

int main(void)
{
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
		if (nsd == 2)
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
			ArrayT<ofstream> Du_out;  // first gradients	
			ArrayT<ofstream> DDu_out; // second gradients	
			OpenStreams(root, nsd, Du_out, DDu_out);

			/* extra streams */
			StringT dd_file = root;
			dd_file.Append(".dd");
			ofstream dd_out(dd_file);
			dd_out.setf(ios::showpoint);
			dd_out.setf(ios::right, ios::adjustfield); 
			dd_out.setf(ios::scientific, ios::floatfield);
			dd_out.precision(dprecision);
			dd_out << "<comment line>\n";

			StringT l1_file = root;
			l1_file.Append(".l1");
			ofstream l1_out(l1_file);
			l1_out.setf(ios::showpoint);
			l1_out.setf(ios::right, ios::adjustfield); 
			l1_out.setf(ios::scientific, ios::floatfield);
			l1_out.precision(dprecision);
			l1_out << "<comment line>\n";

			StringT l2_file = root;
			l2_file.Append(".l2");
			ofstream l2_out(l2_file);
			l2_out.setf(ios::showpoint);
			l2_out.setf(ios::right, ios::adjustfield); 
			l2_out.setf(ios::scientific, ios::floatfield);
			l2_out.precision(dprecision);
			l2_out << "<comment lines>\n";

			StringT x0_file = root;
			x0_file.Append(".x0");
			ofstream x0_out(x0_file);
			x0_out.setf(ios::showpoint);
			x0_out.setf(ios::right, ios::adjustfield); 
			x0_out.setf(ios::scientific, ios::floatfield);
			x0_out.precision(dprecision);

			/* loop over nodes */
			int update = nnd/10;
			int count = 0, wrap = 0;
			dArrayT x(nsd);
			dArray2DT loc_disp(0, nsd);
			nVariArray2DT<double> loc_disp_man(0, loc_disp, nsd);
			dArrayT vec;
			VariArrayT<double> vec_man(0, vec);
			for (int i = 0; i < nnd; i++)
			{
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
				
				/* compute 1st gradients */
				const dArray2DT& Dphi = D2MLS.DFieldAt();

				loc_disp.ColumnCopy(0, vec);
				double u1_1 = Dphi.DotRow(0, vec.Pointer());
				double u1_2 = Dphi.DotRow(1, vec.Pointer());
				
				loc_disp.ColumnCopy(1, vec);
				double u2_1 = Dphi.DotRow(0, vec.Pointer());
				double u2_2 = Dphi.DotRow(1, vec.Pointer());

				/* output */
				Du_out[0] << setw(dwidth) << u1_1;
				Du_out[1] << setw(dwidth) << u1_2;
				Du_out[2] << setw(dwidth) << u2_1;
				Du_out[3] << setw(dwidth) << u2_2;

				/* strain norm */
				double strain_norm = sqrt(u1_1*u1_1 + u2_2*u2_2 + 0.5*(u1_2*u1_2 + u2_1*u2_1 + 2.0*(u1_2*u2_1)));
				double grad_norm = sqrt(u1_1*u1_1 + u1_2*u1_2 + u2_1*u2_1 + u2_2*u2_2);

				/* compute 2nd gradients */
				const dArray2DT& DDphi = D2MLS.DDFieldAt();

				loc_disp.ColumnCopy(0, vec);
				double u1_11 = DDphi.DotRow(0, vec.Pointer());
				double u1_22 = DDphi.DotRow(1, vec.Pointer());
				double u1_12 = DDphi.DotRow(2, vec.Pointer());

				loc_disp.ColumnCopy(1, vec);
				double u2_11 = DDphi.DotRow(0, vec.Pointer());
				double u2_22 = DDphi.DotRow(1, vec.Pointer());
				double u2_12 = DDphi.DotRow(2, vec.Pointer());

				/* strain gradient norm */
				double grad_grad_norm = sqrt(u1_11*u1_11 + 2.0*u1_12*u1_12 + u1_22*u1_22 +
											 u2_11*u2_11 + 2.0*u2_12*u2_12 + u2_22*u2_22);

				/* output */
				DDu_out[0] << setw(dwidth) << u1_11;
				DDu_out[1] << setw(dwidth) << u1_22;
				DDu_out[2] << setw(dwidth) << u1_12;

				DDu_out[3] << setw(dwidth) << u2_11;
				DDu_out[4] << setw(dwidth) << u2_22;
				DDu_out[5] << setw(dwidth) << u2_12;

				dd_out << setw(dwidth) << grad_grad_norm;

				if (grad_grad_norm == 0.0) {
				  l1_out << setw(dwidth) << 1000.0;
				  l2_out << setw(dwidth) << 1000.0;
				}
				else {
				  l1_out << setw(dwidth) << strain_norm/grad_grad_norm;
				  l2_out << setw(dwidth) << grad_norm/grad_grad_norm;				  
				}

				/* data along y =-10.58 (the fracture plane */
				if (fabs(x[1] + 10.58) < 0.5)
				  {
					x0_out << setw(iwidth) << i+1
                         << setw(dwidth) << x[0]
						 << setw(dwidth) << x[1]
						 << setw(dwidth) << displacements(i,0)
						 << setw(dwidth) << displacements(i,1)
						 << setw(dwidth) << u1_1
						 << setw(dwidth) << u1_2
						 << setw(dwidth) << u2_1
						 << setw(dwidth) << u2_2
						 << setw(dwidth) << strain_norm
						 << setw(dwidth) << grad_norm
						 << setw(dwidth) << u1_11
						 << setw(dwidth) << u1_22
						 << setw(dwidth) << u1_12
						 << setw(dwidth) << u2_11
						 << setw(dwidth) << u2_22
						 << setw(dwidth) << u2_12
						 << setw(dwidth) << grad_grad_norm;

					if (grad_grad_norm == 0.0) {
					  x0_out << setw(dwidth) << 1000.0;
					  x0_out << setw(dwidth) << 1000.0;
					}
					else {
					  x0_out << setw(dwidth) << strain_norm/grad_grad_norm;
					  x0_out << setw(dwidth) << grad_norm/grad_grad_norm;				  
					}
					x0_out << '\n';
				  }
				  }
				else
				  {
				Du_out[0] << setw(dwidth) << 0.0;
				Du_out[1] << setw(dwidth) << 0.0;
				Du_out[2] << setw(dwidth) << 0.0;
				Du_out[3] << setw(dwidth) << 0.0;

				DDu_out[0] << setw(dwidth) << 0.0;
				DDu_out[1] << setw(dwidth) << 0.0;
				DDu_out[2] << setw(dwidth) << 0.0;
				DDu_out[3] << setw(dwidth) << 0.0;
				DDu_out[4] << setw(dwidth) << 0.0;
				DDu_out[5] << setw(dwidth) << 0.0;

				dd_out << setw(dwidth) << 0.0;

				l1_out << setw(dwidth) << 0.0;
				l2_out << setw(dwidth) << 0.0;
				  }

				/* wrap values */
				if (++wrap == 6)
				{
					Du_out[0] << '\n';
					Du_out[1] << '\n';
					Du_out[2] << '\n';
					Du_out[3] << '\n';
				
					DDu_out[0] << '\n';
					DDu_out[1] << '\n';
					DDu_out[2] << '\n';
					DDu_out[3] << '\n';				
					DDu_out[4] << '\n';
					DDu_out[5] << '\n';				
	
					dd_out << '\n';				

					l1_out << '\n';				
					l2_out << '\n';
					wrap = 0;
				}
			}

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

void OpenStreams(const StringT& root, int nsd, ArrayT<ofstream>& Du_out, 
	ArrayT<ofstream>& DDu_out)
{
	if (nsd == 2)
	{
		StringT file;	
		Du_out.Allocate(4);
		DDu_out.Allocate(6);

		/* first gradient */
		file = root;
		file.Append(".ux_x");
		Du_out[0].open(file);
		Du_out[0] << " <comment line>\n";

		file = root;
		file.Append(".ux_y");
		Du_out[1].open(file);
		Du_out[1] << " <comment line>\n";

		file = root;
		file.Append(".uy_x");
		Du_out[2].open(file);
		Du_out[2] << " <comment line>\n";

		file = root;
		file.Append(".uy_y");
		Du_out[3].open(file);
		Du_out[3] << " <comment line>\n";

		/* second gradient */
		file = root;
		file.Append(".ux_xx");
		DDu_out[0].open(file);
		DDu_out[0] << " <comment line>\n";

		file = root;
		file.Append(".ux_yy");
		DDu_out[1].open(file);
		DDu_out[1] << " <comment line>\n";

		file = root;
		file.Append(".ux_xy");
		DDu_out[2].open(file);
		DDu_out[2] << " <comment line>\n";

		file = root;
		file.Append(".uy_xx");
		DDu_out[3].open(file);
		DDu_out[3] << " <comment line>\n";

		file = root;
		file.Append(".uy_yy");
		DDu_out[4].open(file);
		DDu_out[4] << " <comment line>\n";

		file = root;
		file.Append(".uy_xy");
		DDu_out[5].open(file);
		DDu_out[5] << " <comment line>\n";
	}
	else
	{
		cout << " no 3D yet" << endl;
		throw;
	}

	/* set prefs */
	for (int i = 0; i < Du_out.Length(); i++)
	{
		Du_out[i].setf(ios::showpoint);
		Du_out[i].setf(ios::right, ios::adjustfield); 
		Du_out[i].setf(ios::scientific, ios::floatfield);
		Du_out[i].precision(dprecision);
	}	

	for (int j = 0; j < DDu_out.Length(); j++)
	{
		DDu_out[j].setf(ios::showpoint);
		DDu_out[j].setf(ios::right, ios::adjustfield); 
		DDu_out[j].setf(ios::scientific, ios::floatfield);
		DDu_out[j].precision(dprecision);
	}	

	return 0;
}
