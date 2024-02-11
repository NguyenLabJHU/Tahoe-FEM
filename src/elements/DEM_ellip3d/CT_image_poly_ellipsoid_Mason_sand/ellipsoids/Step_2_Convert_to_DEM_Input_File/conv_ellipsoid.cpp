/*
 * conv.cpp
 *
 *  Last updated: Oct 31, 2012
 *  By Boning Zhang
 *
 * Original Author: Tom Buzbee
 */

/*
 * conv.cpp
 *
 *  Previous updated: Feb 22, 2012
 *  By Yevgeniy Kaufman
 *
 * Previously Updated: Sep 17, 2008
 * Original Author: Tom Buzbee
 */

/*
//usage:
//
//g++ -c conv.cpp 
//g++ -o convert conv.o 
//
//This will generate the executable "convert" 
//
//Copy and paste the contents of grain.xls into a textfile, such as grain.txt, then run convert as 
//
//convert grain.txt out.txt 
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
using namespace std;

//#define LENGTH_CONVERSION 92165.8986   //1000000 microns to meters; 92165.8986 to convert pixels (10.85 microns/pixel) to meters 
#define LENGTH_CONVERSION 3.5*0.000001 // convert from mm to m, the original data in F-75 excel file is mm
const char TEMP_FILE[] = "temp";

// Reads a single line from an input stream (stops at newline or EOF)
string readLine(istream& stream)
{
	string line;
	char ch;

	do
	{
		stream.get(ch);
		line += ch;
	}
	while(ch != '\n' && ch != '\r' && !stream.eof());

	return line;
}

int main(int argc, char *argv[])
{
	if (argc != 3)
	{
		printf("usage: conv input output\n");
		return 0;
	}

	int num = 0;
	string line, output;

	// Fields from the input file
	float PI = 3.14159;
	int grainNumber, grainCoordination;
	double loc[3], volume, inscribeRadius, surfaceArea, aspectRatio;
	double xLength, yLength, zLength;
	double l1, l2, l3, m1, m2, m3, n1, n2, n3;
	double magnitude_l, magnitude_m, magnitude_n;
	double Radius_a, Radius_b, Radius_c;

	// Open the input and output files
	ifstream ifs(argv[1]);
	ofstream ofs(TEMP_FILE);
	ifstream ifs2(TEMP_FILE);
	ofstream ofs2(argv[2]);

	if(!ifs || !ofs || !ifs || !ofs2)
	{
		cout << "File I/O Error\n";
		return 1;
	}

	/* First pass */
	while(!ifs.eof())
	{
		line = readLine(ifs);
		istringstream iss(line);

		// Get each field
		iss >> grainNumber;
	        iss >> loc[0];
        	iss >> loc[1];
        	iss >> loc[2];
        	iss >> xLength;
		iss >> yLength;
		iss >> zLength;

// previously used input 
/*        iss >> axis[0];
        iss >> axis[1];
        iss >> axis[2];
*/	
        	iss >> l1;
        	iss >> l2;
        	iss >> l3;

		iss >> m1;
        	iss >> m2;
        	iss >> m3;

		iss >> n1;
        	iss >> n2;
        	iss >> n3;

//	iss >> volume; // Not used
//	iss >> surfaceArea; // Not used

		// See if there were any problems with input
		if (iss.fail())
		{
			cout << "Ignoring line: \n" << line << "\n\n";
			iss.clear();
		}
	  /*
		// Check for bad axes - came across a few of these in grain.xls
		else if (axis[0] == 0.0 && axis[1] == 0.0 && axis[2] == 0.0)
		{
			cout << "Ignoring line: \n" << line << "\n\n";
			iss.clear();
		}
	*/
		else
		{
// Previously used conversion by Yevgeniy
	//***********************IMPROPER ORIENTATIONS:
            // Convert the "axis" to the "axle" basis vectors
/*			axle[0][0] = (axis[0]*PI/180);
			axle[0][1] = cos(axis[1]*PI/180);	    // Use the major axis for the first axle
			axle[0][2] = cos(axis[2]*PI/180);

			axle[1][0] = cos(axis[0]*PI/180);
			axle[1][1] = cos(axis[1]*PI/180);	// Use a simple perpendicular for the second axle
			axle[1][2] = cos(axis[2]*PI/180);

			axle[2][0] = cos(axis[0]*PI/180);
			axle[2][1] = cos(axis[1]*PI/180);	// third axle = first X second
			axle[2][2] = cos(axis[2]*PI/180);
*/
// Below is modified
/*
Following is the thought for angle conversion
n1 n2 n3 are known.
Then I construct a vector that is perpendicular to v1 and name it
as v2.
l2 = -n1
m2 = 0
n2 = l1. (Since n2 is known, then from this we can get l1 = n2)

Then we can get v3 based on the perpendicular relations.
l3 = m1*n2-n1*m2
m3 = n1*l2-l1*n2
n3 = l1*m2-m1*l2. (Since n3 is known, then from this we can get m1
= (l1*m2-3)/l2, m2 and l2 are determined before)

*/

			magnitude_l = sqrt(pow(l1, 2) + pow(l2, 2) + pow(l3, 2));
            		magnitude_m = sqrt(pow(m1, 2) + pow(m2, 2) + pow(m3, 2));
           		magnitude_n = sqrt(pow(n1, 2) + pow(n2, 2) + pow(n3, 2));

			l1 = acos(l1/magnitude_l);
			l2 = acos(l2/magnitude_l);
			l3 = acos(l3/magnitude_l);

			m1 = acos(m1/magnitude_m);
			m2 = acos(m2/magnitude_m);
			m3 = acos(m3/magnitude_m);

			n1 = acos(n1/magnitude_n);
			n2 = acos(n2/magnitude_n);
			n3 = acos(n3/magnitude_n);

			// Start writing the current row
			ofs << setw(9) << grainNumber;		//ID
			ofs << setw(10) << 0;				//type
			ofs.setf(ios::scientific);
			ofs.precision(6);

			// Write the radii in order of decreasing magnitude
			//if (aspectRatio >= 1)
			Radius_a = 0.5*xLength;
			Radius_b = 0.5*yLength;
			Radius_c = 0.5*zLength;
			//aspectRatio = majorLength / minorLength;

			ofs << setw(16) << Radius_a*LENGTH_CONVERSION;
			ofs << setw(16) << Radius_b*LENGTH_CONVERSION;
			ofs << setw(16) << Radius_c*LENGTH_CONVERSION;
			//if (aspectRatio < 1)
				//ofs << setw(16) << inscribeRadius * aspectRatio / LENGTH_CONVERSION;

			// Write the position
			ofs << setw(16) << loc[0]*LENGTH_CONVERSION;			//position_x
			ofs << setw(15) << loc[1]*LENGTH_CONVERSION;			//position_y
			ofs << setw(16) << loc[2]*LENGTH_CONVERSION;			//position_z
			ofs.precision(12);

			// Write the axle vectors in the same order as the radii
			ofs << setw(22) << l1;
			ofs << setw(22) << l2;
			ofs << setw(22) << l3;
			ofs << setw(22) << m1;
			ofs << setw(22) << m2;
			ofs << setw(22) << m3;
			ofs << setw(22) << n1;
			ofs << setw(22) << n2;
			ofs << setw(22) << n3;

			// Write rest of the values as zero
			ofs.precision(6);
			ofs << setw(16) << 0.0;				//velocity_x
			ofs << setw(16) << 0.0;				//velocity_y
			ofs << setw(16) << 0.0;				//velocity_z
			ofs << setw(16) << 0.0;				//omga_x
			ofs << setw(16) << 0.0;				//omga_y
			ofs << setw(16) << 0.0;				//omga_z
			ofs << setw(16) << 0.0;				//force_x
			ofs << setw(16) << 0.0;				//force_y
			ofs << setw(16) << 0.0;				//force_z
			ofs << setw(16) << 0.0;				//moment_x
			ofs << setw(16) << 0.0;				//moment_y
			ofs << setw(16) << 0.0 << "\n";		//moment_z

			num++;
		}
	}

	ifs.close();
	ofs.close();
	cout << num << " particles found.\n";

	/* Second pass - adds information to the front of the file */

	// Write simulation info
	ofs2.setf(ios::scientific);
	ofs2.precision(6);
	ofs2 << setw(10) << num << "         1" << "\n";
	ofs2 << setw(16) << 0.0;
	ofs2 << setw(16) << 0.0;		// Centroid of simulation volume
	ofs2 << setw(16) << 0.0;
	ofs2 << setw(16) << 0.5;
	ofs2 << setw(16) << 0.5;		// Simulation dimensions
	ofs2 << setw(16) << 0.5 << "\n";

	// Write column names
	ofs2 << "     ID         type     radius_a        radius_b" <<
			"        radius_c       position_x      position_y" <<
			"      position_z       axle_a_x              axle_a_y" <<
			"              axle_a_z              axle_b_x     " <<
			"         axle_b_y              axle_b_z          " <<
			"    axle_c_x              axle_c_y              " <<
			"axle_c_z             velocity_x      velocity_y  " <<
			"    velocity_z        omga_x          omga_y     " <<
			"     omga_z         force_x         force_y      " <<
			"   force_z         moment_x        moment_y        moment_z\n";

	//Copy the contents of the first file into the second
	while(!ifs2.eof())
	{
		ofs2 << (char) ifs2.get();
	}

	ifs2.close();
	ofs2.close();
	remove(TEMP_FILE);

	return 0;
}
