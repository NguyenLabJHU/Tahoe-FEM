#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

const long double PI = 3.1415927;
const int NODE = 134;
const int ELEM = 144;

int main(int argc, char *argv[])
{
    if(argc < 7)
    {
	cout << endl
	     << "USAGE: plot1p  particle_number  contact_file_prefix  particle_file_prefix  first_suffix  last_suffix  suffix_increment" << endl
	     << "   --example: plot1p 215 dep_contact dep_particle 0 100 5" << endl << endl;
	return -1;
    }	

    int first, last, incre, pnumber;
    first   = atoi(argv[4]);
    last    = atoi(argv[5]);
    incre   = atoi(argv[6]);
    pnumber = atoi(argv[1]);

    ifstream ifs1,ifs2;
    ofstream ofs;
    char filein1[50], filein2[50];
    char fileout[50];
    char num[4], s[20];

    int TotalCntct, ptcl_1, ptcl_2;
    long double point1_x, point1_y, point1_z, point2_x, point2_y, point2_z;
    long double radius_1, radius_2, penetration, tang_dispmt, contact_radius, R0, E0, normal_force, shear_force;

    int ID, type, TotalNum, RORS;
    long double cx, cy, cz, rd, wd, lt, ht;
    long double a, b, c, x0, y0, z0, l1, l2, l3, m1, m2, m3, n1, n2, n3, tmp;
    long double x, y, z, xp, yp, zp, t, theta; // x, y, z are local coodinates, xp, yp, zp are global.
    int n, i, j, ic;

    char pnum[6];
    sprintf(pnum, "%05d", pnumber);

    for(n=first; n<=last; n+=incre)
    {
	sprintf(num, "%03d", n);
	strcpy(filein1, argv[2]);
	strcat(filein1, "_");
	strcat(filein1, num);
	strcpy(filein2, argv[3]);
	strcat(filein2, "_");
	strcat(filein2, num);

	strcpy(fileout, "particle_");
	strcat(fileout, pnum);
	strcat(fileout, "_");
	strcat(fileout, num);
	strcat(fileout, ".dat");
        cout << "generating file " << fileout << " ......" << endl;

	ifs1.open(filein1);
	if(!ifs1)  { cout<<"ifs1 stream error!"<<endl; exit(-1);}
	ifs2.open(filein2);
	if(!ifs2)  { cout<<"ifs2 stream error!"<<endl; exit(-1);}
	ofs.open(fileout);
	if(!ofs)   { cout<<"ofs stream error!"<<endl; exit(-1);}
	ofs.setf(ios::scientific, ios::floatfield);

	ifs1 >> TotalCntct;
	ifs1 >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;

	ifs2 >> TotalNum >> RORS;
	if(RORS==0)
	    ifs2 >> cx >> cy >> cz >> rd >> ht;
	else
	    ifs2 >> cx >> cy >> cz >> wd >> lt >> ht;
	ifs2>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
	    >>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;

	int  flag=0; 
	bool found=false;
	for (ic=0;ic<TotalCntct;++ic) {
	    ifs1>>ptcl_1>>ptcl_2>>point1_x>>point1_y>>point1_z>>point2_x>>point2_y>>point2_z
		>>radius_1>>radius_2>>penetration>>tang_dispmt>>contact_radius>>R0>>E0>>normal_force>>shear_force;
	    if (ptcl_1==pnumber) {
		found=true;
		if (flag == 0) { // flag == 0 means the first time finding the particle
		    // process particle 1
		    do{
			ifs2 >> ID >> type >> a >> b >> c >> x0 >> y0 >> z0 >> l1 >> m1 >> n1 >> l2 >> m2 >> n2 >> l3 >> m3 >> n3
			     >>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp;
		    } while (ID!=ptcl_1);

		    ofs << "ZONE N="<<NODE<<", E="<<ELEM<<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << endl;

		    ofs << setfill(' ') << setw(16) << cosl(l1)*(-a) + x0
			<< setfill(' ') << setw(16) << cosl(m1)*(-a) + y0
			<< setfill(' ') << setw(16) << cosl(n1)*(-a) + z0
			<< endl;
		    
		    for(i=0, x = -a+0.2*a/3; i < 11; ++i)
		    {
			t = sqrtl(1 - powl(x/a, 2));
			for(j=0, theta = 0; j < 12; ++j, theta += PI/6)
			{
			    y = b*t*sinl(theta);
			    z = c*t*cosl(theta);
			    
			    xp  = cosl(l1)*x + cosl(l2)*y + cosl(l3)*z + x0;
			    yp  = cosl(m1)*x + cosl(m2)*y + cosl(m3)*z + y0;
			    zp  = cosl(n1)*x + cosl(n2)*y + cosl(n3)*z + z0;
			    
			    ofs << setfill(' ') << setw(16) << xp
				<< setfill(' ') << setw(16) << yp
				<< setfill(' ') << setw(16) << zp
				<< endl;
			}
			if (i==0 || i==9)
			    x += 0.2*a/3*2;
			else
			    x += 0.2*a;
		    }
		    
		    ofs << setfill(' ') << setw(16) << cosl(l1)*a + x0
			<< setfill(' ') << setw(16) << cosl(m1)*a + y0
			<< setfill(' ') << setw(16) << cosl(n1)*a + z0
			<< endl;
		    
		    ofs <<
			"	1	1	2	3	\n	"
			"	1	1	3	4	\n	"
			"	1	1	4	5	\n	"
			"	1	1	5	6	\n	"
			"	1	1	6	7	\n	"
			"	1	1	7	8	\n	"
			"	1	1	8	9	\n	"
			"	1	1	9	10	\n	"
			"	1	1	10	11	\n	"
			"	1	1	11	12	\n	"
			"	1	1	12	13	\n	"
			"	1	1	13	2	\n	"
			"	3	2	14	15	\n	"
			"	4	3	15	16	\n	"
			"	5	4	16	17	\n	"
			"	6	5	17	18	\n	"
			"	7	6	18	19	\n	"
			"	8	7	19	20	\n	"
			"	9	8	20	21	\n	"
			"	10	9	21	22	\n	"
			"	11	10	22	23	\n	"
			"	12	11	23	24	\n	"
			"	13	12	24	25	\n	"
			"	2	13	25	14	\n	"
			"	15	14	26	27	\n	"
			"	16	15	27	28	\n	"
			"	17	16	28	29	\n	"
			"	18	17	29	30	\n	"
			"	19	18	30	31	\n	"
			"	20	19	31	32	\n	"
			"	21	20	32	33	\n	"
			"	22	21	33	34	\n	"
			"	23	22	34	35	\n	"
			"	24	23	35	36	\n	"
			"	25	24	36	37	\n	"
			"	14	25	37	26	\n	"
			"	27	26	38	39	\n	"
			"	28	27	39	40	\n	"
			"	29	28	40	41	\n	"
			"	30	29	41	42	\n	"
			"	31	30	42	43	\n	"
			"	32	31	43	44	\n	"
			"	33	32	44	45	\n	"
			"	34	33	45	46	\n	"
			"	35	34	46	47	\n	"
			"	36	35	47	48	\n	"
			"	37	36	48	49	\n	"
			"	26	37	49	38	\n	"
			"	39	38	50	51	\n	"
			"	40	39	51	52	\n	"
			"	41	40	52	53	\n	"
			"	42	41	53	54	\n	"
			"	43	42	54	55	\n	"
			"	44	43	55	56	\n	"
			"	45	44	56	57	\n	"
			"	46	45	57	58	\n	"
			"	47	46	58	59	\n	"
			"	48	47	59	60	\n	"
			"	49	48	60	61	\n	"
			"	38	49	61	50	\n	"
			"	51	50	62	63	\n	"
			"	52	51	63	64	\n	"
			"	53	52	64	65	\n	"
			"	54	53	65	66	\n	"
			"	55	54	66	67	\n	"
			"	56	55	67	68	\n	"
			"	57	56	68	69	\n	"
			"	58	57	69	70	\n	"
			"	59	58	70	71	\n	"
			"	60	59	71	72	\n	"
			"	61	60	72	73	\n	"
			"	50	61	73	62	\n	"
			"	63	62	74	75	\n	"
			"	64	63	75	76	\n	"
			"	65	64	76	77	\n	"
			"	66	65	77	78	\n	"
			"	67	66	78	79	\n	"
			"	68	67	79	80	\n	"
			"	69	68	80	81	\n	"
			"	70	69	81	82	\n	"
			"	71	70	82	83	\n	"
			"	72	71	83	84	\n	"
			"	73	72	84	85	\n	"
			"	62	73	85	74	\n	"
			"	75	74	86	87	\n	"
			"	76	75	87	88	\n	"
			"	77	76	88	89	\n	"
			"	78	77	89	90	\n	"
			"	79	78	90	91	\n	"
			"	80	79	91	92	\n	"
			"	81	80	92	93	\n	"
			"	82	81	93	94	\n	"
			"	83	82	94	95	\n	"
			"	84	83	95	96	\n	"
			"	85	84	96	97	\n	"
			"	74	85	97	86	\n	"
			"	87	86	98	99	\n	"
			"	88	87	99	100	\n	"
			"	89	88	100	101	\n	"
			"	90	89	101	102	\n	"
			"	91	90	102	103	\n	"
			"	92	91	103	104	\n	"
			"	93	92	104	105	\n	"
			"	94	93	105	106	\n	"
			"	95	94	106	107	\n	"
			"	96	95	107	108	\n	"
			"	97	96	108	109	\n	"
			"	86	97	109	98	\n	"
			"	99	98	110	111	\n	"
			"	100	99	111	112	\n	"
			"	101	100	112	113	\n	"
			"	102	101	113	114	\n	"
			"	103	102	114	115	\n	"
			"	104	103	115	116	\n	"
			"	105	104	116	117	\n	"
			"	106	105	117	118	\n	"
			"	107	106	118	119	\n	"
			"	108	107	119	120	\n	"
			"	109	108	120	121	\n	"
			"	98	109	121	110	\n      "
			"	111	110	122	123	\n      "
			"	112	111	123	124	\n      "
			"	113	112	124	125	\n      "
			"	114	113	125	126	\n      "
			"	115	114	126	127	\n      "
			"	116	115	127	128	\n      "
			"	117	116	128	129	\n      "
			"	118	117	129	130	\n      "
			"	119	118	130	131	\n      "
			"	120	119	131	132	\n      "
			"	121	120	132	133	\n      "
			"	110	121	133	122	\n      "
			"	123	122	134	134	\n      "
			"	124	123	134	134	\n      "
			"	125	124	134	134	\n      "
			"	126	125	134	134	\n      "
			"	127	126	134	134	\n      "
			"	128	127	134	134	\n      "
			"	129	128	134	134	\n      "
			"	130	129	134	134	\n      "
			"	131	130	134	134	\n      "
			"	132	131	134	134	\n      "
			"	133	132	134	134	\n      "
			"	122	133	134	134	\n";
		    
		    // process particle 2
		    do{
			ifs2 >> ID >> type >> a >> b >> c >> x0 >> y0 >> z0 >> l1 >> m1 >> n1 >> l2 >> m2 >> n2 >> l3 >> m3 >> n3
			     >>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp;
		    } while (ID!=ptcl_2);

		    ofs << "ZONE N="<<NODE<<", E="<<ELEM<<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL, CONNECTIVITYSHAREZONE = 1" << endl;

		    ofs << setfill(' ') << setw(16) << cosl(l1)*(-a) + x0
			<< setfill(' ') << setw(16) << cosl(m1)*(-a) + y0
			<< setfill(' ') << setw(16) << cosl(n1)*(-a) + z0
			<< endl;
		    
		    for(i=0, x = -a+0.2*a/3; i < 11; ++i)
		    {
			t = sqrtl(1 - powl(x/a, 2));
			for(j=0, theta = 0; j < 12; ++j, theta += PI/6)
			{
			    y = b*t*sinl(theta);
			    z = c*t*cosl(theta);
			    
			    xp  = cosl(l1)*x + cosl(l2)*y + cosl(l3)*z + x0;
			    yp  = cosl(m1)*x + cosl(m2)*y + cosl(m3)*z + y0;
			    zp  = cosl(n1)*x + cosl(n2)*y + cosl(n3)*z + z0;
			    
			    ofs << setfill(' ') << setw(16) << xp
				<< setfill(' ') << setw(16) << yp
				<< setfill(' ') << setw(16) << zp
				<< endl;
			}
			if (i==0 || i==9)
			    x += 0.2*a/3*2;
			else
			    x += 0.2*a;
		    }
		    
		    ofs << setfill(' ') << setw(16) << cosl(l1)*a + x0
			<< setfill(' ') << setw(16) << cosl(m1)*a + y0
			<< setfill(' ') << setw(16) << cosl(n1)*a + z0
			<< endl;

		}
		else { // no need to redraw particle 1, only draw particle 2
		    do{
			ifs2 >> ID >> type >> a >> b >> c >> x0 >> y0 >> z0 >> l1 >> m1 >> n1 >> l2 >> m2 >> n2 >> l3 >> m3 >> n3
			     >>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp;
		    } while (ID!=ptcl_2);

		    ofs << "ZONE N="<<NODE<<", E="<<ELEM<<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL, CONNECTIVITYSHAREZONE = 1" << endl;
		    
		    ofs << setfill(' ') << setw(16) << cosl(l1)*(-a) + x0
			<< setfill(' ') << setw(16) << cosl(m1)*(-a) + y0
			<< setfill(' ') << setw(16) << cosl(n1)*(-a) + z0
			<< endl;
		    
		    for(i=0, x = -a+0.2*a/3; i < 11; ++i)
		    {
			t = sqrtl(1 - powl(x/a, 2));
			for(j=0, theta = 0; j < 12; ++j, theta += PI/6)
			{
			    y = b*t*sinl(theta);
			    z = c*t*cosl(theta);
			    
			    xp  = cosl(l1)*x + cosl(l2)*y + cosl(l3)*z + x0;
			    yp  = cosl(m1)*x + cosl(m2)*y + cosl(m3)*z + y0;
			    zp  = cosl(n1)*x + cosl(n2)*y + cosl(n3)*z + z0;
			    
			    ofs << setfill(' ') << setw(16) << xp
				<< setfill(' ') << setw(16) << yp
				<< setfill(' ') << setw(16) << zp
				<< endl;
			}
			if (i==0 || i==9)
			    x += 0.2*a/3*2;
			else
			    x += 0.2*a;
		    }
		    
		    ofs << setfill(' ') << setw(16) << cosl(l1)*a + x0
			<< setfill(' ') << setw(16) << cosl(m1)*a + y0
			<< setfill(' ') << setw(16) << cosl(n1)*a + z0
			<< endl;
		}

		// draw the two points
		ofs << "ZONE I=2, DATAPACKING=POINT" << endl;
		ofs << setfill(' ') << setw(16) << point1_x
		    << setfill(' ') << setw(16) << point1_y
		    << setfill(' ') << setw(16) << point1_z << endl
		    << setfill(' ') << setw(16) << point2_x
		    << setfill(' ') << setw(16) << point2_y
		    << setfill(' ') << setw(16) << point2_z << endl;

		++flag;

	    }
	    else{
		if (flag!=0)
		    break;
	    }
	}
	
	ifs1.close();
	ifs2.close();
	ofs.close();  

	// Yan: Even if the particle is not found (for example, no contacts at all during some steps), we still need to draw 
	//      the particle for animation.
	if (!found){
	    ifs2.open(filein2);
	    if(!ifs2)  { cout<<"ifs2 stream error!"<<endl; exit(-1);}
	    ofs.open(fileout);
	    if(!ofs)   { cout<<"ofs stream error!"<<endl; exit(-1);}
	    ofs.setf(ios::scientific, ios::floatfield);
	    
	    ifs2 >> TotalNum >> RORS;
	    if(RORS==0)
		ifs2 >> cx >> cy >> cz >> rd >> ht;
	    else
		ifs2 >> cx >> cy >> cz >> wd >> lt >> ht;
	    ifs2>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s
		>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s>>s;	    

	    do{
		ifs2 >> ID >> type >> a >> b >> c >> x0 >> y0 >> z0 >> l1 >> m1 >> n1 >> l2 >> m2 >> n2 >> l3 >> m3 >> n3
		     >>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp;
	    } while (ID!=pnumber);

	    ofs << "ZONE N="<<NODE<<", E="<<ELEM<<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << endl;
	    
	    ofs << setfill(' ') << setw(16) << cosl(l1)*(-a) + x0
		<< setfill(' ') << setw(16) << cosl(m1)*(-a) + y0
		<< setfill(' ') << setw(16) << cosl(n1)*(-a) + z0
		<< endl;
	    
	    for(i=0, x = -a+0.2*a/3; i < 11; ++i)
	    {
		t = sqrtl(1 - powl(x/a, 2));
		for(j=0, theta = 0; j < 12; ++j, theta += PI/6)
		{
		    y = b*t*sinl(theta);
		    z = c*t*cosl(theta);
		    
		    xp  = cosl(l1)*x + cosl(l2)*y + cosl(l3)*z + x0;
		    yp  = cosl(m1)*x + cosl(m2)*y + cosl(m3)*z + y0;
		    zp  = cosl(n1)*x + cosl(n2)*y + cosl(n3)*z + z0;
		    
		    ofs << setfill(' ') << setw(16) << xp
			<< setfill(' ') << setw(16) << yp
			<< setfill(' ') << setw(16) << zp
			<< endl;
		}
		if (i==0 || i==9)
		    x += 0.2*a/3*2;
		else
		    x += 0.2*a;
	    }
	    
	    ofs << setfill(' ') << setw(16) << cosl(l1)*a + x0
		<< setfill(' ') << setw(16) << cosl(m1)*a + y0
		<< setfill(' ') << setw(16) << cosl(n1)*a + z0
		<< endl;
	    
	    ofs <<
		"	1	1	2	3	\n	"
		"	1	1	3	4	\n	"
		"	1	1	4	5	\n	"
		"	1	1	5	6	\n	"
		"	1	1	6	7	\n	"
		"	1	1	7	8	\n	"
		"	1	1	8	9	\n	"
		"	1	1	9	10	\n	"
		"	1	1	10	11	\n	"
		"	1	1	11	12	\n	"
		"	1	1	12	13	\n	"
		"	1	1	13	2	\n	"
		"	3	2	14	15	\n	"
		"	4	3	15	16	\n	"
		"	5	4	16	17	\n	"
		"	6	5	17	18	\n	"
		"	7	6	18	19	\n	"
		"	8	7	19	20	\n	"
		"	9	8	20	21	\n	"
		"	10	9	21	22	\n	"
		"	11	10	22	23	\n	"
		"	12	11	23	24	\n	"
		"	13	12	24	25	\n	"
		"	2	13	25	14	\n	"
		"	15	14	26	27	\n	"
		"	16	15	27	28	\n	"
		"	17	16	28	29	\n	"
		"	18	17	29	30	\n	"
		"	19	18	30	31	\n	"
		"	20	19	31	32	\n	"
		"	21	20	32	33	\n	"
		"	22	21	33	34	\n	"
		"	23	22	34	35	\n	"
		"	24	23	35	36	\n	"
		"	25	24	36	37	\n	"
		"	14	25	37	26	\n	"
		"	27	26	38	39	\n	"
		"	28	27	39	40	\n	"
		"	29	28	40	41	\n	"
		"	30	29	41	42	\n	"
		"	31	30	42	43	\n	"
		"	32	31	43	44	\n	"
		"	33	32	44	45	\n	"
		"	34	33	45	46	\n	"
		"	35	34	46	47	\n	"
		"	36	35	47	48	\n	"
		"	37	36	48	49	\n	"
		"	26	37	49	38	\n	"
		"	39	38	50	51	\n	"
		"	40	39	51	52	\n	"
		"	41	40	52	53	\n	"
		"	42	41	53	54	\n	"
		"	43	42	54	55	\n	"
		"	44	43	55	56	\n	"
		"	45	44	56	57	\n	"
		"	46	45	57	58	\n	"
		"	47	46	58	59	\n	"
		"	48	47	59	60	\n	"
		"	49	48	60	61	\n	"
		"	38	49	61	50	\n	"
		"	51	50	62	63	\n	"
		"	52	51	63	64	\n	"
		"	53	52	64	65	\n	"
		"	54	53	65	66	\n	"
		"	55	54	66	67	\n	"
		"	56	55	67	68	\n	"
		"	57	56	68	69	\n	"
		"	58	57	69	70	\n	"
		"	59	58	70	71	\n	"
		"	60	59	71	72	\n	"
		"	61	60	72	73	\n	"
		"	50	61	73	62	\n	"
		"	63	62	74	75	\n	"
		"	64	63	75	76	\n	"
		"	65	64	76	77	\n	"
		"	66	65	77	78	\n	"
		"	67	66	78	79	\n	"
		"	68	67	79	80	\n	"
		"	69	68	80	81	\n	"
		"	70	69	81	82	\n	"
		"	71	70	82	83	\n	"
		"	72	71	83	84	\n	"
		"	73	72	84	85	\n	"
		"	62	73	85	74	\n	"
		"	75	74	86	87	\n	"
		"	76	75	87	88	\n	"
		"	77	76	88	89	\n	"
		"	78	77	89	90	\n	"
		"	79	78	90	91	\n	"
		"	80	79	91	92	\n	"
		"	81	80	92	93	\n	"
		"	82	81	93	94	\n	"
		"	83	82	94	95	\n	"
		"	84	83	95	96	\n	"
		"	85	84	96	97	\n	"
		"	74	85	97	86	\n	"
		"	87	86	98	99	\n	"
		"	88	87	99	100	\n	"
		"	89	88	100	101	\n	"
		"	90	89	101	102	\n	"
		"	91	90	102	103	\n	"
		"	92	91	103	104	\n	"
		"	93	92	104	105	\n	"
		"	94	93	105	106	\n	"
		"	95	94	106	107	\n	"
		"	96	95	107	108	\n	"
		"	97	96	108	109	\n	"
		"	86	97	109	98	\n	"
		"	99	98	110	111	\n	"
		"	100	99	111	112	\n	"
		"	101	100	112	113	\n	"
		"	102	101	113	114	\n	"
		"	103	102	114	115	\n	"
		"	104	103	115	116	\n	"
		"	105	104	116	117	\n	"
		"	106	105	117	118	\n	"
		"	107	106	118	119	\n	"
		"	108	107	119	120	\n	"
		"	109	108	120	121	\n	"
		"	98	109	121	110	\n      "
		"	111	110	122	123	\n      "
		"	112	111	123	124	\n      "
		"	113	112	124	125	\n      "
		"	114	113	125	126	\n      "
		"	115	114	126	127	\n      "
		"	116	115	127	128	\n      "
		"	117	116	128	129	\n      "
		"	118	117	129	130	\n      "
		"	119	118	130	131	\n      "
		"	120	119	131	132	\n      "
		"	121	120	132	133	\n      "
		"	110	121	133	122	\n      "
		"	123	122	134	134	\n      "
		"	124	123	134	134	\n      "
		"	125	124	134	134	\n      "
		"	126	125	134	134	\n      "
		"	127	126	134	134	\n      "
		"	128	127	134	134	\n      "
		"	129	128	134	134	\n      "
		"	130	129	134	134	\n      "
		"	131	130	134	134	\n      "
		"	132	131	134	134	\n      "
		"	133	132	134	134	\n      "
		"	122	133	134	134	\n";

	    ifs2.close();
	    ofs.close();
	}
    }

    return 0;
}
    
