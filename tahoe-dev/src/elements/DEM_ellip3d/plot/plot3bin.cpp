#include "TECIO.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cmath>

#ifndef NULL
#define NULL 0
#endif

const long double PI = 3.1415927;
const std::size_t OWID  = 15;
const std::size_t OPREC = 6;

const int NODE  = 134;
const int ELEM  = 144;
const int CONNECT[ELEM][4] = {
  {1	,	1	,	2	,	3	},
  {1	,	1	,	3	,	4	},
  {1	,	1	,	4	,	5	},
  {1	,	1	,	5	,	6	},
  {1	,	1	,	6	,	7	},
  {1	,	1	,	7	,	8	},
  {1	,	1	,	8	,	9	},
  {1	,	1	,	9	,	10	},
  {1	,	1	,	10	,	11	},
  {1	,	1	,	11	,	12	},
  {1	,	1	,	12	,	13	},
  {1	,	1	,	13	,	2	},
  {3	,	2	,	14	,	15	},
  {4	,	3	,	15	,	16	},
  {5	,	4	,	16	,	17	},
  {6	,	5	,	17	,	18	},
  {7	,	6	,	18	,	19	},
  {8	,	7	,	19	,	20	},
  {9	,	8	,	20	,	21	},
  {10	,	9	,	21	,	22	},
  {11	,	10	,	22	,	23	},
  {12	,	11	,	23	,	24	},
  {13	,	12	,	24	,	25	},
  {2	,	13	,	25	,	14	},
  {15	,	14	,	26	,	27	},
  {16	,	15	,	27	,	28	},
  {17	,	16	,	28	,	29	},
  {18	,	17	,	29	,	30	},
  {19	,	18	,	30	,	31	},
  {20	,	19	,	31	,	32	},
  {21	,	20	,	32	,	33	},
  {22	,	21	,	33	,	34	},
  {23	,	22	,	34	,	35	},
  {24	,	23	,	35	,	36	},
  {25	,	24	,	36	,	37	},
  {14	,	25	,	37	,	26	},
  {27	,	26	,	38	,	39	},
  {28	,	27	,	39	,	40	},
  {29	,	28	,	40	,	41	},
  {30	,	29	,	41	,	42	},
  {31	,	30	,	42	,	43	},
  {32	,	31	,	43	,	44	},
  {33	,	32	,	44	,	45	},
  {34	,	33	,	45	,	46	},
  {35	,	34	,	46	,	47	},
  {36	,	35	,	47	,	48	},
  {37	,	36	,	48	,	49	},
  {26	,	37	,	49	,	38	},
  {39	,	38	,	50	,	51	},
  {40	,	39	,	51	,	52	},
  {41	,	40	,	52	,	53	},
  {42	,	41	,	53	,	54	},
  {43	,	42	,	54	,	55	},
  {44	,	43	,	55	,	56	},
  {45	,	44	,	56	,	57	},
  {46	,	45	,	57	,	58	},
  {47	,	46	,	58	,	59	},
  {48	,	47	,	59	,	60	},
  {49	,	48	,	60	,	61	},
  {38	,	49	,	61	,	50	},
  {51	,	50	,	62	,	63	},
  {52	,	51	,	63	,	64	},
  {53	,	52	,	64	,	65	},
  {54	,	53	,	65	,	66	},
  {55	,	54	,	66	,	67	},
  {56	,	55	,	67	,	68	},
  {57	,	56	,	68	,	69	},
  {58	,	57	,	69	,	70	},
  {59	,	58	,	70	,	71	},
  {60	,	59	,	71	,	72	},
  {61	,	60	,	72	,	73	},
  {50	,	61	,	73	,	62	},
  {63	,	62	,	74	,	75	},
  {64	,	63	,	75	,	76	},
  {65	,	64	,	76	,	77	},
  {66	,	65	,	77	,	78	},
  {67	,	66	,	78	,	79	},
  {68	,	67	,	79	,	80	},
  {69	,	68	,	80	,	81	},
  {70	,	69	,	81	,	82	},
  {71	,	70	,	82	,	83	},
  {72	,	71	,	83	,	84	},
  {73	,	72	,	84	,	85	},
  {62	,	73	,	85	,	74	},
  {75	,	74	,	86	,	87	},
  {76	,	75	,	87	,	88	},
  {77	,	76	,	88	,	89	},
  {78	,	77	,	89	,	90	},
  {79	,	78	,	90	,	91	},
  {80	,	79	,	91	,	92	},
  {81	,	80	,	92	,	93	},
  {82	,	81	,	93	,	94	},
  {83	,	82	,	94	,	95	},
  {84	,	83	,	95	,	96	},
  {85	,	84	,	96	,	97	},
  {74	,	85	,	97	,	86	},
  {87	,	86	,	98	,	99	},
  {88	,	87	,	99	,	100	},
  {89	,	88	,	100	,	101	},
  {90	,	89	,	101	,	102	},
  {91	,	90	,	102	,	103	},
  {92	,	91	,	103	,	104	},
  {93	,	92	,	104	,	105	},
  {94	,	93	,	105	,	106	},
  {95	,	94	,	106	,	107	},
  {96	,	95	,	107	,	108	},
  {97	,	96	,	108	,	109	},
  {86	,	97	,	109	,	98	},
  {99	,	98	,	110	,	111	},
  {100	,	99	,	111	,	112	},
  {101	,	100	,	112	,	113	},
  {102	,	101	,	113	,	114	},
  {103	,	102	,	114	,	115	},
  {104	,	103	,	115	,	116	},
  {105	,	104	,	116	,	117	},
  {106	,	105	,	117	,	118	},
  {107	,	106	,	118	,	119	},
  {108	,	107	,	119	,	120	},
  {109	,	108	,	120	,	121	},
  {98	,	109	,	121	,	110     },
  {111	,	110	,	122	,	123     },
  {112	,	111	,	123	,	124     },
  {113	,	112	,	124	,	125     },
  {114	,	113	,	125	,	126     },
  {115	,	114	,	126	,	127     },
  {116	,	115	,	127	,	128     },
  {117	,	116	,	128	,	129     },
  {118	,	117	,	129	,	130     },
  {119	,	118	,	130	,	131     },
  {120	,	119	,	131	,	132     },
  {121	,	120	,	132	,	133     },
  {110	,	121	,	133	,	122     },
  {123	,	122	,	134	,	134     },
  {124	,	123	,	134	,	134     },
  {125	,	124	,	134	,	134     },
  {126	,	125	,	134	,	134     },
  {127	,	126	,	134	,	134     },
  {128	,	127	,	134	,	134     },
  {129	,	128	,	134	,	134     },
  {130	,	129	,	134	,	134     },
  {131	,	130	,	134	,	134     },
  {132	,	131	,	134	,	134     },
  {133	,	132	,	134	,	134     },
  {122	,	133	,	134	,	134     }  };

int main(int argc, char *argv[]) 
{
  if(argc < 2) {
    std::cout << std::endl 
	      << "-- Plot All Particles --" << std::endl
	      << "Usage:" << std::endl
	      << "1) process a single file:  plot particle_file" << std::endl
	      << "   --example: plot iso_particle_008" << std::endl
	      << "2) process multiple files: plot particle_file_prefix  first_suffix  last_suffix  suffix_increment" << std::endl
	      << "   --example: plot iso_particle  0  100  5" << std::endl << std::endl;
    return -1;
  }	

  int first, last, incre;
  if(argc == 2) {
    first = 0;
    last  = 1;
    incre = 2;
  }
  else {
    first = atoi(argv[2]);
    last  = atoi(argv[3]);
    incre = atoi(argv[4]);
  }

  std::ifstream ifs;
  char filein[50];
  char fileout[50];
  char num[4], str[20];

  for (int snapshot = first; snapshot <= last; snapshot += incre) {
    if(argc == 2)
      strcpy(filein, argv[1]);
    else {
      sprintf(num, "%03d", snapshot);
      strcpy(filein, argv[1]);
      strcat(filein, "_");
      strcat(filein, num);
    }

    INTEGER4 Debug      = 0;
    INTEGER4 VIsDouble  = 1; // 0=Single, 1=Double
    INTEGER4 FileType   = 0; // 0=Full, 1=Grid, 2=Solution
    INTEGER4 FileFormat = 0; // 0 == PLT, 1 == SZPLT
    strcpy(fileout, filein);
    if (FileFormat == 0)
      strcat(fileout, ".plt");
    else if (FileFormat == 1)
      strcat(fileout, ".szplt");
    std::cout << "generating file " << fileout << " ......" <<std::endl;
    TECINI142((char*)"FEQuadrilateral",
	      (char*)"x y z",
	      (char*)fileout,
	      (char*)".",
	      &FileFormat,
	      &FileType,
	      &Debug,
	      &VIsDouble);

    ifs.open(filein);
    if(!ifs)  { std::cout<<"stream error!"<<std::endl; exit(-1);}

    long int totalParticle;
    ifs >> totalParticle;
    ifs >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str 
	>> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str;

    long int id;
    int type;
    long double cx, cy, cz, rd, wd, lt, ht;
    long double a, b, c, x0, y0, z0, l1, l2, l3, m1, m2, m3, n1, n2, n3, tmp;
    long double x, y, z, xp, yp, zp, tt, theta; // x, y, z are local coodinates, xp, yp, zp are global

    INTEGER4 ZoneType   = 3;     /* set the zone type to
                                * FEQuadrilateral
                                */
    INTEGER4 NumPts     = 0;
    INTEGER4 NumElems   = 0;
    INTEGER4 NumFaces   = 0;
    INTEGER4 ICellMax   = 0;     /* not used */
    INTEGER4 JCellMax   = 0;     /* not used */
    INTEGER4 KCellMax   = 0;     /* not used */
    double   SolTime    = 0;
    INTEGER4 StrandID   = 0;     /* StaticZone */
    INTEGER4 ParentZn   = 0;
    INTEGER4 IsBlock    = 1;     /* Block */
    INTEGER4 NFConns    = 0;     /* Specify the number of Face
                                * Neighbor Connections in the
                                * Zone.  When this value is
                                * greater than zero, TECFACE must
                                * be called prior to creating the
                                * next zone or ending the file.
                                */
    NumPts   = NODE * totalParticle;
    NumElems = ELEM * totalParticle;
    SolTime  = snapshot;
    StrandID = 1;

    INTEGER4 FNMode                       = 2;

    INTEGER4 TotalNumFaceNodes            = 1;  /* Not used for
                                               * FEQuad zones*/
    INTEGER4 NumConnectedBoundaryFaces    = 1;  /* Not used for
                                               * FEQuad zones*/
    INTEGER4 TotalNumBoundaryConnections  = 1;  /* Not used for
                                               * FEQuad zones*/
    INTEGER4 ShrConn                      = 0;
    INTEGER4 ValueLocation[3] = {1, 1, 1};  /* Specify the variable
                                         * values at the nodes.
                                         * NOTE:  Because all of
                                         * the variables are
                                         * defined at the nodes,
                                         * we can just pass
                                         * NULL for this array.
                                         * We are providing the
                                         * array for illustration
                                         * purposes.
                                         */
    TECZNE142((char*)"ptcl",
	      &ZoneType,
	      &NumPts,
	      &NumElems,
	      &NumFaces,
	      &ICellMax,
	      &JCellMax,
	      &KCellMax,
	      &SolTime,
	      &StrandID,
	      &ParentZn,
	      &IsBlock,
	      &NFConns,
	      &FNMode,
	      &TotalNumFaceNodes,
	      &NumConnectedBoundaryFaces,
	      &TotalNumBoundaryConnections,
	      NULL,
	      ValueLocation,
	      NULL,
	      &ShrConn);    

    std::vector<double> coordX;
    std::vector<double> coordY;
    std::vector<double> coordZ;
    for (long int numParticle = 0; numParticle < totalParticle; ++numParticle) {
      ifs >> id >> type >> a >> b >> c >> x0 >> y0 >> z0 >> l1 >> m1 >> n1 >> l2 >> m2 >> n2 >> l3 >> m3 >> n3
	  >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp;

      std::vector<long double> point;
      coordX.push_back(cosl(l1)*(-a) + x0);
      coordY.push_back(cosl(m1)*(-a) + y0);
      coordZ.push_back(cosl(n1)*(-a) + z0);

      // C++ requires multiple definitions in for loop must be the same type;
      // if they are moved into the loop, it and jt conflict with theta type, leading to errors.
      // so it and jt must be kept out of the loop.
      int it, jt; 
      for (it = 0, x = -a+0.2*a/3; it < 11; ++it) {
	tt = sqrtl(1 - powl(x/a, 2));
	for (jt = 0, theta = 0; jt < 12; ++jt, theta += PI/6) {
	  y = b * tt * sinl(theta);
	  z = c * tt * cosl(theta);
		    
	  xp  = cosl(l1)*x + cosl(l2)*y + cosl(l3)*z + x0;
	  yp  = cosl(m1)*x + cosl(m2)*y + cosl(m3)*z + y0;
	  zp  = cosl(n1)*x + cosl(n2)*y + cosl(n3)*z + z0;

	  coordX.push_back(xp);
	  coordY.push_back(yp);
	  coordZ.push_back(zp);

	}
	if (it == 0 || it == 9) // two ends
	  x += 0.2*a/3*2;
	else
	  x += 0.2*a;
      }

      coordX.push_back(cosl(l1)*a + x0);
      coordY.push_back(cosl(m1)*a + y0);
      coordZ.push_back(cosl(n1)*a + z0);
    } // end of for numParticle

    INTEGER4 DIsDouble =  1;  /* Set DIsDouble to zero to use
                             * variables in float format.
                             */

    /* Call TECDAT once for each variable */
    TECDAT142(&NumPts, &coordX[0], &DIsDouble);
    TECDAT142(&NumPts, &coordY[0], &DIsDouble);
    TECDAT142(&NumPts, &coordZ[0], &DIsDouble);

    ///*
    std::vector<int> ConnList(NumElems*4);
    int index = 0;
    for (int kt = 0; kt < totalParticle; ++kt) {
      for (int it = 0; it < ELEM; ++it) {
	ConnList[index++] = CONNECT[it][0] + NODE*kt;
	ConnList[index++] = CONNECT[it][1] + NODE*kt;
	ConnList[index++] = CONNECT[it][2] + NODE*kt;
	ConnList[index++] = CONNECT[it][3] + NODE*kt;
      }
    }
    TECNOD142(&ConnList[0]);
    //*/

    /*
    INTEGER4 *ConnList = new INTEGER4[NumElems*4];
    int index = 0;
    for (int kt = 0; kt < totalParticle; ++kt) {
      std::cout << " ptcl=" << kt << " index=" << index;
      for (int it = 0; it < ELEM; ++it) {
	ConnList[index++] = CONNECT[it][0] + NODE*kt;
	ConnList[index++] = CONNECT[it][1] + NODE*kt;
	ConnList[index++] = CONNECT[it][2] + NODE*kt;
	ConnList[index++] = CONNECT[it][3] + NODE*kt;
      }
      std::cout << " " << index << std::endl;
    }
    TECNOD142(ConnList);
    delete []ConnList;
    */

    TECEND142();
    ifs.close();

  }
    
  return 0;
}
