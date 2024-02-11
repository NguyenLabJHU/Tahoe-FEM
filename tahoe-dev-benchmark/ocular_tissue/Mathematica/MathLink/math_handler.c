/* init_mathematica_compare.c

AUTHOR:  Brian M. Adams (briadam@sandia.gov) 
CREATED: 10/13/2005

DESCRIPTION:  C program to 

   1.  read params.in (or other specified DAKOTA params file)
   2.  launch mathematica kernel and call mathematica code for
       fn_val = compare(x1,x2)
   3.  write results to results.out (or other specified outfile)


DISCLAIMER:  THERE IS LITTLE ERROR CHECKING!

INPUTS:  	command-line parameters only
		argv[0] -- executable as called
 		argv[1] -- filename containing params written by DAKOTA
				 (see below for file format)
		argv[2] -- filename to write results for return to DAKOTA

OUTPUTS:	Writes function value to argv[2]

SAMPLE COMPILER DIRECTIVE:  gcc init_mathematica_compare.c -I /usr/local/Wolfram/Mathematica/5.01/AddOns/MathLink/DeveloperKit/Linux/CompilerAdditions -L /usr/local/Wolfram/Mathematica/5.01/AddOns/MathLink/DeveloperKit/Linux/CompilerAdditions -l ML -o init_mathematica_compare

or 

mcc init_mathematica_compare.c -o init_mathematica_compare
*/

#include<stdio.h>
#include"mathlink.h"


/* Number of design variables and length of file I/O buffer */
#define SCRIPT_NAME "compare"
#define NUM_VARS 6			
#define MAX_LINE_LENGTH 512
#define REAL_LIST false
#define DEBUG

int main(int argc, char * argv[ ])
{

	int i;
	
	/* File I/O */
	FILE *infile;
	FILE *outfile;
	char discard_string[MAX_LINE_LENGTH];
	char *line_ptr = NULL;
	int n = 0;
	
	/* input parameters */
	double x[NUM_VARS];
	
	/* Define the connection parameters to the Mathematica kernel, assuming
	   we'll instantiate Mathematica at every fn. eval. */
	int num_strings = 4;
	char *Mparams[ ]= { "-linkname",
			    "math -mathlink",
			    "-linkmode",
			    "launch",
			    NULL };
	
	MLINK link;
	MLEnvironment env;
	
	/* Open input file for reading parameters */
	infile = fopen(argv[1],"r"); 
	if ( infile == NULL ) {
		printf("Error -- couldn't open file %s\n",argv[1]);
		return(1);
	}
	
	/* Discard the first line -- could read number vars */
	getline(&line_ptr, &n, infile);
	
	/* Read value, variable name for x1, x2 */	
	for ( i=0; i<NUM_VARS; i++ ) {
		fscanf( infile,"%le", &x[i]);
		fscanf( infile,"%s", discard_string );
	}
	
	fclose(infile);

#ifdef DEBUG
	printf(SCRIPT_NAME " : ");
	printf("{ %f", x[0]);
	for ( i=1; i<NUM_VARS; i++ ) {
		printf(", %f",x[i]);
	}
	printf("} -> ");
#endif

	/* Now initialize and call Mathematica */
	env = MLInitialize(NULL);
	if (env == NULL) return 1;

	link = MLOpen(num_strings, Mparams);
	if (link == NULL) return 1;

	MLActivate(link);

	/* Load the compare function definition into Mathematica */
	MLPutFunction(link, "EnterTextPacket", 1);
	MLPutString(link, "<<" SCRIPT_NAME ".m");
	MLEndPacket(link);

	/* call the compare function with x[1], x[2]*/
  	MLPutFunction(link, "EvaluatePacket", 1);
#if REAL_LIST
//	MLPutFunction(link, "Last", 1);
  	MLPutFunction(link, SCRIPT_NAME, 1);
	MLPutRealList(link,x,NUM_VARS);
#else
	MLPutFunction(link, SCRIPT_NAME, NUM_VARS);
  	for ( i=0; i<NUM_VARS; i++ ) { MLPutDouble(link, x[i]); }
#endif
	MLEndPacket(link);

	/* skip any packets before the first ReturnPacket 
	   and then get the result from the calculation   */
	int pkt;
	while( (pkt = MLNextPacket(link), pkt) && pkt != RETURNPKT) {
#ifdef DEBUG
//		printf("pkt %d\n",pkt);
#endif
		MLNewPacket(link);
		if (MLError(link)) {
			fprintf( stderr, "2 Error detected by MathLink: %s.\n",
			MLErrorMessage(link)); }
	}
	
#ifdef DEBUG
//		printf(">pkt %d\n",pkt);
#endif
	double fn_val=0.0;
	MLGetDouble(link, &fn_val);
	if (MLError(link)) {
		fprintf( stderr, "3 Error detected by MathLink: %s.\n",
		MLErrorMessage(link)); }
	MLClose(link);
	MLDeinitialize(env);
	
#ifdef DEBUG
	printf("f: %f\n",fn_val);
#endif

	/* Write results.out.X */
	outfile = fopen(argv[2],"w");   // open results.out.X to get parameters
	if ( outfile == NULL ) {
		printf("Error\n");
		return(1);
	}

	fprintf(outfile, "%18.10le  f\n", fn_val);
	fclose(outfile);

	return(0);

}
