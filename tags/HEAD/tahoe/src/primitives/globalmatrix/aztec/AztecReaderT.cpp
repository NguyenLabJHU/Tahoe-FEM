/* $Id: AztecReaderT.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (08/12/1998)                                          */
/* utility to read Aztec options and parameters                           */

#include "AztecReaderT.h"

/* library support options */
#ifdef __AZTEC__

#include <iostream.h>
#include <ctype.h>
#include <string.h>

#include "ExceptionCodes.h"
#include "az_aztec_defs.h"
#include "fstreamT.h"

/* class parameters */
const int kNameLength = 255;

/* option constants and indices */
const int kNumOptions = 13; // AZ_OPTIONS_SIZE is 15, but only support 13
const char* opt_name[kNumOptions] = {
"AZ_solver",
"AZ_scaling",
"AZ_precond",
"AZ_conv",
"AZ_output",
"AZ_pre_calc",
"AZ_max_iter",
"AZ_poly_ord",
"AZ_overlap",
"AZ_kspace",
"AZ_orthog",
"AZ_aux_vec",
"AZ_print_freq"};
const int opt_idex[kNumOptions] = {0,1,2,3,4,5,6,7,8,9,10,11,12};

/* parameter constants and indices */
const int   kNumParams = 2; // AZ_PARAMS_SIZE is 5, but only support 2
const char* param_name[kNumParams] = {"AZ_tol","AZ_drop"};
const int   param_idex[kNumParams] = {0,2};

/* constructor */
AztecReaderT::AztecReaderT(void) { }
	
/* read option and return integer value */
void AztecReaderT::ReadOption(ifstreamT& in, int& dex, int& value)
{
	/* next character */
	char next_char = in.next_char();
	
	/* read parameter index */
	if (isdigit(next_char))
		in >> dex;
	else
	{
		char name[kNameLength];

		/* read up to next white space or "=" */
		int i = 0;
		in.get(next_char);
		while (!isspace(next_char) && next_char != '=')
		{
			if (i >= kNameLength - 1)
			{
				/* terminate */
				name[i] = '\0';
				cout << "\n AztecReaderT::ReadOption: error reading option name: "
				     << name << endl;
				throw eBadInputValue;	
			}

			name[i++] = next_char;
			in.get(next_char);
		}
		
		/* terminate */
		name[i] = '\0';
		
		/* putback "=" */
		if (next_char == '=') in.putback(next_char);

		/* convert to index */
		dex = OptionNameToIndex(name);
	}

	/* clear "=" */
	in.get(next_char);
	while (next_char != '=') in.get(next_char);
		
	/* read parameter value */
	in >> value;
}

/* read parameter and return double value */
void AztecReaderT::ReadParameter(ifstreamT& in, int& dex, double& value)
{
	/* next character */
	char next_char = in.next_char();
	
	/* read parameter index */
	if (isdigit(next_char))
		in >> dex;
	else
	{
		char name[kNameLength];

		/* read up to next white space or "=" */
		int i = 0;
		in.get(next_char);
		while (!isspace(next_char) && next_char != '=')
		{
			if (i >= kNameLength - 1)
			{
				/* terminate */
				name[i] = '\0';
				cout << "\n AztecReaderT::ReadParameter: error reading parameter name: "
				     << name << endl;
				throw eBadInputValue;	
			}
		
			name[i++] = next_char;
			in.get(next_char);
		}
		
		/* terminate */
		name[i] = '\0';
		
		/* putback "=" */
		if (next_char == '=') in.putback(next_char);

		/* convert to index */
		dex = ParamNameToIndex(name);
	}

	/* clear "=" */
	in.get(next_char);
	while (next_char != '=') in.get(next_char);
	
	/* read parameter value */
	in >> value;
}

/*************************************************************************
* Private
*************************************************************************/

/* option name to index conversion */
int AztecReaderT::OptionNameToIndex(const char* name) const
{
	int i, found = 0;
	for (i = 0; i < kNumOptions && !found; i++)
		if (strcmp(name,opt_name[i]) == 0) found = 1;

	if (found)
		return opt_idex[i-1];
	else
	{
		cout << "\n AztecReaderT: unknown option name: " << name << endl;
		throw eBadInputValue;
		return -1;
	}
}

/* parameter name to index conversion */
int AztecReaderT::ParamNameToIndex(const char* name) const
{
	int i, found = 0;
	for (i = 0; i < kNumParams && !found; i++)
		if (strcmp(name,param_name[i]) == 0) found = 1;

	if (found)
		return param_idex[i-1];
	else
	{
		cout << "\n AztecReaderT: unknown parameter name: " << name << endl;
		throw eBadInputValue;
		return -1;
	}
}

/* library support options */
#endif /* __AZTEC__ */
