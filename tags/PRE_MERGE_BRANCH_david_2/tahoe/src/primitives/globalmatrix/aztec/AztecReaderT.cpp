/* $Id: AztecReaderT.cpp,v 1.7 2004-07-19 04:16:45 paklein Exp $ */
/* created: paklein (08/12/1998) */
#include "AztecReaderT.h"

/* library support options */
#ifdef __AZTEC__

#include <iostream.h>
#include <ctype.h>
#include <string.h>

#include "ExceptionT.h"
#include "az_aztec_defs.h"
#include "ifstreamT.h"

using namespace Tahoe;

/* class parameters */
const int kNameLength = 255;

/* option constants and indices */
const int kNumOptions = 26; // AZ_OPTIONS_SIZE is 47, but only support some
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
"AZ_type_overlap",
"AZ_kspace",
"AZ_orthog",
"AZ_aux_vec",
"AZ_reorder",
"AZ_keep_info",
"AZ_recursion_level",
"AZ_print_freq",
"AZ_graph_fill",
"AZ_subdomain_solve",
"AZ_init_guess",
"AZ_keep_kvecs",
"AZ_apply_kvecs",
"AZ_orth_kvecs",
"AZ_ignore_scaling",
"AZ_check_update_size",
"AZ_extreme"
};
const int opt_idex[kNumOptions] = {
 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
20, 21, 22, 23, 24, 25};

/* parameter constants and indices */
const int   kNumParams = 3; // AZ_PARAMS_SIZE is 30, but only support some
const char* param_name[kNumParams] = {"AZ_tol", "AZ_drop", "AZ_ilut_fill"};
const int   param_idex[kNumParams] = {0, 1, 2};

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
				throw ExceptionT::kBadInputValue;	
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
				throw ExceptionT::kBadInputValue;	
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
		throw ExceptionT::kBadInputValue;
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
		throw ExceptionT::kBadInputValue;
		return -1;
	}
}

/* library support options */
#endif /* __AZTEC__ */
