/* $Id: UpLagr_ExternalFieldT.cpp,v 1.9 2004-06-28 22:41:14 hspark Exp $ */

#include "UpLagr_ExternalFieldT.h"

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "ExodusT.h"

using namespace Tahoe;

/* constructor */
UpLagr_ExternalFieldT::UpLagr_ExternalFieldT(const ElementSupportT& support, const FieldT& field):
	UpdatedLagrangianT(support, field),
	fLocExternalField(LocalArrayT::kUnspecified)
{

}

/* initialization */
void UpLagr_ExternalFieldT::Initialize(void)
{
	/* inherited */
	UpdatedLagrangianT::Initialize();

	/* streams */
	ifstreamT&  in = ElementSupport().Input();
	ofstreamT& out = ElementSupport().Output();
	
	/* external file */
	in >> fExternalFieldFormat;
	in >> fExternalFieldFile;

	/* echo */
	out << " External field file format. . . . . . . . . . . = " << fExternalFieldFormat << '\n';
	out << " External field file . . . . . . . . . . . . . . = " << fExternalFieldFile << '\n';
	if (fExternalFieldFormat != IOBaseT::kExodusII) {
		cout << "\n UpLagr_ExternalFieldT::Initialize: external field file format must be ExodusII" 
		     << endl;
		throw ExceptionT::kBadInputValue;
	}

	/* read variable labels */
	int num_variables = -1;
	in >> num_variables; if (num_variables < 0) throw ExceptionT::kBadInputValue;
	fExternalFieldLabels.Dimension(num_variables);
	for (int i = 0; i < num_variables; i++)
		in >> fExternalFieldLabels[i];
	
	/* build relative path to external file */
	fExternalFieldFile.ToNativePathName();
	StringT path;
	path.FilePath(in.filename());
	fExternalFieldFile.Prepend(path);
	
	/* verify file */
	ExodusT exo(cout);
	if (!exo.OpenRead(fExternalFieldFile)) 
	{
		cout << "\n UpLagr_ExternalFieldT::Initialize: error opening external field file: " 
		     << fExternalFieldFile << endl;
		throw ExceptionT::kBadInputValue;
	}
	else /* resolve index of variables in database */
	{
		/* read nodal variable labels stored in file */
		ArrayT<StringT> labels;
		exo.ReadLabels(labels, ExodusT::kNode);
		
		/* index list */
		fFieldVariableIndex.Dimension(fExternalFieldLabels.Length());
		fFieldVariableIndex = -1;

		/* resolve location of each field component */
		for (int i = 0; i < fExternalFieldLabels.Length(); i++)
		{
			bool found = false;
			for (int j = 0; j < labels.Length() && !found; j++)
				if (labels[j] == fExternalFieldLabels[i])
				{
					fFieldVariableIndex[i] = j;
					found = true;					
				}
				
			/* verify */
			if (!found)
			{
				cout << "\n UpLagr_ExternalFieldT::Initialize: variable not found in\n"
				     <<   "     external database: " << fExternalFieldLabels[i] << endl;
				throw ExceptionT::kBadInputValue;
			}
		}
	}

	/* echo labels and indicies */
	out << " Number of external field components . . . . . . = " << fExternalFieldLabels.Length() << '\n';
	out << setw(25) << "variable" << setw(kIntWidth) << "index" << '\n';
	for (int j = 0; j < fExternalFieldLabels.Length(); j++)
		out << setw(25)        << fExternalFieldLabels[j] 
		    << setw(kIntWidth) <<  fFieldVariableIndex[j]+1 << '\n';

	/* read list of times in the database */
	dArrayT time_steps(exo.NumTimeSteps());
	for (int i = 0; i < time_steps.Length(); i++)
		exo.ReadTime(i+1, time_steps[i]);
	fTimeSteps.SetValues(time_steps);

	/* get node number map */
	fNodeMap.Dimension(exo.NumNodes());
	exo.ReadNodeMap(fNodeMap);
	fNodeMap--; /* offset to internal numbering */

	/* allocate global assembly array */
	fExternalField.Dimension(ElementSupport().NumNodes(), fExternalFieldLabels.Length());
//NOTE: This is going to be a big waste of space if this element group
//      only involves a small fraction of the total number of nodes in
//      the model, but it makes it much easier to retrieve the values
//      like other field variables.

	/* allocate local array */
	fLocExternalField.Dimension(NumElementNodes(), fExternalFieldLabels.Length());
	fLocExternalField.SetGlobal(fExternalField);
	
	/* allocate space for reading nodal values */
	fNodalValues.Dimension(exo.NumNodes());
}

/* interpolate external field to the current time step */
void UpLagr_ExternalFieldT::InitStep(void)
{
	/* inherited */
	UpdatedLagrangianT::InitStep();
	
	/* interval containing current time */
	double time = ElementSupport().Time();
	int i_1 = fTimeSteps.Range(time);
	
	/* interpolate database to current time */
	fExternalField = 0.0;
	
	/* interpolation: [1 - beta] u(i_1) + beta u(i_2) */
	int i_2;
	double beta; 

	/* set interpolation parameters */
	if (i_1 == 0) /* before first set in database */
	{
		i_2 = i_1;
		beta = 1.0;	
	}
	else if (i_1 == fTimeSteps.Length()) /* beyond last interval */
	{
		i_2 = --i_1;
		beta = 0.0;	
	}
	else /* interpolate */
	{
		i_1--;
		i_2 = i_1 + 1;
		beta = (time - fTimeSteps[i_1])/
		       (fTimeSteps[i_2] - fTimeSteps[i_1]);	
	}

	/* open results file */	
	ExodusT exo(cout);
	if (!exo.OpenRead(fExternalFieldFile))
	{
		cout << "\n UpLagr_ExternalFieldT::InitStep: error opening file: " 
		     << fExternalFieldFile << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* loop over variables */
	for (int i = 0; i < fExternalField.MinorDim(); i++)
	{
		/* assemble values from time 1 */
		if (1.0 - beta > kSmall)
		{
			/* read data */
			exo.ReadNodalVariable(i_1 + 1, i + 1, fNodalValues);

			/* assemble */
			AssembleField(i, 1.0 - beta, fNodalValues);
		}
		
		/* assemble values from time 2 */
		if (beta > kSmall)
		{
			/* read data */
			exo.ReadNodalVariable(i_2 + 1, i + 1, fNodalValues);

			/* assemble */
			AssembleField(i, beta, fNodalValues);
		}
	}
}

/***********************************************************************
* Protected
***********************************************************************/

bool UpLagr_ExternalFieldT::NextElement(void)
{
	/* inherited */
	bool result = UpdatedLagrangianT::NextElement();

	/* gather element values of the external field */
	if (result) SetLocalU(fLocExternalField);
	
	return result;
}

/***********************************************************************
* Private
***********************************************************************/

/* assemble values into external field using the node map */
void UpLagr_ExternalFieldT::AssembleField(int col, double scale, 
	const dArrayT& values)
{
	/* checks */
	if (col < 0 || col >= fExternalField.MinorDim()) throw ExceptionT::kOutOfRange;
	if (values.Length() != fNodeMap.Length()) throw ExceptionT::kSizeMismatch;

	for (int i = 0; i < fNodeMap.Length(); i++)
	{
		int node = fNodeMap[i];
		fExternalField(node, col) += scale*values[i];
	}
}
