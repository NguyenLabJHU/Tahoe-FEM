/* $Id: main.cpp,v 1.1.2.1 2003-04-28 08:17:08 paklein Exp $ */
#include "house.h"

/* parameter headers */
#include "ParameterListT.h"
#include "ParameterTreeT.h"
#include "XML_Attribute_FormatterT.h"

/* application-specific headers */
#include "house.h"

/* prototypes */
void define_parameters(ParameterInterfaceT& root, const char* out_path);

int main(int argc, char** argv)
{
#pragma unused(argc)
#pragma unused(argv)

	/* dump a parameter description */
	house house1;
	define_parameters(house1, "/Volumes/Uster/USERS/paklein/Code/protected-tahoe/BRANCH_parameter_1/contrib/parameter_example/example/house.dtd");

	return 0;
}

void define_parameters(ParameterInterfaceT& root, const char* out_path)
{
	ParameterTreeT tree;
	tree.BuildDescription(root);

	/* write description */
	ofstreamT out;
	out.open(out_path);
	XML_Attribute_FormatterT attribute;
	attribute.InitDescriptionFile(out);
	
	const ArrayT<ParameterListT*>& branches = tree.Branches();
	for (int i = 0; i < branches.Length(); i++)
		attribute.WriteDescription(out, *(branches[i]));

	attribute.CloseDescriptionFile(out);
	out.close();
}
