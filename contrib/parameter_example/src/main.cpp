/* $Id: main.cpp,v 1.1.2.2 2003-05-03 09:08:27 paklein Exp $ */
#include "Environment.h"
#include "house.h"
#include "ofstreamT.h"

/* parameter headers */
#include "ParameterListT.h"
#include "ParameterTreeT.h"
#include "XML_Attribute_FormatterT.h"
#include "DotLine_FormatterT.h"
#include "expat_ParseT.h"

/* application-specific headers */
#include "house.h"

/* prototypes */
void define_parameters(ParameterInterfaceT& root, const char* out_path);
void validate_list(ParameterInterfaceT& root, ParameterListT& list, 
	const char* file_path);

int main(int argc, char** argv)
{
#pragma unused(argc)
#pragma unused(argv)

//	StringT cwd = "/Volumes/Uster/USERS/paklein/Code/protected-tahoe/BRANCH_parameter_1/contrib/parameter_example/example/";
	StringT cwd = "./";

	/* dump a parameter description */
	StringT dtd_path;
	dtd_path.Append(cwd, "house.dtd");
	house house1;
	define_parameters(house1, dtd_path);

	/* read and validate input from source */
	StringT source_path;
	source_path.Append(cwd, "house2.xml");
	ParameterListT list1;
	house house2;
	validate_list(house2, list1, source_path);
	
	/* dump validated list */
	StringT pp_path;
	pp_path.Append(cwd, "house.out");
	ofstreamT pp_out(pp_path);
	DotLine_FormatterT pp_format;
	pp_format.SetTabWidth(4);
	pp_format.InitParameterFile(pp_out);
	pp_format.WriteParameterList(pp_out, list1);
	pp_format.CloseParameterFile(pp_out);

	/* write validated XML */
	StringT valid_path;
	valid_path.Append(cwd, "house_valid.xml");
	ofstreamT valid_out(valid_path);
	XML_Attribute_FormatterT att_format;
	att_format.InitParameterFile(valid_out);
	att_format.WriteParameterList(valid_out, list1);
	att_format.CloseParameterFile(valid_out);
	
	/* build house */
	house house3;
	house3.SetParameters(list1);

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

void validate_list(ParameterInterfaceT& root, ParameterListT& valid_list, const char* file_path)
{
	/* read values */
	expat_ParseT parser;
	ParameterListT raw_list;
	raw_list.SetDuplicateListNames(true);
	parser.Parse(file_path, raw_list);

	/* generate validated parameter list */
	ParameterTreeT tree;
	ParameterListT* input_list = raw_list.List(root.Name().Pointer());
	if (!input_list) ExceptionT::GeneralFail("validate_list", "list \"%s\" not found", root.Name().Pointer());
	tree.Validate(root, *input_list, valid_list);
}
