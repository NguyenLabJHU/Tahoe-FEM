/* $Id: FEManagerT.ParseInput.cpp,v 1.1 2004-07-22 08:16:51 paklein Exp $ */
/* created: paklein (05/22/1996) */
#include "FEManagerT.h"

#include "ofstreamT.h"
#include "ParameterTreeT.h"
#include "expat_ParseT.h"
#include "XML_Attribute_FormatterT.h"
#include "CommunicatorT.h"

/* element configuration header */
#include "ElementsConfig.h"

#ifdef BRIDGING_ELEMENT
#include "MultiManagerT.h"
#endif

using namespace Tahoe;

/* parse input file and valid */
void FEManagerT::ParseInput(const StringT& path, ParameterListT& params, bool validate,
	bool echo_input, bool echo_valid, const ArrayT<StringT>& argv)
{
	const char caller[] = "FEManagerT::ParseInput";

	/* construct parser */
	expat_ParseT parser;

	/* read values */
	ParameterListT tmp_list;
	ParameterListT& raw_list = (validate) ? tmp_list : params;
	raw_list.SetDuplicateListNames(true);
	parser.Parse(path, raw_list);

	/* echo to XML */
	if (echo_input)  {	
		StringT echo_path;
		echo_path.Root(path);
		echo_path.Append(".echo.xml");
		ofstreamT echo_out(echo_path);
		XML_Attribute_FormatterT att_format(XML_Attribute_FormatterT::DTD);
		att_format.InitParameterFile(echo_out);
		att_format.WriteParameterList(echo_out, raw_list);
		att_format.CloseParameterFile(echo_out);
	}

	/* build validated parameter list */
	if (validate)
	{
		/* parameters currently needed to construct an FEManagerT */
		ofstreamT output;
		CommunicatorT comm;

		ParameterTreeT tree;
		if (raw_list.Name() == "tahoe")
		{
			FEManagerT fe_man(path, output, comm, argv);
			tree.Validate(fe_man, raw_list, params);	
		}
#ifdef BRIDGING_ELEMENT
		else if (raw_list.Name() == "tahoe_multi")
		{
			MultiManagerT multi_man(path, output, comm, argv);
			tree.Validate(multi_man, raw_list, params);	
		}
#endif
		else
			ExceptionT::GeneralFail(caller, "unrecorngized list \"%s\"",
				raw_list.Name().Pointer());
	}

	/* write validated XML */
	if (echo_valid) {
		StringT valid_path;
		valid_path.Root(path);
		valid_path.Append(".valid.xml");
		ofstreamT valid_out(valid_path);
		XML_Attribute_FormatterT att_format(XML_Attribute_FormatterT::DTD);
		att_format.InitParameterFile(valid_out);
		att_format.WriteParameterList(valid_out, params);
		att_format.CloseParameterFile(valid_out);
	}
}
