/* $Id: FEManagerT.ParseInput.cpp,v 1.2.2.2 2004-08-09 20:56:17 paklein Exp $ */
#include "FEManagerT.h"

#include "ofstreamT.h"
#include "ParameterTreeT.h"
#include "expat_ParseT.h"
#include "XML_Attribute_FormatterT.h"
#include "CommunicatorT.h"

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
		TaskT task = kParameters;
		
		FEManagerT* fe = NULL;
		try {

			/* construct FEManagerT */
			fe = FEManagerT::New(raw_list.Name(), path, output, comm, argv, task);
			if (!fe) 
				ExceptionT::GeneralFail(caller, "failed to construct \"%s\"",
					raw_list.Name().Pointer());

			/* validate */
			ParameterTreeT tree;
			tree.Validate(*fe, raw_list, params);
			
			/* clean up */
			delete fe;
		}
		
		catch (ExceptionT::CodeT error) {
			delete fe;
			ExceptionT::Throw(error, caller, "validation failed");
		}
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
