/* $Id: XML_Attribute_FormatterT.cpp,v 1.10 2004-03-06 17:25:55 paklein Exp $ */
#include "XML_Attribute_FormatterT.h"
#include "ParameterListT.h"
#include "ParameterT.h"
#include "BinaryTreeT.h"
#include "LinkedListT.h"

#include <iomanip.h>

using namespace Tahoe;

/* attribute type names */
static const char* a_dtype_names[4] = {"integer", "float", "string", "enumeration"};
static const char* a_dtype(ValueT::TypeT t)
{
	switch (t)
	{
		case ValueT::Integer:
			return a_dtype_names[0];
			break;
		case ValueT::Double:
			return a_dtype_names[1];
			break;
		case ValueT::Enumeration:
			return a_dtype_names[3];
			break;
		default:
			return a_dtype_names[2];
	}
};

static const char* type_names[4] = {"Integer", "Double", "String", "Enumeration"};
static const char* TypeName(ValueT::TypeT t)
{
	switch (t)
	{
		case ValueT::Integer:
			return type_names[0];
			break;
		case ValueT::Double:
			return type_names[1];
			break;
		case ValueT::String:
			return type_names[2];
			break;
	
		default:
			return NULL;
	}
};

static const char* occur_char_list[4] = {"", "?", "+", "*"};
static const char* occur_char(ParameterListT::OccurrenceT o)
{
	switch (o)
	{
		case ParameterListT::ZeroOrOnce:
			return occur_char_list[1];
			break;
		case ParameterListT::OnePlus:
			return occur_char_list[2];
			break;
		case ParameterListT::Any:
			return occur_char_list[3];
			break;		
		default:
			return occur_char_list[0];
	}
};

static const char* order_char_list[4] = {",", "|"};
static const char* order_char(ParameterListT::ListOrderT o)
{
	if (o == ParameterListT::Choice)
		return order_char_list[1];
	else
		return order_char_list[0]; /* sequence by default */
};

static const char* order_indicator_list[] = {"xs:sequence", "xs:choice"};
static const char* order_indicator(ParameterListT::ListOrderT o)
{
	if (o == ParameterListT::Choice)
		return order_indicator_list[1];
	else
		return order_indicator_list[0]; /* sequence by default */
};

static const char* data_type_names[4] = {"xs:integer", "xs:decimal", "xs:string", "xs:boolean"};
static const char* DataTypeName(ValueT::TypeT t)
{
	switch (t)
	{
		case ValueT::Integer:
			return data_type_names[0];
			break;

		case ValueT::Double:
			return data_type_names[1];
			break;

		case ValueT::Boolean:
			return data_type_names[3];
			break;
		
		default: /* string by default */
			return data_type_names[2];
			break;
	}
};

static const char* data_restrictions[5] = {"xs:enumeration", "xs:maxExclusive", "xs:maxInclusive", "xs:minExclusive", "xs:minInclusive"};
static const char* DataRestriction(LimitT::BoundT bound)
{
	switch (bound)
	{
		case LimitT::Only:
			return data_restrictions[0];
			break;

		case LimitT::Upper:
			return data_restrictions[1];
			break;

		case LimitT::UpperInclusive:
			return data_restrictions[2];
			break;

		case LimitT::Lower:
			return data_restrictions[3];
			break;

		case LimitT::LowerInclusive:
			return data_restrictions[4];
			break;

		default:
			return NULL;
	}
};

static const char* data_indicators[4] = {"", "minOccurs='0'", "maxOccurs='unbounded'", "minOccurs='0' maxOccurs='unbounded'"};
static const char* DataIndicator(ParameterListT::OccurrenceT occur)
{
	switch (occur)
	{
		case ParameterListT::ZeroOrOnce:
			return data_indicators[1];
			break;

		case ParameterListT::OnePlus:
			return data_indicators[2];
			break;

		case ParameterListT::Any:
			return data_indicators[3];
			break;

		default:
			return data_indicators[0];
	}
};

/* macros */
static int Max(int a, int b) { return (a > b) ? a : b; };

XML_Attribute_FormatterT::XML_Attribute_FormatterT(DocTypeT doc_type):
	fDocType(doc_type),
	fDocumentRoot("parameter_list")
{

}

void XML_Attribute_FormatterT::SetDocDescription(const StringT& doc_root, const StringT& description_path)
{
	fDocumentRoot = doc_root;
	fPath = description_path;
}

bool XML_Attribute_FormatterT::InitParameterFile(ostream& out) const
{
	/* write XML header */
	out << "<?xml version='1.0' encoding='US-ASCII' standalone='no' ?>" << '\n';

	/* DTD is set */
	if (fPath.StringLength() > 0)
		out << "<!DOCTYPE " << fDocumentRoot << " SYSTEM '" << fPath << "'>" << '\n';
	
	out.flush();
	return true;
}

bool XML_Attribute_FormatterT::WriteParameterList(ostream& out, const ParameterListT& list) const
{
	/* get max attribute name length */
	int name_width = ParameterWidth(list);
	if (list.Description().StringLength() > 0) 
		name_width = Max(name_width, strlen("description"));
	name_width += 1; /* extra space */

	/* list tag name */
	const StringT& tag = list.Name();

	/* open tag */
	out << '\n' << Tab() << "<" << tag;
	
	/* description */
	if (list.Description().StringLength() > 0)
		out << '\n' << Tab() << setw(name_width) << "description" << "='" << list.Description() << "'";
	
	/* plain parameters as attributes */
	const ArrayT<ParameterT>& params = list.Parameters();
	for (int i = 0; i < params.Length(); i++)
		out << '\n' << Tab() << setw(name_width) << params[i].Name() << "='" << params[i] << "'";

	/* non-empty */
	if (list.Lists().Length() > 0)
	{
		/* close tag */
		out << ">";
		TabOut();
	
		/* nested paramater lists */
		const ArrayT<ParameterListT>& nested_lists = list.Lists();
		for (int i = 0; i < nested_lists.Length(); i++)
			WriteParameterList(out, nested_lists[i]);

		TabIn();
		out << '\n' << Tab() << "</" << tag << ">";
	}
	else /* empty element */
		out << "/>";

	out.flush();
	return true;
}

bool XML_Attribute_FormatterT::CloseParameterFile(ostream& out) const
{
	/* add extra newline and flush */
	out << endl;
	return true;
}

bool XML_Attribute_FormatterT::InitDescriptionFile(ostream& out) const
{
	if (fDocType == DTD)
	{
		out << "<?xml version='1.0' encoding='UTF-8'?>" << '\n';
		out << '\n'; 
		out << "<!--Generated by Tahoe::XML_Attribute_FormatterT $Revision: 1.10 $-->";
		return true;
	}
	else if (fDocType == XSD)
	{	
		out << "<?xml version='1.0' encoding='UTF-8'?>" << '\n';
		out << "<xs:schema xmlns:xs='http://www.w3.org/2001/XMLSchema' >" << '\n';
		out << '\n'; 
		out << "<!--Generated by Tahoe::XML_Attribute_FormatterT $Revision: 1.10 $-->";
		return true;	
	}
	else
		return false;
}

bool XML_Attribute_FormatterT::CloseDescriptionFile(ostream& out) const
{
	if (fDocType == XSD) /* close schema */
		out << "</xs:schema>";

	/* add extra newline and flush */
	out << endl;
	return true;
}

/* write the data description */
bool XML_Attribute_FormatterT::WriteDescription(ostream& out, const ParameterListT& list) const
{
	BinaryTreeT<StringT> tags;
	if (fDocType == DTD)
		return DoWriteDTD(out, list, tags);
	else if (fDocType == XSD)
	{
		out.setf(ios::fixed, ios::floatfield);
		return DoWriteXSD(out, list, tags);
	}
	else
		return false;
}

/*************************************************************************
 * Private
 *************************************************************************/

/* write the data description */
bool XML_Attribute_FormatterT::DoWriteDTD(ostream& out, const ParameterListT& list, BinaryTreeT<StringT>& tags) const
{
	const char caller[] = "XML_Attribute_FormatterT::DoWriteDTD";

	/* not for inlined lists */
	if (!list.Inline())
	{
		/* check tag name */
		if (!tags.InsertUnique(list.Name()))
			return false;

		/* open description */
		out << "\n<!ELEMENT " << list.Name();
	
		/* child  */
		const ArrayT<ParameterListT>& nested_lists = list.Lists();
		const ArrayT<StringT>& references = list.References();
		if (nested_lists.Length() > 0 || references.Length() > 0)
		{
			/* open */
			out << " (";
		
			bool first = true;
		
			/* nested lists */
			const ArrayT<ParameterListT::OccurrenceT>& list_occur = list.ListOccurrences();
			for (int i = 0; i < nested_lists.Length(); i++)
			{
				if (first) first = false;
				else out << order_char(list.ListOrder());

				if (!nested_lists[i].Inline())
					out << "\n\t" << nested_lists[i].Name() << occur_char(list_occur[i]);
				else /* inlined list */
				{
					const ParameterListT& list_group = nested_lists[i];
				
					/* open */
					out << "\n\t(";
					bool first_list = true;
				
					const ArrayT<ParameterListT::OccurrenceT>& group_list_occur = list_group.ListOccurrences();
					const ArrayT<ParameterListT>& group_nested_lists = list_group.Lists();
					for (int j = 0; j < group_nested_lists.Length(); j++) {

						if (first_list) first_list = false;
						else out << order_char(list_group.ListOrder());
					
						out << group_nested_lists[j].Name() << occur_char(group_list_occur[j]);
					}

					const ArrayT<StringT>& group_references = list_group.References();
					const ArrayT<ParameterListT::OccurrenceT>& group_ref_occur = list_group.ReferenceOccurrences();
					for (int j = 0; j < group_references.Length(); j++) {

						if (first_list) first_list = false;
						else out << order_char(list_group.ListOrder());
					
						out << group_references[j] << occur_char(group_ref_occur[j]);
					}

					/* close */
					out << ")" << occur_char(list_occur[i]);
					
					/* check */
					if (first_list)
						ExceptionT::GeneralFail(caller, "group list \"%s\" is empty",
							list_group.Name().Pointer());
				}
			}

			/* references */
			const ArrayT<ParameterListT::OccurrenceT>& ref_occur = list.ReferenceOccurrences();
			for (int i = 0; i < references.Length(); i++)
			{
				if (first) first = false;
				else out << ',';
				out << "\n\t" << references[i] << occur_char(ref_occur[i]);		
			}

			/* close */
			out << ")";
		}
		else /* empty element */
			out << " EMPTY";
	
		/* close description */
		out << ">\n";
	}

	/* write parameters as attributes */
	const ArrayT<ParameterT>& params = list.Parameters();
	if (params.Length() > 0 || list.Description().StringLength() > 0)
	{
		/* get max attribute name length */
		int name_width = ParameterWidth(list);
		if (list.Description().StringLength() > 0) 
			name_width = Max(name_width, strlen("description"));
		name_width += 1; /* extra space */

		/* open attribute list */
		out << "<!ATTLIST " << list.Name();
	
		/* description */
		if (list.Description().StringLength() > 0)
			out << "\n" << setw(name_width) << "description" << " CDATA #FIXED '" << list.Description() << "'";
	
		const ArrayT<ParameterListT::OccurrenceT>& param_occur = list.ParameterOccurrences();
		for (int i = 0; i < params.Length(); i++)
		{
			/* the parameter */
			const ParameterT& parameter = params[i];
		
			out << "\n" << setw(name_width) << parameter.Name();

			/* parameters limits */
			const AutoArrayT<LimitT>& limits = parameter.Limits();
			
			/* check for enumeration */
			bool all_limits_only = (limits.Length() > 1);
			for (int j = 0; all_limits_only && j < limits.Length(); j++)
				all_limits_only = limits[j].Bound() == LimitT::Only;

			/* attribute is an enumeration */
			if (all_limits_only && parameter.Type() != ValueT::Boolean)
			{
				out << " (";
				for (int j = 0; j < limits.Length(); j++)
				{
					if (j > 0) out << "|";
					out << limits[j];
				}
				out << ")";
			}
			else if (parameter.Type() == ValueT::Boolean) /* boolean */
				out << " (true|false)";			
			else /* all other attributes */
				out << " CDATA";

			/* default value */
			const ValueT* def = parameter.Default();
			if (def) 
				out << " '" << (*def) << "'";
			else
			{
				switch (param_occur[i])
				{
					case ParameterListT::ZeroOrOnce:
						out << " #IMPLIED";
						break;

					case ParameterListT::Once:
						out << " #REQUIRED";
						break;

					case ParameterListT::OnePlus:
						cout << "\n XML_Attribute_FormatterT: occurrence of " << parameter.Name() << " changed from OnePlus to Once" << endl;
						out << " #REQUIRED";
						break;

					case ParameterListT::Any:
						cout << "\n XML_Attribute_FormatterT: occurrence of " << parameter.Name() << " changed from Any to ZeroOrOnce" << endl;
						out << " #IMPLIED";
						break;		
						
					default:		
						cout << "\n XML_Attribute_FormatterT: occurrence of " << parameter.Name() << " changed to ZeroOrOnce" << endl;
						out << " #IMPLIED";
						break;		
				}
			}
		}
		
		/* attribute data types */
		if (params.Length() > 0)
		{
			bool first = true;
			for (int i = 0; i < params.Length(); i++)
			{
				if (params[i].Type() == ParameterT::Double || params[i].Type() == ParameterT::Integer)
				{
					/* open line */
					if (first) {
						out << "\n" << setw(name_width) << "a-dtype" << " NMTOKENS '";
						first = false;
					}
					else /* space between entries */
						out << " ";

					out << params[i].Name() << " " << a_dtype(params[i].Type());
				}
			}
			
			/* close line */
			if (!first) out << "'";
		}

		/* close attributes */
		out << ">\n";
	}

	/* recursively process nested_lists */
	bool nested_lists_OK = true;
	const ArrayT<ParameterListT>& nested_lists = list.Lists();
	for (int i = 0; i < nested_lists.Length(); i++)
		nested_lists_OK = DoWriteDTD(out, nested_lists[i], tags) && nested_lists_OK;

	/* return */
	return nested_lists_OK;
}

/* write the XML schema data description */
bool XML_Attribute_FormatterT::DoWriteXSD(ostream& out, const ParameterListT& list, BinaryTreeT<StringT>& tags) const
{
//	const char caller[] = "XML_Attribute_FormatterT::DoWriteXSD";

	/* not for inlined lists */
	if (!list.Inline())
	{
		/* check tag name */
		if (!tags.InsertUnique(list.Name()))
			return false;

		/* open description */
		out << "\n<xs:element name='" << list.Name() << "'>\n";
	
		/* nested elements and attributes */
		const ArrayT<ParameterT>& params = list.Parameters();
		const ArrayT<ParameterListT>& nested_lists = list.Lists();
		if (nested_lists.Length() > 0 || params.Length() > 0 || list.Description().StringLength() > 0)
		{
			/* open complex type definition */
			out << "<xs:complexType>\n";

			/* nested lists */
			if (nested_lists.Length() > 0)
			{
				/* open order indicator */
				out << "\t\t<" << order_indicator(list.ListOrder()) << ">\n";
						
				const ArrayT<ParameterListT::OccurrenceT>& list_occur = list.ListOccurrences();
				for (int i = 0; i < nested_lists.Length(); i++)
				{
					/* simple element */
					if (!nested_lists[i].Inline())
						out << "\t\t\t<xs:element ref='" 
						    << nested_lists[i].Name() << "' " 
						    << DataIndicator(list_occur[i]) << "/>\n";

					else /* inlined list */
					{
						const ParameterListT& list_group = nested_lists[i];
						const ArrayT<ParameterListT>& group_nested_lists = list_group.Lists();
						const ArrayT<ParameterListT::OccurrenceT>& group_list_occur = list_group.ListOccurrences();

						/* open order indicator */
						out << "\t\t\t<" << order_indicator(list_group.ListOrder()) << " " << DataIndicator(list_occur[i]) << ">\n";

						/* write nested lists */
						for (int j = 0; j < group_nested_lists.Length(); j++)
							out << "\t\t\t\t<xs:element ref='" 
						    	<< group_nested_lists[j].Name() << "' " 
						    	<< DataIndicator(group_list_occur[j]) << "/>\n";

						/* close order indicator */
						out << "\t\t\t</" << order_indicator(list_group.ListOrder()) << ">\n";
					}
				}

				/* close order indicator */
				out << "\t\t</" << order_indicator(list.ListOrder()) << ">\n";
			}

			/* parameters */
			if (params.Length() > 0 || list.Description().StringLength() > 0)
			{	
				/* description */
				if (list.Description().StringLength() > 0)
					out << "\t\t<xs:attribute name='description' fixed='" << list.Description() << "'/>\n";
	
				const ArrayT<ParameterListT::OccurrenceT>& param_occur = list.ParameterOccurrences();
				for (int i = 0; i < params.Length(); i++)
				{
					/* the parameter */
					const ParameterT& parameter = params[i];

					/* open attribute */
					out << "\t\t<xs:attribute name='" << parameter.Name() << "'";

					/* default value */
					const ValueT* def = parameter.Default();
					if (def) 
						out << " default='" << (*def) << "'";
					else
					{
						/* required or optional */				
						switch (param_occur[i])
						{
							case ParameterListT::ZeroOrOnce:
								out << " use='optional'";
								break;

							case ParameterListT::Once:
								out << " use='required'";
								break;

							case ParameterListT::OnePlus:
								cout << "\n XML_Attribute_FormatterT: occurrence of " << parameter.Name() << " changed from OnePlus to Once" << endl;
								out << " use='required'";
								break;

							case ParameterListT::Any:
								cout << "\n XML_Attribute_FormatterT: occurrence of " << parameter.Name() << " changed from Any to ZeroOrOnce" << endl;
								out << " use='optional'";
								break;		
						
							default:		
								cout << "\n XML_Attribute_FormatterT: occurrence of " << parameter.Name() << " changed to ZeroOrOnce" << endl;
								out << " use='optional'";
								break;
						}
					}

					/* parameters limits */
					const AutoArrayT<LimitT>& limits = parameter.Limits();
					if (limits.Length() == 0 || parameter.Type() == ValueT::Boolean) /* no restrictions */
						/* include data type and close description */
						out << " type='" << DataTypeName(parameter.Type()) << "'/>\n";
					else
					{
						/* close tag */
						out << ">\n";
					
						/* open limits */
						out << "\t\t\t<xs:simpleType>\n";	
						out << "\t\t\t<xs:restriction base='" << DataTypeName(parameter.Type()) << "'>\n";					

						/* check for enumeration */
						bool all_limits_only = (limits.Length() > 1);
						for (int j = 0; all_limits_only && j < limits.Length(); j++)
							all_limits_only = limits[j].Bound() == LimitT::Only;

						/* attribute is an enumeration */
						if (all_limits_only)
						{
							for (int j = 0; j < limits.Length(); j++)
								out << "\t\t\t\t<xs:enumeration value='" << limits[j] << "'/>\n";
						}
						else /* other limits */
						{
							for (int j = 0; j < limits.Length(); j++)
								out << "\t\t\t\t<" << DataRestriction(limits[j].Bound()) << " value='" << limits[j] << "'/>\n";
						}

						/* close limits */
						out << "\t\t\t</xs:restriction>\n";	
						out << "\t\t\t</xs:simpleType>\n";	
						
						/* close description */
						out << "\t\t</xs:attribute>\n";
					}
				}
			}

			/* close complex type definition */
			out << "</xs:complexType>\n";
		}
		else /* empty element */
			out << "<xs:complexType/>\n";
	
		/* close description */
		out << "</xs:element>\n";
	}

	/* recursively process nested_lists */
	bool nested_lists_OK = true;
	const ArrayT<ParameterListT>& nested_lists = list.Lists();
	for (int i = 0; i < nested_lists.Length(); i++)
		nested_lists_OK = DoWriteXSD(out, nested_lists[i], tags) && nested_lists_OK;

	/* return */
	return nested_lists_OK;
}

int XML_Attribute_FormatterT::ParameterWidth(const ParameterListT& list) const
{
	/* get max attribute name length */
	const ArrayT<ParameterT>& params = list.Parameters();
	int name_width = 0;
	for (int i = 0; i < params.Length(); i++)
	{
		int len = params[i].Name().StringLength();
		name_width = (len > name_width) ? len : name_width;
	}
	return name_width;
}

int XML_Attribute_FormatterT::ListWidth(const ParameterListT& list) const
{
	/* get max attribute name length */
	const ArrayT<ParameterListT>& lists = list.Lists();
	int name_width = 0;
	for (int i = 0; i < lists.Length(); i++)
	{
		int len = lists[i].Name().StringLength();
		name_width = (len > name_width) ? len : name_width;
	}
	return name_width;
}

int XML_Attribute_FormatterT::ReferenceWidth(const ParameterListT& list) const
{
	/* get max attribute name length */
	const ArrayT<StringT>& refs = list.References();
	int name_width = 0;
	for (int i = 0; i < refs.Length(); i++)
	{
		int len = refs[i].StringLength();
		name_width = (len > name_width) ? len : name_width;
	}
	return name_width;
}
