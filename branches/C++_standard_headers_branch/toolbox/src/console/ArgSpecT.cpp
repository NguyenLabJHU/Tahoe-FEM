/* $Id: ArgSpecT.cpp,v 1.13 2003-11-10 22:14:15 cjkimme Exp $ */
#include "ArgSpecT.h"
#include <ctype.h>

using namespace Tahoe;

/* array copy behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ArgSpecT*>::fByteCopy = true; 
DEFINE_TEMPLATE_STATIC const bool ArrayT<ArgSpecT>::fByteCopy = false;
} /* namespace Tahoe */

/* type names */
static const char* type_names[] = {
	"integer",
	"double",
	"string",
	"boolean",
	"float"};

/* constructor */
ArgSpecT::ArgSpecT(ArgTypeT t):
	fType(t),
	fDefault(NULL),
	fValue(NULL)
{

}

/* constructor */
ArgSpecT::ArgSpecT(ArgTypeT type, const StringT& name):
	fName(name),
	fType(type),
	fDefault(NULL),
	fValue(NULL)
{

}

/* copy constructor */
ArgSpecT::ArgSpecT(const ArgSpecT& arg):
	fDefault(NULL),
	fValue(NULL)
{
	operator=(arg);
}

/* destructor */
ArgSpecT::~ArgSpecT(void) 
{
	ClearDefault();
	ClearValue();
}

/* clear the default value */
void ArgSpecT::ClearDefault(void)
{
	if (fType == int_)
		delete (int*) fDefault;
	else if (fType == float_)
		delete (float*) fDefault;
	else if (fType == double_)
		delete (double*) fDefault;
	else if (fType == string_)
		delete (StringT*) fDefault;
	else if (fType == bool_)
		delete (bool*) fDefault;
	else {
		cout << "\n ArgSpecT::ClearDefault: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
	fDefault = NULL;
}

/* return the type name */
const char* ArgSpecT::TypeName(void) const { return type_names[Type()]; }

/* clear the value */
void ArgSpecT::ClearValue(void)
{
	if (fType == int_)
		delete (int*) fValue;
	else if (fType == float_)
		delete (float*) fValue;
	else if (fType == double_)
		delete (double*) fValue;
	else if (fType == string_)
		delete (StringT*) fValue;
	else if (fType == bool_)
		delete (bool*) fValue;
	else {
		cout << "\n ArgSpecT::ClearValue: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
	fValue = NULL;
}

/* set default value */
void ArgSpecT::SetDefault(const int& a)
{
	if (fType == int_) {
		
		int* tmp = new int;
		*tmp = a;
		fDefault = (void*) tmp;
	}
	else if (fType == float_)
		SetDefault(float(a));
	else if (fType == double_)
		SetDefault(double(a));
	else {
		cout << "\n ArgSpecT::SetDefault: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::SetDefault(const double& a)
{
	if (fType == double_) {
		ClearDefault();
		double* tmp = new double;
		*tmp = a;
		fDefault = (void*) tmp;
	}
	else {
		cout << "\n ArgSpecT::SetDefault: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}	
}

void ArgSpecT::SetDefault(const StringT& a)
{
	if (fType == string_) {
		ClearDefault();
		StringT* tmp = new StringT;
		*tmp = a;
		fDefault = (void*) tmp;
	}
	else {
		cout << "\n ArgSpecT::SetDefault: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::SetDefault(const bool& a)
{
	if (fType == bool_) {
		ClearDefault();
		bool* tmp = new bool;
		*tmp = a;
		fDefault = (void*) tmp;
	}
	else {
		cout << "\n ArgSpecT::SetDefault: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::SetDefault(const float& a)
{
	if (fType == float_) {
		ClearDefault();
		float* tmp = new float;
		*tmp = a;
		fDefault = (void*) tmp;
	}
	else if (fType == double_)
		SetDefault(double(a));
	else {
		cout << "\n ArgSpecT::SetDefault: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/* set default value */
void ArgSpecT::GetDefault(int& a) const
{
	if (fType == int_)
		a = *((int*) fDefault);
	else {
		cout << "\n ArgSpecT::GetDefault: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::GetDefault(double& a) const
{
	if (fType == double_)
		a = *((double*) fDefault);
	else {
		cout << "\n ArgSpecT::GetDefault: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::GetDefault(StringT& a) const
{
	if (fType == string_)
		a = *((StringT*) fDefault);
	else {
		cout << "\n ArgSpecT::GetDefault: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::GetDefault(bool& a) const
{
	if (fType == bool_)
		a = *((bool*) fDefault);
	else {
		cout << "\n ArgSpecT::GetDefault: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::GetDefault(float& a) const
{
	if (fType == float_)
		a = *((float*) fDefault);
	else {
		cout << "\n ArgSpecT::GetDefault: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/* set the value to its default */
void ArgSpecT::SetValueToDefault(void)
{	
	/* set default */
	if (HasDefault())
	{
		switch (fType)
		{
			case int_:
			{
				int a;
				GetDefault(a);
				SetValue(a);
				break;
			}
			case double_:
			{
				double a;
				GetDefault(a);
				SetValue(a);
				break;
			}
			case string_:
			{
				StringT a;
				GetDefault(a);
				SetValue(a);
				break;
			}
			case bool_:
			{
				bool a;
				GetDefault(a);
				SetValue(a);
				break;
			}
			case float_:
			{
				float a;
				GetDefault(a);
				SetValue(a);
				break;
			}
			default:
				cout << "\n ArgSpecT::SetValueToDefault: unknown type: " 
				     << fType << endl;
				throw ExceptionT::kGeneralFail;
		}	
	}	
}

/* assignment operator */
ArgSpecT& ArgSpecT::operator=(const ArgSpecT& rhs)
{
	fName = rhs.Name();
	fType = rhs.Type();
	fPrompt = rhs.Prompt();

	/* set default */
	if (rhs.HasDefault())
	{
		switch (fType)
		{
			case int_:
			{
				int a;
				rhs.GetDefault(a);
				SetDefault(a);
				break;
			}
			case double_:
			{
				double a;
				rhs.GetDefault(a);
				SetDefault(a);
				break;
			}
			case string_:
			{
				StringT a;
				rhs.GetDefault(a);
				SetDefault(a);
				break;
			}
			case bool_:
			{
				bool a;
				rhs.GetDefault(a);
				SetDefault(a);
				break;
			}
			case float_:
			{
				float a;
				rhs.GetDefault(a);
				SetDefault(a);
				break;
			}
			default:
				cout << "\n ArgSpecT::operator=: unknown default type" << endl;
				throw ExceptionT::kGeneralFail;
		}	
	}	
	else
		ClearDefault();

	/* set value */
	if (rhs.HasValue())
	{
		switch (fType)
		{
			case int_:
			{
				int a;
				rhs.GetValue(a);
				SetValue(a);
				break;
			}
			case double_:
			{
				double a;
				rhs.GetValue(a);
				SetValue(a);
				break;
			}
			case string_:
			{
				StringT a;
				rhs.GetValue(a);
				SetValue(a);
				break;
			}
			case bool_:
			{
				bool a;
				rhs.GetValue(a);
				SetValue(a);
				break;
			}
			case float_:
			{
				float a;
				rhs.GetValue(a);
				SetValue(a);
				break;
			}
			default:
				cout << "\n ArgSpecT::operator=: unknown value type" << endl;
				throw ExceptionT::kGeneralFail;
		}	
	}	
	else
		ClearValue();

	return *this;
}

/* read value from stream */
bool ArgSpecT::ReadValue(istream& in)
{
	switch (fType)
	{
		case int_:
		{
			int a = -991991;
			in >> a;
			if (a != -991991) {
				SetValue(a);
				return true;
			}
		}
		case double_:
		{
			double a = -991.991;
			in >> a;
			if (a != -991.991) {
				SetValue(a);
				return true;
			}
		}
		case string_:
		{
			StringT a;
			a.GetLineFromStream(in);
			if (a.StringLength() > 0) {
				SetValue(a);
				return true;
			}
		}
		case bool_:
		{
			/* get next non-white space character */
			char a;	
			in.get(a);
			while (in.good() && isspace(a)) 
				in.get(a);	

			/* resolve */
			switch (a)
			{
				case '1':
				case 't':
				case 'T':
					SetValue(true);
					return true;
				case '0':
				case 'f':
				case 'F':
					SetValue(false);
					return true;
				default:
					return false;
			}
		}
		case float_:
		{
			float a = -91.91;
			in >> a;
			if (a != -91.91) {
				SetValue(a);
				return true;
			}
		}
		default:
			cout << "\n ArgSpecT::operator=: unknown value type" << endl;
			throw ExceptionT::kGeneralFail;
	}

	/* dummy */
	return false;
}

/* set default value */
void ArgSpecT::SetValue(const int& a)
{
	if (fType == int_) {
		ClearValue();
		int* tmp = new int;
		*tmp = a;
		fValue = (void*) tmp;
	}
	else if (fType == float_)
		SetValue(float(a));
	else if (fType == double_)
		SetValue(double(a));
	else {
		cout << "\n ArgSpecT::SetValue: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::SetValue(const double& a)
{
	if (fType == double_) {
		ClearValue();
		double* tmp = new double;
		*tmp = a;
		fValue = (void*) tmp;
	}
	else {
		cout << "\n ArgSpecT::SetValue: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}	
}

void ArgSpecT::SetValue(const StringT& a)
{
	if (fType == string_) {
		ClearValue();
		StringT* tmp = new StringT;
		*tmp = a;
		fValue = (void*) tmp;
	}
	else {
		cout << "\n ArgSpecT::SetValue: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::SetValue(const bool& a)
{
	if (fType == bool_) {
		ClearValue();
		bool* tmp = new bool;
		*tmp = a;
		fValue = (void*) tmp;
	}
	else {
		cout << "\n ArgSpecT::SetValue: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::SetValue(const float& a)
{
	if (fType == float_) {
		ClearValue();
		float* tmp = new float;
		*tmp = a;
		fValue = (void*) tmp;
	}
	else if (fType == double_)
		SetValue(double(a));
	else {
		cout << "\n ArgSpecT::SetValue: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/* set default value */
void ArgSpecT::GetValue(int& a) const
{
	if (fType == int_)
		a = *((int*) fValue);
	else {
		cout << "\n ArgSpecT::GetValue: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::GetValue(double& a) const
{
	if (fType == double_)
		a = *((double*) fValue);
	else {
		cout << "\n ArgSpecT::GetValue: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::GetValue(StringT& a) const
{
	if (fType == string_)
		a = *((StringT*) fValue);
	else {
		cout << "\n ArgSpecT::GetValue: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::GetValue(bool& a) const
{
	if (fType == bool_)
		a = *((bool*) fValue);
	else {
		cout << "\n ArgSpecT::GetValue: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void ArgSpecT::GetValue(float& a) const
{
	if (fType == float_)
		a = *((float*) fValue);
	else {
		cout << "\n ArgSpecT::GetValue: type mismatch" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/* write spec to output */
void ArgSpecT::Write(ostream& out) const
{
	out << '\t';
		
	/* name */
	if (Name().StringLength() > 0) 
		out << Name() << " ";
	else
		out << "(unnamed) ";

	/* type */
	out << "(" << TypeName() << ")";

	/* current value */		
	if (HasValue())
	{
		out << " = ";
		WriteValue(out);
	}
	
	/* default value */
	if (HasDefault())
	{
		out << " (default = ";
		switch (fType)
		{
			case int_:
			{
				int a;
				GetDefault(a);
				out << a;
				break;
			}
			case double_:
			{
				double a;
				GetDefault(a);
				out << a;
				break;
			}
			case string_:
			{
				StringT a;
				GetDefault(a);
				out << '\"' << a << '\"';
				break;
			}
			case bool_:
			{
				bool a;
				GetDefault(a);
				if (a)
					out << "true";
				else
					out << "false";
				break;
			}
			case float_:
			{
				float a;
				GetDefault(a);
				out << a;
				break;
			}
			default:
				cout << "\n ArgSpecT::Write: unknown type: " << fType << endl;
				throw ExceptionT::kGeneralFail;	
		}
		out << ')';
	}
	else
		out << " (no default)";
		
	/* prompt */
	if (Prompt().StringLength() > 0) out << " : " << Prompt();
}

/* write value to the output stream */
void ArgSpecT::WriteValue(ostream& out) const
{
	/* no value */
	if (!HasValue()) 
		out << "<unset>";
	else
	{
		switch (fType)
		{
			case int_:
			{
				int a;
				GetValue(a);
				out << a;
				break;
			}
			case double_:
			{
				double a;
				GetValue(a);
				out << a;
				break;
			}
			case string_:
			{
				StringT a;
				GetValue(a);
				out << '\"' << a << '\"';
				break;
			}
			case bool_:
			{
				bool a;
				GetValue(a);
				if (a)
					out << "true";
				else
					out << "false";
				break;
			}
			case float_:
			{
				float a;
				GetValue(a);
				out << a;
				break;
			}
				default:
				cout << "\n ArgSpecT::WriteValue: unknown type: " << fType << endl;
				throw ExceptionT::kGeneralFail;	
		}
	}
}
