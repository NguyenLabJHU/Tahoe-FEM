/* $Id: iConsoleBaseT.cpp,v 1.1.1.1 2001-01-25 20:56:28 paklein Exp $ */
/* created: paklein (12/21/2000)                                          */
/* iConsoleBaseT.cpp                                                      */

#include "iConsoleBaseT.h"
#include <strstream.h>
#include <iomanip.h>
#include <ctype.h>

/* array behavior */
const bool ArrayT<iConsoleBaseT::VariableType>::fByteCopy = true;

/* constructor */
iConsoleBaseT::iConsoleBaseT(void):
	fCommands(0, false),
	fVariables(0, false),
	fVariableTypes(0, true),
	fVariableValues(0, true),
	fVariableIsConst(0, true)
{

}

/* write scope variables */
void iConsoleBaseT::iWriteVariables(ostream& out) const
{
	if (fVariables.Length() == 0)
		out << setw(4) << " " << "<none>" << endl;
	else
	{
		const char* type_names[] = {"integer", "double", "string", "boolean"};
		for (int i = 0; i < fVariables.Length(); i++)
		{
			out << setw(4) << " " << fVariables[i] << " = ";
			WriteVariable(out, i);
			
			/* type */
			out << " (";
			if (fVariableIsConst[i]) out << "const ";
			out << type_names[fVariableTypes[i]] << ")\n";
		}
	}
}

/* execute given command - returns false on fail */
bool iConsoleBaseT::iDoCommand(const StringT& command, StringT& line)
{
#pragma unused(line)
	cout << "unrecognized command: \"" << command << "\"" << endl;
	return false;
}

/* operate on given variable */
bool iConsoleBaseT::iDoVariable(const StringT& variable, StringT& line)
{
	/* resolve variable */
	int dex = -1;
	for (int i = 0; i < fVariables.Length() && dex == -1; i++)
		if (fVariables[i] == variable)
			dex = i;
	if (dex == -1)
		return false;
	else if (fVariableIsConst[dex])
	{
		cout << '\"' << variable << "\" cannot be modified" << endl;
		return false;
	}
	else
	{
		/* resolve operator */
		VariableOperator op = ResolveOperator(line);
		bool OK = true;
		if (op == kFail)
		{
			WriteVariable(cout, dex);
			cout << '\n';
		}
		else if (fVariableTypes[dex] == bool_)
			OK = Operate(*((bool*) fVariableValues[dex]), op, line);
		else if (fVariableTypes[dex] == int_)
			OK = Operate(*((int*) fVariableValues[dex]), op, line);
		else if (fVariableTypes[dex] == double_)
			OK = Operate(*((double*) fVariableValues[dex]), op, line);
		else if (fVariableTypes[dex] == string_)
			OK = Operate(*((StringT*) fVariableValues[dex]), op, line);
		else
		{
			cout << "iConsoleBaseT::DoVariable: unsupported type: "
			     << fVariableTypes[dex] << endl;
			return kFail;
		}
		
		/* not successful */
		if (!OK)
			cout << "could not operate on \"" << fVariables[dex] << "\" with \""
			     << line << '\"' << endl;
		return OK;
	}
}

/************************************************************************
* Protected
************************************************************************/

/* add command to the dictionary - true if added */
bool iConsoleBaseT::iAddCommand(const StringT& command)
{
	if (fCommands.AppendUnique(command))
	{
		Sort(fCommands);
		return true;
	}
	else
	{
		cout << " iConsoleBaseT::iAddCommand: duplicate command not added: "
		     << command << endl;
		return false;
	}
}

/* add variable */
bool iConsoleBaseT::iAddVariable(const StringT& name, bool& variable)
{
	return AddVariable(name, bool_, (void*) &variable, false);
}

bool iConsoleBaseT::iAddVariable(const StringT& name, const bool& variable)
{
	return AddVariable(name, bool_, (void*) &variable, true);
}

bool iConsoleBaseT::iAddVariable(const StringT& name, int& variable)
{
	return AddVariable(name, int_, (void*) &variable, false);
}

bool iConsoleBaseT::iAddVariable(const StringT& name, const int& variable)
{
	return AddVariable(name, int_, (void*) &variable, true);
}

bool iConsoleBaseT::iAddVariable(const StringT& name, double& variable)
{
	return AddVariable(name, double_, (void*) &variable, false);
}

bool iConsoleBaseT::iAddVariable(const StringT& name, const double& variable)
{
	return AddVariable(name, double_, (void*) &variable, true);
}

bool iConsoleBaseT::iAddVariable(const StringT& name, StringT& variable)
{
	return AddVariable(name, string_, (void*) &variable, false);
}

bool iConsoleBaseT::iAddVariable(const StringT& name, const StringT& variable)
{
	return AddVariable(name, string_, (void*) &variable, true);
}

/* alphabetize the list */
void iConsoleBaseT::Sort(ArrayT<StringT>& list) const
{
	if (list.Length() < 2) return;
	bool changed;
	do {
		changed = false;
		for (int i = 1; i < list.Length(); i++)
			if (strcmp(list[i], list[i-1]) < 0)
			{
				list[i].Swap(list[i-1]);
				changed = true;
			}
	} while (changed == true);
}

/* resolving arguments in () */
bool iConsoleBaseT::ResolveArgument(StringT& source, bool& arg, bool* default_arg)
{
	char* str = source;

	/* find start and end */
	int start, end;
	if (!Position(str, '(', start) || !Position(str + start, ')', end))
		return false;
	else
	{
		/* empty */
		if (end == start + 1)
		{
			if (default_arg == NULL)
				return false;
			else
			{	
				arg = *default_arg;
				source.Drop(end + 1);
				return true;
			}
		}
		else
		{
			/* advance passed white space */
			str += start + 1;
			while (*str != ')' && isspace(*str))
				str++;
			
			/* fail */
			if (*str == ')') return false;
			
			/* resolve */
			switch (*str)
			{
				case '1':
				case 't':
				case 'T':
					arg = true;
					break;
				case '0':
				case 'f':
				case 'F':
					arg = false;
					break;
				default:
					return false;
			}
			
			/* pull */
			source.Drop(end + 1);
			return true;
		}
	}
}

bool iConsoleBaseT::ResolveArgument(StringT& source, int& arg, int* default_arg)
{
	char* str = source;

	/* find start and end */
	int start, end;
	if (!Position(str, '(', start) || !Position(str + start, ')', end))
		return false;
	else
	{
		/* empty */
		if (end == start + 1)
		{
			if (default_arg == NULL)
				return false;
			else
			{	
				arg = *default_arg;
				source.Drop(end + 1);
				return true;
			}
		}
		else
		{
			/* construct string stream */
			istrstream in(str + start + 1);
			int test_val = -99199199;
			arg = test_val;
			in >> arg;
			if (arg == test_val)
			{
				arg = 0;
				return false;
			}
			else
			{
				source.Drop(end + 1);
				return true;
			}
		}
	}
}

bool iConsoleBaseT::ResolveArgument(StringT& source, double& arg, double* default_arg)
{
	char* str = source;

	/* find start and end */
	int start, end;
	if (!Position(str, '(', start) || !Position(str + start, ')', end))
		return false;
	else
	{
		/* empty */
		if (end == start + 1)
		{
			if (default_arg == NULL)
				return false;
			else
			{	
				arg = *default_arg;
				source.Drop(end + 1);
				return true;
			}
		}
		else
		{
			/* construct string stream */
			istrstream in(str + start + 1);
			double test_val = -99199199;
			arg = test_val;
			in >> arg;
			if (arg == test_val)
			{
				arg = 0;
				return false;
			}
			else
			{
				source.Drop(end + 1);
				return true;
			}
		}
	}
}

bool iConsoleBaseT::ResolveArgument(StringT& source, StringT& arg, StringT* default_arg)
{
	char* str = source;

	/* find start and end */
	int start, end;
	if (!Position(str, '(', start) || !Position(str + start, ')', end))
		return false;
	else
	{
		/* empty */
		if (end == 1)
		{
			if (default_arg == NULL)
				return false;
			else
			{	
				arg = *default_arg;
				source.Drop(start + end + 1);
				return true;
			}
		}
		else
		{
			/* string start and end */
			int q_start, q_end;
			if (!Position(str, '"', q_start) || !Position(str + q_start + 1, '"', q_end))
				return false;
			else
			{
				/* copy string */
				arg.Take(source, q_start + 1, q_start + q_end);
				
				/* pull */
				source.Drop(start + end + 1);
				return true;
			}		
		}
	}
}

/* write list of strings with tab and wrap */
void iConsoleBaseT::WriteList(ostream& out, const ArrayT<StringT>& list,
	int tab, int wrap) const
{
	/* checks */
	if (tab < 0 || wrap < 0) throw eGeneralFail;
	if (tab > 0) out << setw(tab-1) << " ";
	int count = tab;
	for (int i = 0; i < list.Length(); i++)
	{
		int length = strlen(list[i]);
		if (count + length + 1 > wrap)
		{
			out << '\n';
			if (tab > 0) out << setw(tab-1) << " ";
			count = tab;
		}
		out << " " << list[i];
		count += length + 1;
	}
}

/* write single variable */
void iConsoleBaseT::WriteVariable(ostream& out, int i) const
{
	switch (fVariableTypes[i])
	{
		case bool_:
			out << ((*((bool*) fVariableValues[i]) == true) ? "true" : "false");
			break;
		case int_:
			out << *((int*) fVariableValues[i]);
			break;
		case double_:
			out << *((double*) fVariableValues[i]);
			break;
		case string_:
			out << "\"" << *((StringT*) fVariableValues[i]) << "\"";
			break;
		default:
			cout << "unrecognized variable type: "
			    << fVariableTypes[i] << endl;
			throw eGeneralFail;
	}
}

iConsoleBaseT::VariableOperator iConsoleBaseT::ResolveOperator(StringT& line) const
{
	/* shift */
	line.DropLeadingSpace();

	/* resolve */
	if (line[0] == '=')
	{
		line.Drop(1);
		return kEQ;
	}
	else if (line[1] == '=')
	{
		VariableOperator op = kFail;
		switch (line[0])
		{
			case '+':
				op = kPlusEQ;
				break;
			case '-':
				op = kMinusEQ;
				break;
			case '*':
				op = kTimesEQ;
				break;
			case '/':
				op = kDivEQ;
				break;		
		}
		if (op != kFail) line.Drop(2);
		return op;
	}
	else
		return kFail;
}

/************************************************************************
* Private
************************************************************************/

/* find first position - returns false if not found */
bool iConsoleBaseT::Position(char* str, char a, int& position)
{
	position = 0;
	int length = strlen(str);
	while (position < length && *str != a)
	{
		str++;
		position++;
	}
	
	if (position == length)
	{
		position = 0;
		return false;
	}
	else
		return true;
}

bool iConsoleBaseT::AddVariable(const StringT& name, VariableType type,
	void* variable, bool is_const)
{
	if (fVariables.AppendUnique(name))
	{
		fVariableTypes.Append(type);
		fVariableValues.Append((void*) variable);
		fVariableIsConst.Append(is_const);
		return true;
	}
	else
	{
		cout << " iConsoleBaseT::AddVariable: duplicate command not added: "
		     << name << endl;
		return false;
	}
}

/* variable operators - return false on fail */
bool iConsoleBaseT::Operate(bool& variable, VariableOperator op, StringT& line) const
{
	if (op != kEQ)
		return kFail;
	else
	{
		int count;
		StringT rhs;
		rhs.FirstWord(line, count, true);
		if (strlen(rhs) == 0) return false;

		/* resolve */
		switch (rhs[0])
		{
			case '1':
			case 't':
			case 'T':
				variable = true;
				break;
			case '0':
			case 'f':
			case 'F':
				variable = false;
				break;
			default:
				return false;
		}
		
		/* successful */
		line.Drop(count);
		return true;
	}
}

bool iConsoleBaseT::Operate(int& variable, VariableOperator op, StringT& line) const
{
	/* convert first word to integer */
	int count;
	StringT rhs_str;
	rhs_str.FirstWord(line, count, false);
	
	istrstream rhs_in(rhs_str.Pointer());
	int test_val = -99199199;
	int rhs = test_val;
	rhs_in >> rhs;
	if (rhs == test_val)
	{
		cout << "could not resolve integer from: \"" << line << "\"" << endl;
		return false;
	}
	else
	{
		line.Drop(count);
		switch (op)
		{
			case kEQ:
				variable = rhs;
				break;
			case kPlusEQ:
				variable += rhs;
				break;
			case kMinusEQ:
				variable -= rhs;
				break;
			case kTimesEQ:
				variable *= rhs;
				break;
			case kDivEQ:
				variable /= rhs;
				break;
			default:				
				cout << "iConsoleBaseT::Operate: unsupported operator: "
				     << op << endl;
				return false;
		}
		return true;
	}
}

bool iConsoleBaseT::Operate(double& variable, VariableOperator op, StringT& line) const
{
	/* convert first word to integer */
	int count;
	StringT rhs_str;
	rhs_str.FirstWord(line, count, false);
	
	istrstream rhs_in(rhs_str.Pointer());
	double test_val = -99199199;
	double rhs = test_val;
	rhs_in >> rhs;
	if (rhs == test_val)
	{
		cout << "could not resolve integer from: \"" << line << "\"" << endl;
		return false;
	}
	else
	{
		line.Drop(count);
		switch (op)
		{
			case kEQ:
				variable = rhs;
				break;
			case kPlusEQ:
				variable += rhs;
				break;
			case kMinusEQ:
				variable -= rhs;
				break;
			case kTimesEQ:
				variable *= rhs;
				break;
			case kDivEQ:
				variable /= rhs;
				break;
			default:				
				cout << "iConsoleBaseT::Operate: unsupported operator: "
				     << op << endl;
				return false;
		}
		return true;
	}
}

bool iConsoleBaseT::Operate(StringT& variable, VariableOperator op, StringT& line) const
{
	if (op != kEQ || op != kPlusEQ)
		return kFail;
	else
	{
		int count;
		StringT rhs;
		rhs.FirstWord(line, count, true);
		line.Drop(count);
		if (op == kEQ)	
			variable = rhs;
		else
			variable.Append(rhs);
		return true;
	}
}
