/* $Id: iConsoleBaseT.cpp,v 1.17 2003-11-10 22:14:15 cjkimme Exp $ */
/* created: paklein (12/21/2000) */
#include "iConsoleBaseT.h"
#include "CommandSpecT.h"
#include "ArgSpecT.h"

#ifdef _MSC_VER
#include <strstrea.h>
#else
#include <strstream.h>
#endif
#include <iomanip.h>
#include <ctype.h>

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<iConsoleBaseT::VariableType>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
iConsoleBaseT::iConsoleBaseT(void):
	fCommands(0),
	fVariables(0),
	fVariableTypes(0),
	fVariableValues(0),
	fVariableIsConst(0)
{

}

/* destructor */
iConsoleBaseT::~iConsoleBaseT(void)
{
	/* free command specs */
	for (int i = 0; i < fCommands.Length(); i++)
		delete fCommands[i];
}

/* write scope variables */
void iConsoleBaseT::iWriteVariables(ostream& out) const
{
	if (fVariables.Length() == 0)
		out << setw(4) << " " << "<none>" << endl;
	else
	{
		const char* type_names[] = {"integer", "double", "string", "boolean", "float"};
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
bool iConsoleBaseT::iDoCommand(const CommandSpecT& command, StringT& line)
{
#pragma unused(line)

	cout << "unrecognized command: \"" << command.Name() << "\"" << endl;
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
		else if (fVariableTypes[dex] == float_)
			OK = Operate(*((float*) fVariableValues[dex]), op, line);
		else if (fVariableTypes[dex] == double_)
			OK = Operate(*((double*) fVariableValues[dex]), op, line);
		else if (fVariableTypes[dex] == string_)
			OK = Operate(*((StringT*) fVariableValues[dex]), op, line);
		else
		{
#if __option(extended_errorcheck)
			cout << "iConsoleBaseT::DoVariable: ";
#endif
			cout << "unsupported type: "
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

/* resolve name into function specification */
const CommandSpecT* iConsoleBaseT::iResolveCommand(const StringT& command_name, 
	StringT& line) const
{
	/* find spec */
	CommandSpecT* command_spec = iCommand(command_name);

	/* resolve arguments */
	if (command_spec && ResolveArguments(*command_spec, line, cout, cin))
		return command_spec;
	else 
		return NULL;
}

/* return the command specification with the given name. */
CommandSpecT* iConsoleBaseT::iCommand(const StringT& command_name) const
{
	CommandSpecT* command_spec = NULL;
	for (int i = 0; !command_spec && i < fCommands.Length(); i++)
		if (fCommands[i]->Name() == command_name)
			command_spec = fCommands[i];
	return command_spec;
}

/************************************************************************
* Protected
************************************************************************/

/* write prompt for the specific argument */
void iConsoleBaseT::ValuePrompt(const CommandSpecT& command, int index, 
	ostream& out) const
{
#pragma unused(command)
#pragma unused(index)
#pragma unused(out)
}

/* add command to the dictionary - true if added */
bool iConsoleBaseT::iAddCommand(const CommandSpecT& command)
{
	/* check for duplicate name */
	bool dup = false;
	for (int i = 0; !dup && i < fCommands.Length(); i++)
		dup = command.Name() == fCommands[i]->Name(); 

	if (!dup)
	{
		/* copy command spec */
		CommandSpecT* new_command = new CommandSpecT(command);

		/* add */
		fCommands.Append(new_command);
		
		/* alphabetize */
		SortCommands(fCommands);
		return true;
	}
	else
	{
		cout << " iConsoleBaseT::iAddCommand: duplicate command not added: "
		     << command.Name() << endl;
		return false;
	}
}

/* resolve command arguments */
bool iConsoleBaseT::ResolveArguments(CommandSpecT& command, StringT& line, 
	ostream& out, istream& in) const
{
	/* clear all values */
	command.ClearValues();
	
	/* ordered arguments */
	if (command.Ordered()) {
	
		/* look for arguments in order */
		const ArrayT<ArgSpecT*>& args = command.Arguments();
		bool scan_OK = true;
		for (int i = 0; scan_OK && i < args.Length(); i++)
			scan_OK = ResolveNamedValue(command, i, line, out, in, true);
		
		/* result */
		return scan_OK;
			
	} else { /* arguments in random order */

		/* argument list */
		const ArrayT<ArgSpecT*>& args = command.Arguments();
	
		/* scan line */
		StringT first;
		int count;
		bool scan_OK = true;
//		bool set_defaults = false;
		bool done = false;
		
		first.FirstWord(line, count, true);
		while (!done && scan_OK)
		{
			if (first == ".")
			{
				out << "interrupt." << endl;
				scan_OK = false;
			}
			else if (first.StringLength() == 0)
			{
				/* use defaults for any unset values */
				for (int i = 0; i < args.Length(); i++)
					if (!args[i]->HasValue())
						args[i]->SetValueToDefault();
			
				/* prompt since no argument found */
				bool prompt = true;
				for (int i = 0; prompt && i < args.Length(); i++)
					if (!args[i]->HasValue())
					{
						scan_OK = ResolveNamedValue(command, i, line, out, in, true);
						prompt = false;
					}
					
				/* all filled */
				if (prompt && scan_OK) done = true;
			}
			else
			{
				/* pull word */
				line.Drop(count);
		
				/* find unset argument */
				bool found = false;
				for (int i = 0; !found && scan_OK && i < args.Length(); i++)
				{
					if (args[i]->Name() == first)
					{
						found = true;
						if (args[i]->HasValue())
						{
#if __option(extended_errorcheck)
							out << "iConsoleBaseT::ResolveArguments: ";
#endif
							out << "argument \"" << args[i]->Name() << "\" is already set" << endl;
							scan_OK = false;
						}
						else
							scan_OK = ResolveValue(command, i, line, out, in, true);
					}
				}
			}
			
			/* next */		
			if (scan_OK) first.FirstWord(line, count, true);
		}

		/* result */
		return scan_OK;
	}
}

/* resolve named argument value */
bool iConsoleBaseT::ResolveNamedValue(CommandSpecT& command, int index, StringT& line, 
	ostream& out, istream& in, bool prompt) const
{
	/* the argument */
	ArgSpecT& arg = command.Argument(index);

	/* check name */
	if (arg.Name().StringLength() > 0)
	{
		/* grab first word */
		StringT first;
		int count;
		first.FirstWord(line, count, true);
		
		/* found exit character - minus sign */
		if (first == ".") {
			cout << "interrupt." << endl;
			return false;
		}
	
		/* check name */
		if (first != arg.Name())
		{
			if (!prompt) {
#if __option(extended_errorcheck)
				cout << "iConsoleBaseT::ResolveNamedValue: ";
#endif
				out << "expecting name \"" << arg.Name() << "\", found \"" 
				    << first << '\"' << endl;
				return false;
			}
			else /* prompt */
			{
				/* user-defined prompt */
				const iConsoleBaseT* prompter = command.Prompter();
				if (prompter) prompter->ValuePrompt(command, index, out);
			
				out << "?" << arg.Name() << " ";
				out.flush();
	
				/* read line */
				StringT new_line;
				new_line.GetLineFromStream(in);
				
				/* add name */
				new_line.Prepend(arg.Name(), " ");
	
				/* append the remaining line */
				new_line.Append(" ", line);
				line.Swap(new_line);
		
				/* recurse */
				return ResolveNamedValue(command, index, line, out, in, false);
			}
		}
		else /* shorten line */
			line.Drop(count);
	}

	/* get value */
	return ResolveValue(command, index, line, out, in, prompt);
}

/* resolve named argument value */
bool iConsoleBaseT::ResolveValue(CommandSpecT& command, int index, StringT& line, 
	ostream& out, istream& in, bool prompt) const
{
	/* the argument */
	ArgSpecT& arg = command.Argument(index);

	/* grab first word */
	StringT first;
	int count;
	first.FirstWord(line, count, false);
	
	/* try to resolve argument from first word */
	bool found_arg = false;
	if (first.StringLength() > 0) {
	
		/* found exit character - minus sign */
		if (first == ".") {
			out << "interrupt." << endl;
			return false;
		}

		/* resolve */
#ifdef _MSC_VER
		istrstream s((char *) first);
#else
		istrstream s((const char *) first);
#endif
		try { found_arg = arg.ReadValue(s); }
		catch (ExceptionT::CodeT) { return false; }
		
		/* eat line */
		if (found_arg) line.Drop(count);
		
		/* fail message */
		if (!found_arg && !prompt)
		{
#if __option(extended_errorcheck)
			out << "iConsoleBaseT::ResolveValue: ";
#endif
			out << "failed to resolve " << arg.TypeName() 
			     << " from " << first << endl;
		}
	}
	/* use default value */
	else if (arg.HasDefault()) {
		arg.SetValueToDefault();
		return true;
	}
	
	/* prompt for value */
	if (found_arg)
		return true;
	else if (prompt && !found_arg) {
	
		/* user-defined prompt */
		const iConsoleBaseT* prompter = command.Prompter();
		if (prompter) prompter->ValuePrompt(command, index, out);
	
		/* prompt */
		if (arg.Prompt().StringLength() > 0)
			out << arg.Prompt() << ": ";
		else if (arg.Name().StringLength() > 0)
			out << "?" << arg.Name() << " ";
		else
			out << "(" << arg.TypeName() << "): ";
	
		/* read line */
		StringT new_line;
		new_line.GetLineFromStream(in);
	
		/* append the remaining line */
		new_line.Append(" ", line);
		line.Swap(new_line);
		
		/* recurse */
		return ResolveValue(command, index, line, out, in, false);
	}
	else /* fail */
	{
#if __option(extended_errorcheck)
		out << "iConsoleBaseT::ResolveValue: ";
#endif
		out << "failed to resolve ";
		if (arg.Name().StringLength() > 0) 
			out << '\"' << arg.Name() << "\" ";
		out << '(' << arg.TypeName() << ')' << endl;
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

bool iConsoleBaseT::iAddVariable(const StringT& name, float& variable)
{
	return AddVariable(name, float_, (void*) &variable, false);
}

bool iConsoleBaseT::iAddVariable(const StringT& name, const float& variable)
{
	return AddVariable(name, float_, (void*) &variable, true);
}

bool iConsoleBaseT::iAddVariable(const StringT& name, StringT& variable)
{
	return AddVariable(name, string_, (void*) &variable, false);
}

bool iConsoleBaseT::iAddVariable(const StringT& name, const StringT& variable)
{
	return AddVariable(name, string_, (void*) &variable, true);
}

/* remove named variable */
bool iConsoleBaseT::iDeleteVariable(const StringT& name)
{
	int index = fVariables.PositionOf(name);
	if (index > -1)
	{
		fVariables.DeleteAt(index);
		fVariableTypes.DeleteAt(index);
		fVariableValues.DeleteAt(index);
		fVariableIsConst.DeleteAt(index);
		return true;
	}
	else
		return false;
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

void iConsoleBaseT::SortCommands(ArrayT<CommandSpecT*>& list) const
{
	if (list.Length() < 2) return;
	bool changed;
	do {
		changed = false;
		for (int i = 1; i < list.Length(); i++)
			if (strcmp(list[i]->Name(), list[i-1]->Name()) < 0)
			{
				/* swap */
				CommandSpecT* tmp = list[i];
				list[i] = list[i-1];
				list[i-1] = tmp;

				changed = true;
			}
	} while (changed == true);
}

/* write list of strings with tab and wrap */
void iConsoleBaseT::WriteList(ostream& out, const ArrayT<StringT>& list,
	int tab, int wrap) const
{
	/* checks */
	if (tab < 0 || wrap < 0) throw ExceptionT::kGeneralFail;
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

void iConsoleBaseT::WriteList(ostream& out, const ArrayT<CommandSpecT*>& list,
	int tab, int wrap) const
{
	/* checks */
	if (tab < 0 || wrap < 0) throw ExceptionT::kGeneralFail;
	if (tab > 0) out << setw(tab-1) << " ";
	int count = tab;
	for (int i = 0; i < list.Length(); i++)
	{
		int length = strlen(list[i]->Name());
		if (count + length + 1 > wrap)
		{
			out << '\n';
			if (tab > 0) out << setw(tab-1) << " ";
			count = tab;
		}
		out << " " << list[i]->Name();
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
		case float_:
			out << *((float*) fVariableValues[i]);
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
			throw ExceptionT::kGeneralFail;
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

/* copy variables from the source. \return true if all variables added. */
bool iConsoleBaseT::AddVariables(const iConsoleBaseT& source)
{
	/* add all variables */
	int count = 0;
	for (int i = 0; i < source.fVariables.Length(); i++)
		if (AddVariable(source.fVariables[i], source.fVariableTypes[i], source.fVariableValues[i], source.fVariableIsConst[i]))
			count++;
	return count == source.fVariables.Length();
}

/* remove all variables */
void iConsoleBaseT::DeleteVariables(void)
{
	fVariables.Dimension(0);
	fVariableTypes.Dimension(0);
	fVariableValues.Dimension(0);
	fVariableIsConst.Dimension(0);
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
		return false;
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

bool iConsoleBaseT::Operate(float& variable, VariableOperator op, StringT& line) const
{
	/* convert first word to integer */
	int count;
	StringT rhs_str;
	rhs_str.FirstWord(line, count, false);
	
	istrstream rhs_in(rhs_str.Pointer());
	float test_val = -919;
	float rhs = test_val;
	rhs_in >> rhs;
	if (rhs == test_val)
	{
		cout << "could not resolve float from: \"" << line << "\"" << endl;
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
		cout << "could not resolve double from: \"" << line << "\"" << endl;
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
	if (op != kEQ && op != kPlusEQ)
		return false;
	else
	{
		int count;
		StringT rhs;
		rhs.FirstWord(line, count, false);
		line.Drop(count);
		if (op == kEQ)	
			variable = rhs;
		else
			variable.Append(rhs);
		return true;
	}
}
