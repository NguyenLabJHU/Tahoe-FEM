/* $Id: CommandSpecT.cpp,v 1.12 2003-11-10 22:14:15 cjkimme Exp $ */
#include "CommandSpecT.h"
#include "ArgSpecT.h"

using namespace Tahoe;

/* array copy behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<CommandSpecT*>::fByteCopy = true; 
DEFINE_TEMPLATE_STATIC const bool ArrayT<CommandSpecT>::fByteCopy = false; 
} /* namespace Tahoe */

CommandSpecT::CommandSpecT(const StringT& name, bool ordered_args):
	fName(name),
	fPrompter(NULL),
	fOrdered(ordered_args),
	fArguments(0)
{

}

/* copy constructor */
CommandSpecT::CommandSpecT(const CommandSpecT& command):
	fName(command.Name()),
	fPrompter(command.Prompter()),
	fOrdered(command.Ordered()),
	fArguments(0)
{
	/* copy argument list */
	const ArrayT<ArgSpecT*>& args = command.Arguments();
	for (int i = 0; i < args.Length(); i++)
		AddArgument(*(args[i]));
}

/* destructor */
CommandSpecT::~CommandSpecT(void)
{
	for (int i = 0; i < fArguments.Length(); i++)
		delete fArguments[i];
}

/* return a reference to a specific argument */
ArgSpecT& CommandSpecT::Argument(const char* name)
{
	int index = -1;
	for (int i = 0; index == -1 && i < fArguments.Length(); i++)
		if (fArguments[i]->Name() == name)
			index = i;
	if (index == -1) {
		cout << "CommandSpecT::Argument: not found \"" << name << '\"' << endl;
		throw ExceptionT::kGeneralFail;
	}
	return Argument(index);
}

/* return a reference to a specific argument */
const ArgSpecT& CommandSpecT::Argument(const char* name) const
{
	int index = -1;
	for (int i = 0; index == -1 && i < fArguments.Length(); i++)
		if (fArguments[i]->Name() == name)
			index = i;
	if (index == -1) {
		cout << "CommandSpecT::Argument: not found \"" << name << '\"' << endl;
		throw ExceptionT::kGeneralFail;
	}
	return Argument(index);
}

/* add an argument to the function */
void CommandSpecT::AddArgument(const ArgSpecT& arg)
{
	/* check - unnamed must be ordered */
	if (!fOrdered && arg.Name().StringLength() == 0)
	{
		cout << "\n CommandSpecT::AddArgument: unordered arguments must be named" << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* new argument spec */
	ArgSpecT* new_arg = new ArgSpecT(arg);

	/* store */
	fArguments.Append(new_arg);
}

/* clear all argument values */
void CommandSpecT::ClearValues(void)
{
	for (int i = 0; i < fArguments.Length(); i++)
	fArguments[i]->ClearValue();
}

/* write function spec to output stream */
void CommandSpecT::Write(ostream& out) const
{
	/* function name */
	out << fName << ": " << fArguments.Length();
	if (fArguments.Length() == 1)
		out << " argument";
	else
		out << " arguments";
		
	if (fArguments.Length() > 1)
	{
		if (fOrdered)
			out << ": ordered\n";
		else
			out << ": not ordered\n";
	}
	else out << '\n';	
		
	/* arguments */
	for (int i = 0; i < fArguments.Length(); i++)
	{
		fArguments[i]->Write(out);
		out << '\n';
	}

	out.flush();
}

/* write command statement to the output stream */
void CommandSpecT::WriteCommand(ostream& out) const
{
	/* function name */
	out << fName << " ";
	for (int i = 0; i < fArguments.Length(); i++)
	{
		/* argument name */
		if (fArguments[i]->Name().StringLength() > 0)
			out << fArguments[i]->Name() << " ";
			
		/* value */
		fArguments[i]->WriteValue(out);
		out << " ";
	}
}
