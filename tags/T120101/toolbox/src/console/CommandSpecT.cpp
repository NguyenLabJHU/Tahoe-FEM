/* $Id: CommandSpecT.cpp,v 1.1 2001-11-28 22:05:44 paklein Exp $ */

#include "CommandSpecT.h"
#include "ArgSpecT.h"

/* array copy behavior */
const bool ArrayT<CommandSpecT*>::fByteCopy = true; 
const bool ArrayT<CommandSpecT>::fByteCopy = false; 

CommandSpecT::CommandSpecT(const StringT& name, bool ordered_args):
	fName(name),
	fOrdered(ordered_args),
	fArguments(0)
{

}

/* copy constructor */
CommandSpecT::CommandSpecT(const CommandSpecT& command):
	fName(command.Name()),
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

/* add an argument to the function */
void CommandSpecT::AddArgument(const ArgSpecT& arg)
{
	/* check - unnamed must be ordered */
	if (!fOrdered && arg.Name().StringLength() == 0)
	{
		cout << "\n CommandSpecT::AddArgument: unordered arguments must be named" << endl;
		throw eGeneralFail;
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
		out << " argument\n";
	else
		out << " arguments\n";
		
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
