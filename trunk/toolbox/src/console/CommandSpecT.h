/* $Id: CommandSpecT.h,v 1.6 2002-07-05 22:26:24 paklein Exp $ */

#ifndef _COMMAND_SPEC_T_H_
#define _COMMAND_SPEC_T_H_

/* direct members */
#include "StringT.h"
#include "AutoArrayT.h"

namespace Tahoe {

/* forward declarations */
class ArgSpecT;
class iConsoleBaseT;

/** definition of console commands */
class CommandSpecT
{
  public:

	/** constructor */
	explicit CommandSpecT(const StringT& name, bool ordered_args = true);

	/** copy constructor */
	CommandSpecT(const CommandSpecT& command);

	/** destructor */
	~CommandSpecT(void);

	/** return the function name */
	const StringT& Name(void) const { return fName; };

	/** return ordered arguments flag */
	bool Ordered(void) const { return fOrdered; };

	/** add an argument to the function */
	void AddArgument(const ArgSpecT& arg);
	
	/** return the argument list */
	const ArrayT<ArgSpecT*>& Arguments(void) const { return fArguments; };

	/** return a reference to a specific argument. The indeces
	 * of arguments is set by the order in which they are added
	 * using CommandSpecT::AddArgument. */
	ArgSpecT& Argument(int index) { return *fArguments[index]; };

	/** constant reference to command argument. The indeces
	 * of arguments is set by the order in which they are added
	 * using CommandSpecT::AddArgument. */
	const ArgSpecT& Argument(int index = 0) const { return *fArguments[index]; };

	/** return a reference to a specific argument. The indeces
	 * of arguments is set by the order in which they are added
	 * using CommandSpecT::AddArgument. Works with named arguments
	 * only. Throws exception if argument not found. */
	ArgSpecT& Argument(const char* name);

	/** return a const reference to a specific argument. The indeces
	 * of arguments is set by the order in which they are added
	 * using CommandSpecT::AddArgument. Works with named arguments
	 * only. Throws exception if argument not found. */
	const ArgSpecT& Argument(const char* name) const;
	
	/** set the argument prompter */
	void SetPrompter(const iConsoleBaseT* prompter) { fPrompter = prompter; };

	/** return the prompter */
	const iConsoleBaseT* Prompter(void) const { return fPrompter; };

	/** clear all argument values */
	void ClearValues(void);

	/** write function spec to output stream */
	void Write(ostream& out) const;

	/** write command statement to the output stream */
	void WriteCommand(ostream& out) const;

  private:
  
  	/** private to disallow */
  	CommandSpecT& operator=(const CommandSpecT&);
	
  private:
  
  	/** command name */
  	const StringT fName;

	/** command prompter */
	const iConsoleBaseT* fPrompter;
  	
  	/** true if function arguments have fixed order */
  	const bool fOrdered;
  	
  	/** list of function arguments */
  	AutoArrayT<ArgSpecT*> fArguments;
};

} // namespace Tahoe 
#endif /* _COMMAND_SPEC_T_H_ */
