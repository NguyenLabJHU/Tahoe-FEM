/* $Id: iConsoleBaseT.h,v 1.13 2003-11-10 22:14:15 cjkimme Exp $ */
/* created: paklein (12/21/2000) */

#ifndef _I_CONSOLE_BASE_T_H_
#define _I_CONSOLE_BASE_T_H_

/* direct members */
#include "AutoArrayT.h"
#include "StringT.h"

namespace Tahoe {

/* forward declarations */
class CommandSpecT;
class ArgSpecT;

/** base class for interactive console and console objects */
class iConsoleBaseT
{
public:

	/** console variable types */
	enum VariableType {int_ = 0, double_ = 1, string_ = 2, bool_ = 3, float_ = 4};

	/** constructor */
	iConsoleBaseT(void);
	
	/** destructor */
	virtual ~iConsoleBaseT(void);

	/** command list */
	const ArrayT<CommandSpecT*>& iCommands(void) const;

	/** variable list */
	const ArrayT<StringT>& iVariables(void) const;

	/** write variables */
	virtual void iWriteVariables(ostream& out) const;

	/** execute given command.
	 * \return true if command executes normally, false otherwise */
	virtual bool iDoCommand(const CommandSpecT& command, StringT& line);

	/** operate on given variable
	 * \return true if executed normally, false otherwise */
	virtual bool iDoVariable(const StringT& variable, StringT& line);

	/** resolve name into function specification. Map function name onto
	 * the list of function specifications, resolve all function arguments
	 * and return a complete function specification. See CommandSpecT for
	 * more information about function specifications.
	 * \param command_name name of function to resolve
	 * \param line pointer command line to probe for function arguments. 
	 * \return a pointer command specification if the function and
	 *         all of its arguments were resolved without problems,
	 *         NULL otherwise. */
	const CommandSpecT* iResolveCommand(const StringT& command_name, 
		StringT& line) const;

	/** return the command specification with the given name. Returns
	 * NULL if the name is not found */
	CommandSpecT* iCommand(const StringT& command_name) const;

#ifdef _MSC_VER
	/** operator types */
	enum VariableOperator {kEQ = 0, kPlusEQ, kMinusEQ, kTimesEQ, kDivEQ, kFail};
#endif

protected:

	/** write prompt for the specific argument of the command
	 * to the output stream. A simple value prompt always appears 
	 * <i>after</i> any information written during this call. By
	 * default, no additional information is written to the output
	 * stream. */
	virtual void ValuePrompt(const CommandSpecT& command, int index, ostream& out) const;

	/** resolve command arguments. Look in line passed in for required
	 * function arguments. If not present in line, use default argument
	 * values. Otherwise, prompt for argument values interactively.
	 * \param command command specification 
	 * \param line command line to probe for arguments 
	 * \param out output stream for prompts
	 * \param in input stream for interactive input 
	 * \return true if all arguments resolved correctly, false otherwise */
	bool ResolveArguments(CommandSpecT& command, StringT& line, ostream& out, 
		istream& in) const;

	/** resolve named argument value. Just checks name and then uses 
	 * iConsoleBaseT::ResolveValue to resolve the argument. */
	bool ResolveNamedValue(CommandSpecT& command, int index, StringT& line, 
		ostream& out, istream& in, bool prompt) const;
		
	/** resolve unnamed argument value.
	 * \param command command being resolved
	 * \param index index of the argument being resolved
	 * \param line source string to probe for values
	 * \param out output stream for prompts
	 * \param in input stream for interactive input 
	 * \param prompt pass true to produce prompt if value not found in line */
	bool ResolveValue(CommandSpecT& command, int index, StringT& line, ostream& out, 
		istream& in, bool prompt) const;

	/** add command to the dictionary.
	 * \return true if added, false otherwise */
	bool iAddCommand(const CommandSpecT& command);
	
	/* adding variables */
	bool iAddVariable(const StringT& name, bool& variable);
	bool iAddVariable(const StringT& name, const bool& variable);

	bool iAddVariable(const StringT& name, int& variable);
	bool iAddVariable(const StringT& name, const int& variable);

	bool iAddVariable(const StringT& name, float& variable);
	bool iAddVariable(const StringT& name, const float& variable);

	bool iAddVariable(const StringT& name, double& variable);
	bool iAddVariable(const StringT& name, const double& variable);

	bool iAddVariable(const StringT& name, StringT& variable);
	bool iAddVariable(const StringT& name, const StringT& variable);
	
	/** remove named variable. \return true if found and removed,
	 * false otherwise */
	bool iDeleteVariable(const StringT& name);

	/** alphabetize the list */
	void Sort(ArrayT<StringT>& list) const;

	/** sort command list by command name */
	void SortCommands(ArrayT<CommandSpecT*>& list) const;

	/** write list of strings with tab and wrap */
	void WriteList(ostream& out, const ArrayT<StringT>& list, int tab,
		int wrap) const;

	/** write list of command names with tab and wrap */
	void WriteList(ostream& out, const ArrayT<CommandSpecT*>& list, int tab,
		int wrap) const;

	/** write single variable
	 * \param out output stream
	 * \param i index of the variable to write */
	void WriteVariable(ostream& out, int i) const;

#ifndef _MSC_VER
	/** operator types */
	enum VariableOperator {kEQ = 0, kPlusEQ, kMinusEQ, kTimesEQ, kDivEQ, kFail};
#endif

	/** resolve the operator type from the input line */
	VariableOperator ResolveOperator(StringT& line) const;

	/** copy variables from the source. \return true if all variables added. */
	bool AddVariables(const iConsoleBaseT& source);

	/** remove all variables */
	void DeleteVariables(void);
	
private:

	/** find first position.
	 * \return false if not found, true otherwise */
	bool Position(char* str, char a, int& position);

	/* add variable */
	bool AddVariable(const StringT& name, VariableType type, void* variable, bool is_const);

	/* variable operators - return false on fail */
	bool Operate(bool& variable, VariableOperator op, StringT& line) const;
	bool Operate(int& variable, VariableOperator op, StringT& line) const;
	bool Operate(float& variable, VariableOperator op, StringT& line) const;
	bool Operate(double& variable, VariableOperator op, StringT& line) const;
	bool Operate(StringT& variable, VariableOperator op, StringT& line) const;

protected:

	/** commands */
	AutoArrayT<CommandSpecT*> fCommands;

	/** variables */
	AutoArrayT<StringT>      fVariables;
	AutoArrayT<VariableType> fVariableTypes;
	AutoArrayT<void*>        fVariableValues;
	AutoArrayT<bool>         fVariableIsConst;
};

/* inlines */
inline const ArrayT<CommandSpecT*>& iConsoleBaseT::iCommands(void) const
{
	return fCommands;
}

inline const ArrayT<StringT>& iConsoleBaseT::iVariables(void) const
{
	return fVariables;
}

} // namespace Tahoe 
#endif /* _I_CONSOLE_BASE_T_H_ */
