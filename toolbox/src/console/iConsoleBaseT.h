/* $Id: iConsoleBaseT.h,v 1.3 2001-11-07 02:33:23 paklein Exp $ */
/* created: paklein (12/21/2000) */

#ifndef _I_CONSOLE_BASE_T_H_
#define _I_CONSOLE_BASE_T_H_

/* direct members */
#include "AutoArrayT.h"
#include "StringT.h"

/** base class for interactive console and console objects */
class iConsoleBaseT
{
public:

	/* constructor */
	iConsoleBaseT(void);

	/* command list */
	const ArrayT<StringT>& iCommands(void) const;

	/* variable specifications */
	enum VariableType {int_ = 0, double_ = 1, string_ = 2, bool_ = 3, float_ = 4};
	const ArrayT<StringT>& iVariables(void) const;

	/* write variables */
	virtual void iWriteVariables(ostream& out) const;

	/* execute given command - returns false on fail */
	virtual bool iDoCommand(const StringT& command, StringT& line);

	/* operate on given variable */
	virtual bool iDoVariable(const StringT& variable, StringT& line);

protected:

	/** clear the input stream. Remove the next 254 characters from the
	 * stream including any trailing newline. This is useful for clearing
	 * any leftovers from the command line when values are read using
	 * >>, which does not grab the trailing newline. */
	void Clean(istream& in) const;

	/* add command to the dictionary - true if added */
	bool iAddCommand(const StringT& command);
	
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

	/* alphabetize the list */
	void Sort(ArrayT<StringT>& list) const;

	/* resolving arguments in () */
	bool ResolveArgument(StringT& source, bool& arg, bool* default_arg);
	bool ResolveArgument(StringT& source, int& arg, int* default_arg);
	bool ResolveArgument(StringT& source, double& arg, double* default_arg);
	bool ResolveArgument(StringT& source, StringT& arg, StringT* default_arg);
//	bool ResolveArgument(StringT& source, ArrayT<ArgumentType>& types, pArrayT<void*>& args);

	/* write list of strings with tab and wrap */
	void WriteList(ostream& out, const ArrayT<StringT>& list, int tab,
		int wrap) const;

	/* write single variable */
	void WriteVariable(ostream& out, int i) const;

	/* operators */
	enum VariableOperator {kEQ = 0, kPlusEQ, kMinusEQ, kTimesEQ, kDivEQ, kFail};
	VariableOperator ResolveOperator(StringT& line) const;

private:

	/* find first position - returns false if not found */
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

	/* commands */
	AutoArrayT<StringT> fCommands;

	/* variables */
	AutoArrayT<StringT>      fVariables;
	AutoArrayT<VariableType> fVariableTypes;
	AutoArrayT<void*>        fVariableValues;
	AutoArrayT<bool>         fVariableIsConst;
};

/* inlines */
inline const ArrayT<StringT>& iConsoleBaseT::iCommands(void) const
{
	return fCommands;
}

inline const ArrayT<StringT>& iConsoleBaseT::iVariables(void) const
{
	return fVariables;
}

/* inlines */
inline void iConsoleBaseT::Clean(istream& in) const
{
  char line[255];
  in.getline(line, 254);
}

#endif /* _I_CONSOLE_BASE_T_H_ */
