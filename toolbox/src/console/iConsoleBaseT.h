/* $Id: iConsoleBaseT.h,v 1.1.1.1 2001-01-25 20:56:28 paklein Exp $ */
/* created: paklein (12/21/2000)                                          */
/* iConsoleBaseT.h                                                        */

#ifndef _I_CONSOLE_BASE_T_H_
#define _I_CONSOLE_BASE_T_H_

/* direct members */
#include "AutoArrayT.h"
#include "StringT.h"

class iConsoleBaseT
{
public:

	/* constructor */
	iConsoleBaseT(void);

	/* command list */
	const ArrayT<StringT>& iCommands(void) const;

	/* variable specifications */
	enum VariableType {int_ = 0, double_ = 1, string_ = 2, bool_ = 3};
	const ArrayT<StringT>& iVariables(void) const;

	/* write variables */
	virtual void iWriteVariables(ostream& out) const;

	/* execute given command - returns false on fail */
	virtual bool iDoCommand(const StringT& command, StringT& line);

	/* operate on given variable */
	virtual bool iDoVariable(const StringT& variable, StringT& line);

protected:

	/* add command to the dictionary - true if added */
	bool iAddCommand(const StringT& command);
	
	/* adding variables */
	bool iAddVariable(const StringT& name, bool& variable);
	bool iAddVariable(const StringT& name, const bool& variable);

	bool iAddVariable(const StringT& name, int& variable);
	bool iAddVariable(const StringT& name, const int& variable);

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

#endif /* _I_CONSOLE_BASE_T_H_ */
