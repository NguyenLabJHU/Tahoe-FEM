/* $Id: iConsoleT.h,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (12/21/2000)                                          */
/* iConsoleT.h                                                            */

#ifndef _I_CONSOLE_T_H_
#define _I_CONSOLE_T_H_

/* base class */
#include "iConsoleBaseT.h"

/* direct members */
#include "fstreamT.h"

/* forward declaration */
class iConsoleObjectT;

class iConsoleT: public iConsoleBaseT
{
public:

	/* constructor */
	iConsoleT(const StringT& log_file, iConsoleObjectT& current);

	/* destructor */
	~iConsoleT(void);

	/* execute given command - returns false on fail */
	virtual bool iDoCommand(const StringT& command, StringT& line);

	/* operate on given variable */
	virtual bool iDoVariable(const StringT& variable, StringT& line);

private:

	/* main event loop */
	void DoInteractive(void);

	/* get command line */
	void GetCommandLine(StringT& line);

	/* change scope */
	void SetScope(iConsoleObjectT& scope);

	/* resolve scope pointer - returns NULL if not found */
	iConsoleObjectT* GetScope(iConsoleObjectT& start, StringT& line) const;

	/* pulls the first word from the line and resolves it into
	 * a command from the console or current scope, or returns
	 * kNone if the word could not be resolved */
	enum CommandScope {kNone = 0,
		     kConsoleCommand = 1,
	        kConsoleVariable = 2,
	           kScopeCommand = 3,
	          kScopeVariable = 4,
	                  kAlias = 5};
	CommandScope ResolveNextWord(StringT& line, StringT& command) const;
	CommandScope ResolveCommand(StringT& command) const;
	
	/* reset dictionary - scope_only sets only scope commands
	 * and variables */
	void BuildDictionary(bool scope_only);

private:

	/* commands */
	void ListCommand(ostream& out) const;
	bool HistoryCommand(StringT& line);
	
	/* flush the command line and all input streams */
	void FlushInput(StringT& line);

	/* make an alias - returns false on fail */
	bool MakeAlias(const StringT& alias, StringT& line);

	/* manipulating the history stack */
	void PushHistory(const StringT& line);
	void PopHistory(void);
	void TopHistory(StringT& line);

private:

	/* log file */
	ofstreamT flog;

	/* parameters */
	int fmax_recursion_depth;
	int fhistory_size;

	/* current console object */
	iConsoleObjectT* fCurrent;
	iConsoleObjectT* fLastCurrent;
		
	/* scope */
	StringT fScope;
	
	/* runtime */
	int frecursion_depth;
	AutoArrayT<ifstreamT*> fInputStack;
	AutoArrayT<StringT>    fDanglingInput;
	AutoArrayT<StringT*>   fHistory;
	

	/* dictionary */
	AutoArrayT<StringT*>     fWord;
	AutoArrayT<CommandScope> fWordScope;
	
	/* aliases */
	AutoArrayT<StringT> fAlias;
	AutoArrayT<StringT> fAliasCommand;
};

#endif /* _I_CONSOLE_T_H_ */
