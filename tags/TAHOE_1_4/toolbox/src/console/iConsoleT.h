/* $Id: iConsoleT.h,v 1.15 2004-06-17 06:38:16 paklein Exp $ */
/* created: paklein (12/21/2000) */
#ifndef _I_CONSOLE_T_H_
#define _I_CONSOLE_T_H_

/* base class */
#include "iConsoleBaseT.h"

/* direct members */
#include "ofstreamT.h"

namespace Tahoe {

/* forward declaration */
class iConsoleObjectT;

/** interactive console. */
class iConsoleT: public iConsoleBaseT
{
  public:

	/* constructor */
	iConsoleT(const StringT& log_file, iConsoleObjectT& current,
		const ArrayT<StringT>* arguments = NULL,
		bool do_interactive = true);

	/* destructor */
	virtual ~iConsoleT(void);

	/* execute given command - returns false on fail */
	virtual bool iDoCommand(const CommandSpecT& command, StringT& line);

	/* operate on given variable */
	virtual bool iDoVariable(const StringT& variable, StringT& line);

	/* console flags */
	enum CommandScope {kNone = 0,
		     kConsoleCommand = 1,
	        kConsoleVariable = 2,
	           kScopeCommand = 3,
	          kScopeVariable = 4,
	                  kAlias = 5};

	/** return the name of the current scope */
	const StringT& Scope(void) { return fScope; };
	
	/** reference to the object which is the current scope */
	iConsoleObjectT& Current(void) { return *fCurrent; };

	/** main event loop */
	void DoInteractive(void);

  private:

	/* get command line */
	void GetCommandLine(StringT& line);

	/* change scope */
	void SetScope(iConsoleObjectT& scope);

	/* resolve scope pointer - returns NULL if not found */
	iConsoleObjectT* GetScope(iConsoleObjectT& start, StringT& line) const;

	/* pulls the first word from the line and resolves it into
	 * a command from the console or current scope, or returns
	 * kNone if the word could not be resolved */
	CommandScope ResolveNextWord(StringT& line, StringT& command) const;
	CommandScope ResolveCommandName(StringT& command) const;
	
	/* reset dictionary - scope_only sets only scope commands
	 * and variables */
	void BuildDictionary(bool scope_only);

	/* commands */
	void ListCommand(ostream& out) const;
	
	/* flush the command line and all input streams */
	void FlushInput(StringT& line);

	/* make an alias - returns false on fail */
	bool MakeAlias(const StringT& alias, StringT& line);

	/** \name command history
	 * manipulating the history stack */
	/*@{*/
	void PushHistory(const StringT& line);
	void PopHistory(void);
	void TopHistory(StringT& line);
	/*@}*/

	/** pull the next command from the line. Commands are separated by ';',
	 * but separators contained in quoted strings are ignored */
	void NextCommand(const StringT& source, StringT& next) const;

private:

	/** log file */
	ofstreamT flog;

	/** console variables */
	/*@{*/
	int fmax_recursion_depth;
	int fhistory_size;
	/*@}*/

	/** current console object */
	/*@{*/
	iConsoleObjectT* fCurrent;
	iConsoleObjectT* fLastCurrent;
	/*@}*/
		
	/** scope */
	StringT fScope;
	
	/** runtime */
	/*@{*/
	int  frecursion_depth;
	bool fstop_read_on_error;
	AutoArrayT<ifstreamT*> fInputStack;
	AutoArrayT<StringT>    fDanglingInput;
	int fHistoryCount;
	AutoArrayT<StringT*>   fHistory;
	/*@}*/

	/** dictionary */
	/*@{*/
	AutoArrayT<const StringT*> fWord;
	AutoArrayT<CommandScope>   fWordScope;
	/*@}*/

	/** aliases */
	/*@{*/
	AutoArrayT<StringT> fAlias;
	AutoArrayT<StringT> fAliasCommand;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _I_CONSOLE_T_H_ */
