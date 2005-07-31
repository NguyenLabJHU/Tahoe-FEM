/* $Id: iConsoleT.cpp,v 1.1.1.1 2001-01-25 20:56:28 paklein Exp $ */
/* created: paklein (12/21/2000)                                          */
/* iConsoleT.cpp                                                          */

#include "iConsoleT.h"

#include <ctype.h>
#include <time.h>
#include <iomanip.h>

#include "iConsoleObjectT.h"
#include "ifstreamT.h"

/* array behavior */
const bool ArrayT<iConsoleT::CommandScope>::fByteCopy = true;

/* constructor */
iConsoleT::iConsoleT(const StringT& log_file, iConsoleObjectT& current):
	flog(log_file, true),
	fmax_recursion_depth(25),
	fhistory_size(10),
	fCurrent(NULL),
	fInputStack(0, true),
	fHistory(fhistory_size + 1, 0), /* shift size by 1 */
	
	/* dictionary */
	fWord(20, true),
	fWordScope(20, true),

	/* aliases */
	fAlias(20, false),
	fAliasCommand(20, false)
{
	/* set commands */
	iAddCommand("end");
	iAddCommand("scope");
	iAddCommand("help");
	iAddCommand("list"); // add option to list N levels of heirarchy?
	iAddCommand("repeat");
	iAddCommand("read");
	iAddCommand("back");
	iAddCommand("alias");
	iAddCommand("history");
	
	/* set variables */
	iAddVariable("max_recursion_depth", fmax_recursion_depth);
	iAddVariable("history_size", fhistory_size);
		
	/* initialize dictionary with scope */
	BuildDictionary(false);

	/* set scope */
	SetScope(current);
	
	/* add default alias */
	if (fCurrent->iSuper() != NULL)
	{
		/* build absolute path to current scope */
		StringT home_path;
		const iConsoleObjectT* obj = fCurrent;
		while (obj->iSuper() != NULL)
		{
			home_path.Prepend(":", obj->iName());
			obj = obj->iSuper();
		}
		home_path.Prepend(":root");
		
		/* add alias */
		StringT alias("home");
		MakeAlias(alias, home_path);
	}
	
	/* clear history */
	fHistory = NULL;
	
	/* run */
	DoInteractive();
}

/* destructor */
iConsoleT::~iConsoleT(void)
{
	/* free remaining streams */
	for (int i = 0; i < fInputStack.Length(); i++)
		delete fInputStack[i];

	/* free history */
	for (int j = 0; j < fHistory.Length(); j++)
		delete fHistory[j];
}

/* execute given command */
bool iConsoleT::iDoCommand(const StringT& command, StringT& line)
{
	/* dispatch */
	if (command == "scope")
	{
		cout << fScope << endl;
		return true;
	}
	else if (command == "help")
	{
		cout << "console commands:\n";
		WriteList(cout, fCommands, 4, 80);
		cout << '\n' << setw(4) << " " << ":[scope] :root" << '\n';
		cout << "scope commands:\n";
		const ArrayT<StringT>& commands = fCurrent->iCommands();
		if (commands.Length() > 0)
		{
			WriteList(cout, commands, 4, 80);
			cout << '\n';
		}
		else
			cout << setw(4) << " " << "<none>\n";
		cout << "variable operators:\n";
		cout << setw(4) << " " << "= += -= *= /=\n";
		cout << "aliases:\n";
		if (fAlias.Length() == 0)
			cout << setw(4) << " " << "<none>\n";		
		else
		{
			for (int i = 0; i < fAlias.Length(); i++)
				cout << setw(4) << " " << fAlias[i] << " -> "
				     << fAliasCommand[i] << '\n';
		}
		cout.flush();
		return true;
	}
	else if (command == "list")
	{
		ListCommand(cout);
		return true;
	}
	else if (command == "repeat")
	{
		/* resolve argument */
		int repeat = 1;
		if (line[0] == '(')
			if (!ResolveArgument(line, repeat, &repeat))
			{
				cout << "could not resolve integer argument from: \"" << line
				     << "\"" << endl;
				line.Drop(strlen(line));
				return false;
			}
			
		/* execute */	
		if (repeat > 1)
		{
			StringT source = line;
			int count = 1;
			while (count++ < repeat)
				line.Append(" ", source);			
		}
		return true;
	}
	else if (command == "read")
	{
		StringT file_name;
		if (!ResolveArgument(line, file_name, NULL))
		{
			cout << "could not resolve string argument from: \"" << line
			     << "\"" << endl;
			line.Drop(strlen(line));
			return false;
		}
		else
		{
			ifstreamT* new_stream = new ifstreamT('#');
			new_stream->open(file_name);
			if (new_stream->is_open())
			{
				cout << "opened input stream: \"" << file_name << "\"" << endl;
				fInputStack.Append(new_stream);
				
				/* store rest of line */
				fDanglingInput.Append(line);
				line.Drop(strlen(line));
				return true;
			}
			else
			{
				cout << "could not open input stream: \"" << file_name << "\"" << endl;
				delete new_stream;
				line.Drop(strlen(line));
				return false;
			}
		}
	}
	else if (command == "back")
	{
		if (fLastCurrent != NULL)
			SetScope(*fLastCurrent);
		return true;
	}
	else if (command == "alias")
	{
		int count;
		StringT alias;
		alias.FirstWord(line, count, true);
		
		/* empty */
		if (strlen(alias) == 0)
			return true;
		else
		{
			/* pull alias name */
			line.Drop(count);
			
			/* check for new definition */
			if (strncmp(line, "->", 2) == 0)
			{
				/* shift */
				line.Drop(2);
				line.DropLeadingSpace();
				
				return MakeAlias(alias, line);
			}
			/* looking for alias definition */
			else
			{
				int index = fAlias.PositionOf(alias);
				if (index == -1)
					cout << "alias \"" << alias << "\" not defined" << endl;
				else
					cout << alias << " -> " << fAliasCommand[index] << endl;
				return true;
			}
		}
	}
	else if (command == "history")
		return HistoryCommand(line);
	else
		return iConsoleBaseT::iDoCommand(command, line);
}

/* operate on given variable */
bool iConsoleT::iDoVariable(const StringT& variable, StringT& line)
{
	/* inherited */
	bool result = iConsoleBaseT::iDoVariable(variable, line);
	if (result)
	{
		/* safe resizing */
		if (variable == "history_size")
		{
			/* shift size */
			fhistory_size++;
		
			if (fhistory_size > fHistory.Length())
				fHistory.Resize(fhistory_size, NULL);
			else if (fhistory_size < fHistory.Length())
			{
				for (int i = fhistory_size; i < fHistory.Length(); i++)
				{
					delete fHistory[i];
					fHistory[i] = NULL;
				}
				fHistory.Resize(fhistory_size);
			}

			/* shift back */
			fhistory_size--;
		}
	}
	return result;
}

/************************************************************************
* Protected
************************************************************************/

/* main event loop */
void iConsoleT::DoInteractive(void)
{
	/* open log stream */
	time_t the_time;
	time(&the_time);
	flog << "\n###################################################\n"
	     << "# open: " << ctime(&the_time)
	     << "###################################################\n";
	
	StringT line_copy, line, command;
	GetCommandLine(line);
	bool end = false;
	bool line_OK = true;
	while (!end)
	{
		/* shift */
		line.DropLeadingSpace();

		/* log (but not if read from external) */
		if (fInputStack.Length() == 0) flog << line << '\n';
		
		/* consume command line */
		frecursion_depth = 0;
		while (strlen(line) > 0 && !end)
		{
			/* fail safe */
			if (frecursion_depth++ > fmax_recursion_depth)
			{
				cout << "exceeded maximum recursion depth: " << fmax_recursion_depth << endl;
				line_OK = false;
			}
			
			/* scope change vs command */
			else if (line[0] == ':')
			{
				/* take scope specifier */
				int count;
				StringT scope_line;
				scope_line.FirstWord(line, count, false);
				line.Drop(count);
			
				/* try scope change */
				iConsoleObjectT* scope = GetScope(*fCurrent, scope_line);
				
				/* change scope */
				if (scope != NULL)
					SetScope(*scope);
				else
					line_OK = false;
			}
			else
			{
				/* resolve next word in line */
				CommandScope scope = ResolveNextWord(line, command);
				
				/* dispatch */
				switch (scope)
				{
					case kConsoleCommand:
					{
						if (command == "end")
							end = true;
						else	
							line_OK = iDoCommand(command, line);
						break;
					}					
					case kConsoleVariable:
						line_OK = iDoVariable(command, line);
						break;
					case kScopeCommand:
						line_OK = fCurrent->iDoCommand(command, line);
						break;
					case kScopeVariable:
						line_OK = fCurrent->iDoVariable(command, line);
						break;
					case kAlias:
					{
						/* locate */
						int index = fAlias.PositionOf(command);
						if (index == -1)
							throw eGeneralFail;
						else
							line.Prepend(fAliasCommand[index], " ");
						break;
					}
					default:
					{
						/* message */
						if (strlen(line) > 0)
							cout << "unresolved: \"" << line << "\"" << endl;
					
						/* not OK */
						line_OK = false;
					}
				}
			}
			
			/* check remaining command line */
			if (line_OK)
				line.DropLeadingSpace();
			else
				FlushInput(line);
		}
		
		/* manage history */
		if (!line_OK) PopHistory();
		
		/* read next command */
		if (!end)
		{
			GetCommandLine(line);
			line_copy = line;
			line_OK = true;
		}
	}

	flog << "###################################################\n"
	     << "# close: " << ctime(&the_time)
	     << "###################################################\n";
}

/* get command line */
void iConsoleT::GetCommandLine(StringT& line)
{
	bool done = false;
	while (!done)
	{
		if (fInputStack.Length() == 0)
		{
			/* prompt */
			cout << fCurrent->iName();
			if (fHistory.Length() > 0 && fHistory[0] != NULL)
				cout << ": " << *(fHistory[0]) << " ";
			cout << "> ";
			line.GetLineFromStream(cin);
			
			/* store history */
			if (strlen(line) > 0)
				PushHistory(line);
			/* retrieve */
			else
				TopHistory(line);

			/* exit */
			done = true;
		}
		else
		{
			ifstreamT* in = fInputStack.Last();
			if (in->good())
			{
				line.GetLineFromStream(*in);
				done = true;
			}
			else
			{
				cout << "closing stream: \"" << in->filename()  << "\"" << endl;
				
				/* free stream */
				delete in;
				in = NULL;
				fInputStack.DeleteAt(fInputStack.Length() - 1);
				
				/* restore rest of line */
				line = fDanglingInput.Last();
				fDanglingInput.DeleteAt(fDanglingInput.Length() - 1);
				done = true;
			}
		}
	}
}

/* change scope */
void iConsoleT::SetScope(iConsoleObjectT& scope)
{
	/* store last */
	fLastCurrent = fCurrent;

	/* scope pointer */
	fCurrent = &scope;

	/* build scope name */
	fScope = fCurrent->iName();
	const iConsoleObjectT* obj = fCurrent;
	while (obj->iSuper() != NULL)
	{
		/* up */
		obj = obj->iSuper();
	
		/* build string */
		fScope.Prepend(obj->iName(), ":");
	}
	
	/* set dictionary (scope symbols only) */
	BuildDictionary(true);
}

/* resolve scope pointer - returns NULL if not found */
iConsoleObjectT* iConsoleT::GetScope(iConsoleObjectT& start,
	StringT& line) const
{
	/* drop leading ":" */
	if (strlen(line) > 0) line.Drop(1);

	/* done */
	if (strlen(line) == 0)
		return &start;
	/* up */
	else if (line[0] == ':')
	{
		/* recurse */
		if (start.iSuper() != NULL)
			return GetScope(*(start.iSuper()), line);
		else
		{
			cout << "could not resolve scope: \""
			     << line << "\"" << endl;
			return NULL;
		}
	}
	/* go to root */
	else if (strncmp(line, "root", 4) == 0)
	{
		/* pull "root" */
		line.Drop(4);
		
		/* find root */
		iConsoleObjectT* scope = &start;
		while (scope->iSuper() != NULL)
			scope = scope->iSuper();
		
		/* recurse */
		return GetScope(*scope, line);
	}
	else /* down */
	{
		/* partial scope length */
		int length = 0;
		char* str = line;
		int max = strlen(str);
		while (!isspace(*str) && *str != ':' && length < max)
		{
			length++;
			str++;
		}
		
		/* sub-objects */
		const ArrayT<iConsoleObjectT*>& subs = start.iSubs();
		int match_count = 0;
		iConsoleObjectT* match = NULL;
		for (int i = 0; i < subs.Length(); i++)
		{
			/* match */
			if (strncmp(line, subs[i]->iName(), length) == 0)
			{
				/* keep match */
				match = subs[i];

				/* exact match */
				if (strlen(subs[i]->iName()) == length)
				{
					match_count = 1;
					break;
				}		
				else
					match_count++;
			}
		}

		/* match */
		if (match_count == 1)
		{
			line.Drop(length);
			if (strlen(line) == 0)
				return match;
			else	
				return GetScope(*match, line);
		}
		/* ambiguous */
		else if (match_count > 1)
		{
			cout << "multiple completions (" << match_count << "):\n";
			int count = 0;
			for (int i = 0; i < subs.Length(); i++)
				if (strncmp(line, subs[i]->iName(), length) == 0)
					cout << setw(5) << ++count << ": "
					     << subs[i]->iName() << '\n';
			cout.flush();
			return NULL;
		}
		/* no match */
		else
		{
			cout << "could not resolve scope: \""
			     << line << "\"" << endl;
			return NULL;
		}
	}
}

iConsoleT::CommandScope iConsoleT::ResolveNextWord(StringT& line,
	StringT& command) const
{
	/* pull first word */
	int count;
	command.FirstWord(line, count, true);
	line.Drop(count);

	/* resolve */
	return ResolveCommand(command);
}

/* pulls the first word from the line and resolves it into
* a command from the console or current scope, or returns
* kNone if the word could not be resolved */
iConsoleT::CommandScope iConsoleT::ResolveCommand(StringT& command) const
{
	int word_length = strlen(command);
	if (word_length == 0) return kNone;

	const StringT* match = NULL;
	CommandScope scope = kNone;
	int match_count = 0;
	int exact_match_count = 0;
	bool stored_exact = false;
	for (int i = 0; i < fWord.Length(); i++)
	{
		/* partial match */
		if (strncmp(command, *(fWord[i]), word_length) == 0)
		{
			match_count++;
			
			/* exact match */
			if (strlen(*(fWord[i])) == word_length)
				exact_match_count++;
				
			/* store */
			if ((match_count == 1 || exact_match_count == 1) && !stored_exact)
			{
				command = *(fWord[i]);
				scope = fWordScope[i];
				if (exact_match_count == 1) stored_exact = true;
			}						
		}
	}
	
	/* no match */
	if (match_count == 0)
	{
		cout << "unrecognized command: \"" << command << "\"" << endl;
		command = "";
		scope = kNone;
	}
	/* multiple matches */
	else if (match_count > 1)
	{
		const char* scope_names[] = {"none",
		                  "console command",
		                 "console variable",
		                    "scope command",
		                   "scope variable",
		                            "alias"};
	
		/* no exact matches */
		if (exact_match_count == 0)
		{
			cout << "multiple completions (" << match_count << "):\n";
			int count = 0;
			for (int i = 0; i < fWord.Length(); i++)
				if (strncmp(command, *(fWord[i]), word_length) == 0)
					cout << setw(5) << ++count
					     << ": " << *(fWord[i]) << ": "
					     << scope_names[fWordScope[i]] << '\n';
			cout.flush();
			scope = kNone;
		}
		/* multiple exact matches */
		else if (exact_match_count > 1)
		{
			cout << "multiple exact matches (" << exact_match_count << "):\n";
			ArrayT<StringT*> exact_match_command(exact_match_count);
			ArrayT<CommandScope> exact_match_scope(exact_match_count);
			exact_match_count = 0;
			for (int i = 0; i < fWord.Length(); i++)
			{
				StringT* pcommand = fWord[i];
				if (command == *pcommand)
				{
					/* write */
					cout << setw(5) << exact_match_count + 1 << ": "
					     << scope_names[fWordScope[i]] << '\n';
				
					/* store */
					exact_match_command[exact_match_count] = pcommand;
					exact_match_scope[exact_match_count] = fWordScope[i];
					exact_match_count++;
				}
			}
			
			/* get user selection */
			int select = -99;
			while (select < 1 || select > exact_match_command.Length())
			{
				cout << " select: ";
				cin >> select;
				if (select == -99)
					cin.clear();
				else if (select < 1 || select > exact_match_command.Length())
					select = -99;
			}
			
			/* resolve */
			select--;
			command = *(exact_match_command[select]);
			scope = exact_match_scope[select];
			
			/* clear (new)line */
			char tmp[255];
			cin.getline(tmp, 254);
		}
	}
	return scope;
}

/* reset dictionary - scope_only sets only scope commands
* and variables */
void iConsoleT::BuildDictionary(bool scope_only)
{
	/* console symbols */
	const ArrayT<StringT>& console_commands = iCommands();
	const ArrayT<StringT>& console_variables = iVariables();
		
	/* reset console symbols */
	if (!scope_only)
	{
		/* resize */
		fWord.Allocate(0);
		fWordScope.Allocate(0);
	
		/* add commands */
		for (int i = 0; i < console_commands.Length(); i++)
		{
			fWord.Append(console_commands.Pointer(i));
			fWordScope.Append(kConsoleCommand);
		}
		
		/* add variables */
		for (int j = 0; j < console_variables.Length(); j++)
		{
			fWord.Append(console_variables.Pointer(j));
			fWordScope.Append(kConsoleVariable);
		}
		
		/* add aliases */
		for (int k = 0; k < fAlias.Length(); k++)
		{
			fWord.Append(fAlias.Pointer(k));
			fWordScope.Append(kAlias);
		}
	}
	else
	{
		/* resize */
		int console_symbols = console_commands.Length() +
		                      console_variables.Length() +
		                      fAlias.Length();
		fWord.Resize(console_symbols);
		fWordScope.Resize(console_symbols);	
	}

	/* scope symbols */
	if (fCurrent != NULL)
	{
		const ArrayT<StringT>& scope_commands = fCurrent->iCommands();
		const ArrayT<StringT>& scope_variables = fCurrent->iVariables();
	
		/* add commands */
		for (int i = 0; i < scope_commands.Length(); i++)
		{
			fWord.Append(scope_commands.Pointer(i));
			fWordScope.Append(kScopeCommand);
		}
		
		/* add variables */
		for (int j = 0; j < scope_variables.Length(); j++)
		{
			fWord.Append(scope_variables.Pointer(j));
			fWordScope.Append(kScopeVariable);
		}
	}
}

/************************************************************************
* Private
************************************************************************/

/* commands */
void iConsoleT::ListCommand(ostream& out) const
{
	/* list console variables */
	out << "console variables:\n";
	iWriteVariables(out);

	/* list scopes */
	out << "scopes:\n";
	const ArrayT<iConsoleObjectT*>& subs = fCurrent->iSubs();
	if (subs.Length() == 0)
		out << setw(4) << " " << "<none>" << endl;
	else
	{
		for (int i = 0; i < subs.Length(); i++)
			out << setw(5) << i+1 << ": "
			     << subs[i]->iName() << '\n';
	}
	
	/* list scope variables */
	out << "scope variables:\n";
	fCurrent->iWriteVariables(out);

	/* flush stream */
	out.flush();
}

bool iConsoleT::HistoryCommand(StringT& line)
{
	/* resolve argument */
	int selection = -99199199;
	if (line[0] == '(')
		if (!ResolveArgument(line, selection, &selection))
		{
			cout << "could not resolve integer argument from: \"" << line
			     << "\"" << endl;
			line.Drop(strlen(line));
			return false;
		}

	/* list history */
	if (selection == -99199199)
	{
		for (int i = fHistory.Length() - 2; i > -1; i--)
			if (fHistory[i] != NULL)
				cout << "history(" << i << ") -> " << *(fHistory[i]) << '\n';
		cout.flush();
		return true;
	}
	/* range check */
	else if (selection < 0 || selection >= fHistory.Length() - 1)
	{
		cout << "out of range" << endl;
		return false;	
	}
	/* repeat command from history */
	else
	{
		StringT* command = fHistory[selection + 1];
		if (command != NULL)
		{
			line.Prepend(*command, " ");
			PushHistory(line);
		}
		return true;
	}
}

/* flush the command line and all input streams */
void iConsoleT::FlushInput(StringT& line)
{
	/* flush line */
	line.Drop(strlen(line));

	/* flush input streams */
	for (int i = 0; i < fInputStack.Length(); i++)
	{
		cout << "closing stream: \"" << fInputStack[i]->filename()  << "\"" << endl;
		delete fInputStack[i];
		fInputStack[i] = NULL;
	}
	fInputStack.Allocate(0);
	fDanglingInput.Allocate(0);
}

/* make an alias - returns false on fail */
bool iConsoleT::MakeAlias(const StringT& alias, StringT& line)
{				
	/* remove alias */
	if (strlen(line) == 0)
	{
		int index = fAlias.PositionOf(alias);
		if (index == -1) throw eGeneralFail;
		
		/* remove */
		fAlias.DeleteAt(index);
		fAliasCommand.DeleteAt(index);
		
		/* message */
		cout << alias << " ->" << endl;
	}
	/* define new alias */
	else
	{
		/* add definition */
		if (!fAlias.AppendUnique(alias))
		{
			cout << "alias \"" << alias << "\" already defined" << endl;
			return false;
		}
		else
		{
			fAliasCommand.Append(line);
			line.Drop(strlen(line));

			/* message */
			cout << alias << " -> " << fAliasCommand.Last() << endl;
		}
	}	
	
	/* rebuild dictionary */
	BuildDictionary(false);
	return true;
}

/* manipulating the history stack */
void iConsoleT::PushHistory(const StringT& line)
{
	if (fHistory.Length() == 0)
		return;
	else if (fHistory[0] == NULL || *(fHistory[0]) != line)
	{
		/* free tail */
		StringT** last = &(fHistory.Last());
		if (*last != NULL)
		{
			delete *last;
			*last = NULL;
		}
		
		/* shorten */
		fHistory.Resize(fHistory.Length() - 1);

		/* add string to stack */
		StringT* new_line = new StringT(line);
		new_line->DropTrailingSpace();
		fHistory.Push(new_line);
	}
}

void iConsoleT::PopHistory(void)
{
	if (fHistory.Length() > 0 && fHistory[0] != NULL)
	{
		/* free */
		delete fHistory[0];
		fHistory[0] = NULL;
	
		/* shift down */
		fHistory.Pop();

		/* fill end */
		fHistory.Append(NULL);
	}
}

void iConsoleT::TopHistory(StringT& line)
{
	if (fHistory.Length() == 0)
		return;
	else if (fHistory[0] != NULL)
		line = *(fHistory[0]);

}
