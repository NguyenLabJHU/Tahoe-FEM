/* $Id: iConsoleT.cpp,v 1.27 2003-11-10 22:14:15 cjkimme Exp $ */
/* created: paklein (12/21/2000) */
#include "iConsoleT.h"

#include <ctype.h>
#include <time.h>
#include <iomanip.h>
#ifdef _MSC_VER
#include <strstrea.h>
#else
#include <strstream.h>
#endif

#include "iConsoleObjectT.h"
#include "ifstreamT.h"
#include "CommandSpecT.h"
#include "ArgSpecT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<iConsoleT::CommandScope>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
iConsoleT::iConsoleT(const StringT& log_file, iConsoleObjectT& current,
	const ArrayT<StringT>* arguments,
	bool do_interactive):
	flog(log_file, true),
	fmax_recursion_depth(25),
	fhistory_size(10),
	fCurrent(NULL),
	fInputStack(0),
	fHistoryCount(0),
	fHistory(fhistory_size + 1, 0), /* shift size by 1 */
	
	/* dictionary */
	fWord(20),
	fWordScope(20),

	/* aliases */
	fAlias(20),
	fAliasCommand(20)
{
	/* check log file is open */
	if (!flog.is_open()) {
		cout << "\n iConsoleT::iConsoleT: unable to open log file: \"" << log_file << "\"\n"
		     <<   "     will not log commands\n" << endl;
	}

	/* set commands */
	CommandSpecT the_end("end");
	ArgSpecT confirm(ArgSpecT::string_);
	confirm.SetPrompt("really quit (y/n)? ");
	the_end.AddArgument(confirm);
	iAddCommand(the_end);

	iAddCommand(CommandSpecT("scope"));
	iAddCommand(CommandSpecT("list"));
	iAddCommand(CommandSpecT("history"));
	iAddCommand(CommandSpecT("back"));

	CommandSpecT repeat_command("repeat");
	ArgSpecT repeat_count(ArgSpecT::int_);
	repeat_count.SetPrompt("repeat count");
	repeat_command.AddArgument(repeat_count);
	iAddCommand(repeat_command);

	CommandSpecT read_command("read");
	ArgSpecT read_file(ArgSpecT::string_);
	read_file.SetPrompt("command file name");
	ArgSpecT stop_on_error(ArgSpecT::bool_);
	stop_on_error.SetPrompt("terminate read on error");
	stop_on_error.SetDefault(true);
	read_command.AddArgument(read_file);
	read_command.AddArgument(stop_on_error);
	iAddCommand(read_command);

	CommandSpecT alias_command("alias");
	ArgSpecT alias_name(ArgSpecT::string_);
	alias_name.SetPrompt("alias name");
	alias_command.AddArgument(alias_name);
	ArgSpecT alias_target(ArgSpecT::string_);
	alias_target.SetPrompt("target of alias, <DELETE> to remove");
	alias_target.SetDefault("<ShowMe>");
	alias_command.AddArgument(alias_target);	
	iAddCommand(alias_command);

	CommandSpecT help_command("help");
	ArgSpecT arg(ArgSpecT::string_);
	arg.SetDefault("all");
	help_command.AddArgument(arg);
	iAddCommand(help_command);

	CommandSpecT wait_command("wait");
	ArgSpecT delay(ArgSpecT::double_);
	delay.SetDefault(0.0);
	delay.SetPrompt("pause length in seconds");
	wait_command.AddArgument(delay);
	iAddCommand(wait_command);
	
	CommandSpecT echo_command("echo");
	ArgSpecT mesg(ArgSpecT::string_);
	mesg.SetDefault(" ");
	mesg.SetPrompt("message to echo to standard out");
	echo_command.AddArgument(mesg);
	iAddCommand(echo_command);

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
	
	/* pull read command line arguments */
	StringT command_line;
	if (arguments) {
		for (int i = arguments->Length() - 2;  i > -1; i--)
			if ((*arguments)[i] == "-r")
			{
				StringT command;
				command.Append("read ", (*arguments)[i+1], ";");	
				fDanglingInput.Append(command);
			}
	}
	
	/* run */
	if (do_interactive) DoInteractive();
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
bool iConsoleT::iDoCommand(const CommandSpecT& command, StringT& line)
{
	/* dispatch */
	if (command.Name() == (const char*) "scope")
	{
		cout << fScope << endl;
		return true;
	}
	else if (command.Name() == "help")
	{
		StringT arg;
		command.Argument(0).GetValue(arg);
		
		/* write help */
		if (arg == "all") 
		{
			cout << "console commands:\n";
			WriteList(cout, fCommands, 4, 80);
			cout << '\n' << setw(4) << " " << ":[scope] :root" << '\n';
			cout << "scope commands:\n";
			const ArrayT<CommandSpecT*>& commands = fCurrent->iCommands();
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
					cout << setw(4) << " " << fAlias[i] << " -> \""
					     << fAliasCommand[i] << "\"\n";
			}
			cout.flush();
			return true;
		}
		else
		{
			/* resolve next word in line */
			StringT command_name;
			CommandScope scope = ResolveNextWord(arg, command_name);
				
			/* dispatch */
			switch (scope)
			{
				case kConsoleCommand:
				{
					int dex = -1;
					for (int i = 0; dex == -1 && i < fCommands.Length(); i++)
						if (fCommands[i]->Name() == command_name)
							dex = i;
					if (dex == -1) return false;

					cout << "console command:\n";
					fCommands[dex]->Write(cout);
					break;
				}
				case kScopeCommand:
				{
					/* scope command list */
					const ArrayT<CommandSpecT*>& commands = fCurrent->iCommands();
					int dex = -1;
					for (int i = 0; dex == -1 && i < commands.Length(); i++)
						if (commands[i]->Name() == command_name)
							dex = i;
					if (dex == -1) return false;

					cout << "scope command:\n";
					commands[dex]->Write(cout);
					break;
				}
				case kConsoleVariable:
					cout << "console variable" << endl;
					break;
				case kScopeVariable:
					cout << "scope variable" << endl;
					break;
				case kAlias:
					cout << "scope variable" << endl;
					break;
				default:
					cout << "could not resolve \"" << arg << '\"' << endl;
					return false;
			}
			return true;
		}
	}
	else if (command.Name() == "list")
	{
		ListCommand(cout);
		return true;
	}
	else if (command.Name() == "repeat")
	{
		int repeat;
		command.Argument(0).GetValue(repeat);

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
	else if (command.Name() == "read")
	{
		StringT file_name;
		command.Argument(0).GetValue(file_name);
		command.Argument(1).GetValue(fstop_read_on_error);		

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
	else if (command.Name() == "back")
	{
		if (fLastCurrent != NULL)
			SetScope(*fLastCurrent);
		return true;
	}
	else if (command.Name() == "alias")
	{
		StringT alias, target;
		command.Argument(0).GetValue(alias);
		command.Argument(1).GetValue(target);
		
		/* show/create/delete alias */
		if (target == "<ShowMe>")
		{
			int index = fAlias.PositionOf(alias);
			if (index == -1)
				cout << "alias \"" << alias << "\" not defined" << endl;
			else
				cout << alias << " -> \"'" << fAliasCommand[index] << '\"' << endl;
			return true;
		}
		else return MakeAlias(alias, target);
	}
	else if (command.Name() == "history")
	{
		for (int i = fHistory.Length() - 2; i > -1; i--)
			if (fHistory[i] != NULL)
				cout << setw(5) << fHistoryCount-i << ": " << *(fHistory[i]) << '\n';
		cout.flush();
		return true;
	}
	else if (command.Name() == "wait")
	{
		double delay;
		command.Argument(0).GetValue(delay);
		delay = (delay < 0.0) ? 0.0 : delay;
		
		/* time delay */
		clock_t start_time;
		start_time = clock();
		while((clock() - start_time) < delay * CLOCKS_PER_SEC) { }

		return true;
	}
	else if (command.Name() == "echo")
	{
		StringT message;
		command.Argument(0).GetValue(message);
		cout << " " << message << endl;
		return true;
	}
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

/* main event loop */
void iConsoleT::DoInteractive(void)
{
	/* open log stream */
	time_t the_time;
	time(&the_time);
	flog << "\n###################################################\n"
	     << "# open: " << ctime(&the_time)
	     << "###################################################\n";

	/* strings that act as command stacks */
	StringT line, log_line;

	/* run interactive */
	bool end = false;
	bool line_OK = true;
	bool do_log = false;
	while (!end)
	{
		/* shift */
		line.DropLeadingSpace();

		/* log if not executing from external */
		do_log = fInputStack.Length() == 0;
		
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
			else if (line[0] == ';') /* command separator */
			{
				/* remove separator */
				line.Drop(1);
			}
			else if (line[0] == '!') /* history */
			{
				/* get selection */
				line.Drop(1);
				StringT first;
				int count;
				first.FirstWord(line, count, false);
				
				/* find by name */
				if (first.StringLength() > 0)
				{
					line.Drop(count);

					/* find by index */
					if (isdigit(first[0]))
					{
#ifdef _MSC_VER
						istrstream s((char*) first);
#else
						istrstream s((const char*) first);
#endif
						int dex = -99199;
						s >> dex;
						dex = fHistoryCount - dex;
						if (dex < 0 || dex >= fHistory.Length() - 1)
						{
							cout << "out of range: " << dex << endl;
							line_OK = false;	
						}
						else
						{
							StringT* command = fHistory[dex];
							if (command != NULL)
								line.Prepend(*command, " ");
							line_OK = true;
						}
					}
					else /* find by name */
					{
						StringT* command = NULL;
						int hint_length = first.StringLength();
						for (int i = 0; !command && i < fHistory.Length(); i++)
						{
							StringT* test = fHistory[i];
							if (test != NULL) {
								int test_length = test->StringLength();
								int n = (hint_length < test_length) ? hint_length : test_length;
								if (strncmp(first, *test, n) == 0)  command = test;
							}
						}
						
						/* match */
						if (command) {
							line.Prepend(*command, " ");
							line_OK = true;
						}
						else {
							cout << "no match in history: \"" << first << '\"' << endl;
							line_OK = false;
						}
					}
				}
				else
					line_OK = false;
			}
			else if (line[0] == ':') /* scope change vs command */
			{
				/* look for separator */
				StringT command_line;
				NextCommand(line, command_line);
			
				/* take scope specifier */
				int count;
				StringT scope_line;
				scope_line.FirstWord(command_line, count, false);
				line.Drop(count);
			
				/* try scope change */
				iConsoleObjectT* scope = GetScope(*fCurrent, scope_line);
				
				/* change scope */
				if (scope != NULL)
				{
					SetScope(*scope);
				
					/* log (but not if read from external) */
					if (do_log) flog << ':' << fScope << '\n';
				}
				else
					line_OK = false;
			}
			else
			{
				/* resolve next word in line */
				StringT command_name;
				CommandScope scope = ResolveNextWord(line, command_name);
				
				/* dispatch */
				switch (scope)
				{
					case kConsoleCommand:
					{
						if (command_name == "end")
						{
							/* fetch command specification */
							const CommandSpecT* command_spec = iResolveCommand(command_name, line);
							if (!command_spec) {
								cout << "missing reply (y/n)" << endl;
								line_OK = false;
							} else {
								/* verify */
								StringT confirm;
								command_spec->Argument(0).GetValue(confirm);
								if (confirm[0] == 'y')
									end = true;
							}
						}
						else
						{
							/* look for separator */
							StringT command_line;
							NextCommand(line, command_line);
							line.Drop(command_line.StringLength());

							/* fetch command specification */
							const CommandSpecT* command_spec = iResolveCommand(command_name, command_line);

							/* execute */
							if (command_spec) 
							{
								line_OK = iDoCommand(*command_spec, command_line);
								
								/* log */
								if (line_OK && do_log) {
									command_spec->WriteCommand(flog);
									flog << '\n';
								}
							}
							else
								line_OK = false;							
						}
						break;
					}					
					case kConsoleVariable:
					{
						/* keep for log */
						if (do_log) log_line = line;

						/* operate */
						line_OK = iDoVariable(command_name, line);
						
						/* log */
						if (line_OK && do_log) {
							log_line.Drop(-line.StringLength());
							flog << command_name << " " << log_line << '\n';
						}
						break;
					}
					case kScopeCommand:
					{
						/* look for separator */
						StringT command_line;
						command_line.Root(line, ';');
						line.Drop(command_line.StringLength());

						/* fetch command specification */
						const CommandSpecT* command_spec = fCurrent->iResolveCommand(command_name, command_line);

						/* execute */
						if (command_spec) {
							line_OK = fCurrent->iDoCommand(*command_spec, command_line);
						
							/* log */
							if (line_OK && do_log) {
								command_spec->WriteCommand(flog);
								flog << '\n';
							}
						}
						else
							line_OK = false;
							
						break;
					}
					case kScopeVariable:
					{
						/* keep for log */
						if (do_log) log_line = line;

						/* operate */
						line_OK = fCurrent->iDoVariable(command_name, line);

						/* log */
						if (line_OK && do_log) {
							log_line.Drop(-line.StringLength());
							flog << command_name << " " << log_line << '\n';
						}
						break;	
					}
					case kAlias:
					{
						/* locate */
						int index = fAlias.PositionOf(command_name);
						if (index == -1)
							throw ExceptionT::kGeneralFail;
						else
							line.Prepend(fAliasCommand[index], " ");
						break;
					}
					default:
					{
						/* message */
						if (strlen(line) > 0)
							cout << "unresolved: \"" << line << '\"' << endl;
					
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
		if (fInputStack.Length() == 0 && fDanglingInput.Length() == 0)
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
			/* try read from stream */
			if (fInputStack.Length() > 0)
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
				}	
			}
			else /* grab from dangling input */
			{
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
		size_t length = 0;
		char* str = line;
		size_t max = strlen(str);
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
			/* try root */
			iConsoleObjectT* scope = &start;
			while (scope->iSuper() != NULL)
				scope = scope->iSuper();
			if (strncmp(line, scope->iName(), length) == 0)
			{
				line.Drop(length);
				line.Prepend(":root");
				return GetScope(start, line);
			}
			else /* failed to find any matches */ 
			{
				cout << "could not resolve scope: \""
				     << line << "\"" << endl;
				return NULL;
			}
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
	return ResolveCommandName(command);
}

/* pulls the first word from the line and resolves it into
* a command from the console or current scope, or returns
* kNone if the word could not be resolved */
iConsoleT::CommandScope iConsoleT::ResolveCommandName(StringT& command) const
{
	size_t word_length = strlen(command);
	if (word_length == 0) return kNone;

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
			ArrayT<const StringT*> exact_match_command(exact_match_count);
			ArrayT<CommandScope> exact_match_scope(exact_match_count);
			exact_match_count = 0;
			for (int i = 0; i < fWord.Length(); i++)
			{
				const StringT* pcommand = fWord[i];
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
			fstreamT::ClearLine(cin);
		}
	}
	return scope;
}

/* reset dictionary - scope_only sets only scope commands
* and variables */
void iConsoleT::BuildDictionary(bool scope_only)
{
	/* console symbols */
	const ArrayT<CommandSpecT*>& console_commands = iCommands();
	const ArrayT<StringT>& console_variables = iVariables();
		
	/* reset console symbols */
	if (!scope_only)
	{
		/* resize */
		fWord.Dimension(0);
		fWordScope.Dimension(0);
	
		/* add commands */
		for (int i = 0; i < console_commands.Length(); i++)
		{
			const StringT& name = console_commands[i]->Name();
			fWord.Append(&name);
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
		const ArrayT<CommandSpecT*>& scope_commands = fCurrent->iCommands();
		const ArrayT<StringT>& scope_variables = fCurrent->iVariables();
	
		/* add commands */
		for (int i = 0; i < scope_commands.Length(); i++)
		{
			fWord.Append(&(scope_commands[i]->Name()));
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

/* flush the command line and all input streams */
void iConsoleT::FlushInput(StringT& line)
{
	/* flush line */
	line.Drop(strlen(line));

	/* flush input streams */
	if (fstop_read_on_error) {
		for (int i = 0; i < fInputStack.Length(); i++)
		{
			cout << "closing stream: \"" << fInputStack[i]->filename()  << "\"" << endl;
			delete fInputStack[i];
			fInputStack[i] = NULL;
		}
		fInputStack.Dimension(0);
		fDanglingInput.Dimension(0);
	}
}

/* make an alias - returns false on fail */
bool iConsoleT::MakeAlias(const StringT& alias, StringT& line)
{				
	/* remove alias */
	if (line == "<DELETE>")
	{
		/* clear line */
		line.Drop(strlen(line));
	
		/* find alias */
		int index = fAlias.PositionOf(alias);
		if (index == -1) throw ExceptionT::kGeneralFail;
		
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
			cout << alias << " -> \"" << fAliasCommand.Last() << '\"' << endl;
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
	else if (fHistory[0] == NULL || /* history is empty */
	     (*(fHistory[0]) != line && line[0] != '!')) /* no repeats, no '!' commands */
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
		
		/* increment count */
		fHistoryCount++;
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

		/* increment count */
		fHistoryCount--;
	}
}

void iConsoleT::TopHistory(StringT& line)
{
	if (fHistory.Length() == 0)
		return;
	else if (fHistory[0] != NULL)
		line = *(fHistory[0]);

}

/* pull the next command from the line */
void iConsoleT::NextCommand(const StringT& source, StringT& next) const
{
	const char* a = source;
	int len = source.StringLength();
	bool cut = false, in_quote = false;
	int dex = -1;
	while (dex < len && !cut) {

		/* next */
		dex++;
	
		if (in_quote) {
			if (a[dex] == '"')
				in_quote = false;
		}
		else if (a[dex] == '"')
			in_quote = true;
		else if (a[dex] == ';')
			cut = true;
	}
	
	/* copy */
	if (cut)
		next.Take(source, dex);
	else if (dex > 0)
		next = source;
	else
		next.Clear();
}
