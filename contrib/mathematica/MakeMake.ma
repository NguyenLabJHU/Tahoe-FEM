(*^
::[	Information =

	"This is a Mathematica Notebook file.  It contains ASCII text, and can be
	transferred by email, ftp, or other text-file transfer utility.  It should
	be read or edited using a copy of Mathematica or MathReader.  If you 
	received this as email, use your mail application or copy/paste to save 
	everything from the line containing (*^ down to the line containing ^*)
	into a plain text file.  On some systems you may have to give the file a 
	name ending with ".ma" to allow Mathematica to recognize it as a Notebook.
	The line below identifies what version of Mathematica created this file,
	but it can be opened using any other version as well.";

	FrontEndVersion = "Macintosh Mathematica Notebook Front End Version 2.2";

	MacintoshStandardFontEncoding; 
	
	fontset = title, inactive, noPageBreakBelow, nohscroll, preserveAspect, groupLikeTitle, center, M7, bold, e8,  24, "Times"; 
	fontset = subtitle, inactive, noPageBreakBelow, nohscroll, preserveAspect, groupLikeTitle, center, M7, bold, e6,  18, "Times"; 
	fontset = subsubtitle, inactive, noPageBreakBelow, nohscroll, preserveAspect, groupLikeTitle, center, M7, italic, e6,  14, "Times"; 
	fontset = section, inactive, noPageBreakBelow, nohscroll, preserveAspect, groupLikeSection, grayBox, M22, bold, a20,  18, "Times"; 
	fontset = subsection, inactive, noPageBreakBelow, nohscroll, preserveAspect, groupLikeSection, blackBox, M19, bold, a15,  14, "Times"; 
	fontset = subsubsection, inactive, noPageBreakBelow, nohscroll, preserveAspect, groupLikeSection, whiteBox, M18, bold, a12,  12, "Times"; 
	fontset = text, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  10, "Times"; 
	fontset = smalltext, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  10, "Times"; 
	fontset = input, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeInput, M42, N23, bold, L0,  10, "Courier"; 
	fontset = output, output, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23, L0,  10, "Courier"; 
	fontset = message, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23, R32768, L0,  10, "Courier"; 
	fontset = print, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23, R32768, B32768, L0,  10, "Courier"; 
	fontset = info, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23, B32768, L0,  10, "Courier"; 
	fontset = postscript, PostScript, formatAsPostScript, output, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeGraphics, M7, l34, w225, h227,  12, "Courier"; 
	fontset = name, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7, italic,  10, "Geneva"; 
	fontset = header, inactive, noKeepOnOnePage, preserveAspect, M7,  12, "Times"; 
	fontset = leftheader, inactive, L2,  12, "Times"; 
	fontset = footer, inactive, noKeepOnOnePage, preserveAspect, center, M7,  12, "Times"; 
	fontset = leftfooter, inactive, L2,  12, "Times"; 
	fontset = help, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  10, "Times"; 
	fontset = clipboard, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, "Times"; 
	fontset = completions, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, "Times"; 
	fontset = special1, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, "Times"; 
	fontset = special2, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, "Times"; 
	fontset = special3, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, "Times"; 
	fontset = special4, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, "Times"; 
	fontset = special5, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, "Times"; 
	automaticGrouping; currentKernel; 
]
:[font = subsubsection; inactive; preserveAspect]
$Id: MakeMake.ma,v 1.1 2001-09-17 18:20:38 paklein Exp $
:[font = subsection; inactive; Cclosed; preserveAspect; startGroup]
functions:
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
notes:
:[font = text; inactive; preserveAspect]
(1) $(PWD):
On wallaby2, this macro gives the current directly running the makefile commands while on the DEC, this macro returns the top level working directory.
;[s]
3:0,0;4,1;10,0;163,-1;
2:2,12,9,Times,0,10,0,0,0;1,12,9,Courier,1,10,0,0,0;
:[font = text; inactive; preserveAspect]
(2) MACROS as targets
GNU make does not like targets that are defined by MACROS passed in as "MACRO = $(MACRO)".
;[s]
3:0,0;4,1;21,0;113,-1;
2:2,12,9,Times,0,10,0,0,0;1,12,9,Courier,1,10,0,0,0;
:[font = text; inactive; preserveAspect]
(3) don't "touch":
"touch" produces a funny time stamp (DEC, ALASKA, SIBERIA) that is many minutes different from the "current" time. Confused file dependencies.
;[s]
3:0,0;3,1;18,0;162,-1;
2:2,12,9,Times,0,10,0,0,0;1,12,9,Courier,1,10,0,0,0;
:[font = text; inactive; preserveAspect; endGroup]
(4) subdirectories as targets
A number of platforms did not show reliable behavior when processing subdirectories as targets. Use for loops instead.
;[s]
5:0,0;4,1;29,0;130,1;133,0;149,-1;
2:3,12,9,Times,0,10,0,0,0;2,12,9,Courier,1,10,0,0,0;
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
basics:
:[font = input; preserveAspect; startGroup]
Root[string_String] := Module[

	{pos},
	
	pos = Flatten[StringPosition[string, "."]];
	If[Length[pos] == 0,
		string,
		StringDrop[string, Last[pos] - StringLength[string] - 1]
	]
];
:[font = message; inactive; preserveAspect; endGroup]
General::spell1: 
   Possible spelling error: new symbol name "Root" is similar to existing symbol 
    "Roots".
:[font = input; preserveAspect; endGroup]
Extension[string_String] := StringDrop[string, StringLength[Root[string]]]
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
insert make files into directory structure:
:[font = text; inactive; preserveAspect]
 "homedir" gets the top level driver makefile, while all subdirectories get compile makefiles with recursive calls to nested directories. object files are assumed to be generated from every ".cpp" file.
:[font = input; preserveAspect]
MakeMake[sourcetop_String] := MakeMake[Directory[], sourcetop]; 
:[font = input; preserveAspect]
MakeMake[homedir_String, sourcetop_String] := Module[

	{directories, src, fp}, 
	
	src = Sort[MakeMakeDir[sourcetop]]
	
	(*
	fp = OpenWrite["object_list"];
	WriteObjectList[fp, src];
	Close[fp]
	*)
	]; 
:[font = input; preserveAspect]
SourceExtensions = {".c", ".cpp", ".C", ".f", ".F"};
:[font = input; preserveAspect]
HeaderExtensions = {".h", ".inc"};
:[font = input; preserveAspect]
HasExtension[name_String, extensions_List] := 
	Length[Position[extensions, Extension[name]]] > 0
:[font = input; preserveAspect]
MakeMakeDir[dir_String] := Module[

	{olddir, hdr, src, subdirs, files, fp, subdirstr, spacenames}, 
    
    olddir = Directory[];
    SetDirectory[dir]; Print[dir];
    files = FileNames[]; 
    subdirs = Select[files, FileType[#1] == Directory & ];
    
    src = Select[Complement[files, subdirs], HasExtension[#, SourceExtensions]&];
    hdr = Select[Complement[files, subdirs], HasExtension[#, HeaderExtensions]&]; 
    
    subdirs = Select[subdirs, (# != "CVS" &&
                               # != "unused" && 
                               # != "UNUSED")& ]; 
    
    fp = OpenWrite["makefile"];
    WriteString[fp, "# $Id: MakeMake.ma,v 1.1 2001-09-17 18:20:38 paklein Exp $\n"];
    WriteString[fp, "# object directory makefile - compilation instructions\n"];
    WriteString[fp, "\ninclude $(MACRO_DIR)/suffix\n"];
	WriteString[fp, "include $(MACRO_DIR)/$(ARCH).macros\n"];
	 
    WriteObjectList[fp, src];
    WriteSourceList[fp, src];
    WriteHeaderList[fp, hdr];
 
    (* warn about spaces in directory names *)
    spacenames = Select[subdirs, (Length[Flatten[StringPosition[#," "]]] > 0)& ];
    If[Length[spacenames] > 0,
    	Print[""];
    	Print["WARNING: directory names with spaces"];
    	Map[Print["\t", #]&,spacenames];
    	Print[""]
    	];    
    
    subdirstr = (StringReplace[#1, " " -> "\\ "] & ) /@ subdirs;
    WriteSubDirList[fp, subdirstr];
  
  	WriteString[fp, "\n# instructions for subdirectories\n"];
  	WriteString[fp, "include $(MACRO_DIR)/subdir.targets\n"];
   
    WriteString[fp, "\n# dependencies\n"];
    WriteString[fp, "DEPEND = /dev/null\n"];
    WriteString[fp, "include $(DEPEND)\n"];
    Close[fp]; 
    
    src = Flatten[{src, MakeMakeDir /@ subdirs}];
    SetDirectory[olddir];
    src
    ]; 
:[font = input; preserveAspect]
TakeRoot[str_String] := StringTake[str, Last[StringPosition[str,"."]][[1]] - 1]
:[font = input; preserveAspect]
WriteObjectList[fp_OutputStream, src_List] := Module[
	
	{dotofiles, i}, 
	
	dotofiles = Map[StringJoin[TakeRoot[#], ".o"]&, src];
	            
	WriteString[fp, "\n# objects\n"];
	If[Length[dotofiles] > 0,
	
		WriteString[fp, "OBJ = \\\n"];
		Do[
			WriteString[fp, 
				StringJoin["\t", dotofiles[[i]], 
					If[i < Length[dotofiles], " \\\n", "\n"]
					]
				], {i, Length[dotofiles]}];
		WriteString[fp, "OBJ_LINK = $(OBJ)\n"];
		WriteString[fp, "DEP = $(OBJ:.o=.d)\n"];
		
		(* empty lists *)
		WriteString[fp, "#OBJ = dummy\n"];
		WriteString[fp, "#OBJ_LINK =\n"];
		WriteString[fp, "#DEP = /dev/null\n"]
				
		,(* no object files *)
		WriteString[fp, "#OBJ = \n"];
		WriteString[fp, "#OBJ_LINK = $(OBJ)\n"];
		WriteString[fp, "#DEP = $(OBJ:.o=.d)\n"];
		WriteString[fp, "OBJ = dummy\n"];
		WriteString[fp, "OBJ_LINK =\n"];
		WriteString[fp, "DEP = /dev/null\n"]
		]
	];
:[font = input; preserveAspect]
WriteSourceList[fp_OutputStream, src_List] := Module[
	
	{i}, 
            
	WriteString[fp, "\n# sources\n"];
	If[Length[src] > 0,
		WriteString[fp, "SRC = \\\n"];
		Do[
			WriteString[fp, 
				StringJoin["\t", src[[i]], 
					If[i < Length[src], " \\\n", "\n"]
					]
				], {i, Length[src]}]
				
		,(* no source files *)
		
		WriteString[fp, "#SRC = \n"]
		]
	];
:[font = input; preserveAspect]
WriteHeaderList[fp_OutputStream, hdr_List] := Module[
	
	{i}, 
            
	WriteString[fp, "\n# headers\n"];
	If[Length[hdr] > 0,
		WriteString[fp, "HDR = \\\n"];
		Do[
			WriteString[fp, 
				StringJoin["\t", hdr[[i]], 
					If[i < Length[hdr], " \\\n", "\n"]
					]
				], {i, Length[hdr]}];
				
		WriteString[fp, "HDR_LINK = $(HDR:.h=.h_link)\n"]		
				
		,(* no header files *)
		
		WriteString[fp, "#HDR = \n"];
		WriteString[fp, "#HDR_LINK = $(HDR:.h=.h_link)\n"]
		]
	];
:[font = input; preserveAspect]
WriteSubDirList[fp_OutputStream, subdir_List] := Module[
	
	{i}, 
            
	WriteString[fp, "\n# subdirectories\n"];
	If[Length[subdir] > 0,
		WriteString[fp, "SUB_DIR = \\\n"];
		Do[
			WriteString[fp, 
				StringJoin["\t", subdir[[i]], 
					If[i < Length[subdir], " \\\n", "\n"]
					]
				], {i, Length[subdir]}];
			WriteString[fp, "subdir_driver: subdir_loop\n"];
			WriteString[fp, "# SUB_DIR is empty\n"];
			WriteString[fp, "#subdir_driver: \n"]
							
		,(* no subdir files *)
		
		WriteString[fp, "SUB_DIR = \n"];
		WriteString[fp, "#subdir_driver: subdir_loop\n"];
		WriteString[fp, "# SUB_DIR is empty\n"];
		WriteString[fp, "subdir_driver: \n"]
		]
	];
:[font = input; preserveAspect]
WriteMacros[fp_OutputStream] := Module[

	{},
	
	WriteString[fp, "\nMACROS = \\\n"];
	WriteString[fp, "\t\"MAKE = $(MAKE)\" \\\n"];
	WriteString[fp, "\t\"MAKE_OPTS = $(MAKE_OPTS)\" \\\n"];
	WriteString[fp, "\t\"INC_DIR = $(INC_DIR)\" \\\n"];
	WriteString[fp, "\t\"OBJ_DIR = $(OBJ_DIR)\" \\\n"];
	WriteString[fp, "\t\"LIB_DIR = $(LIB_DIR)\" \\\n"];
	WriteString[fp, "\t\"TR = $(TR)\" \\\n"];
	WriteString[fp, "\t\"COMP_CC = $(COMP_CC)\" \\\n"];
	WriteString[fp, "\t\"COMP_C = $(COMP_C)\" \\\n"];
	WriteString[fp, "\t\"COMP_F = $(COMP_F)\" \\\n"];
	WriteString[fp, "\t\"CFLAGS_C = $(CFLAGS_C)\" \\\n"];
	WriteString[fp, "\t\"CFLAGS_CC = $(CFLAGS_CC)\" \\\n"];
	WriteString[fp, "\t\"CFLAGS_F = $(CFLAGS_F)\" \\\n"];
	WriteString[fp, "\t\"LIBRARY = $(LIBRARY)\" \\\n"];
	WriteString[fp, "\t\"AR = $(AR)\" \\\n"];
	WriteString[fp, "\t\"ARFLAGS = $(ARFLAGS)\" \\\n"];
	WriteString[fp, "\t\"RM = $(RM)\" \\\n"];
	WriteString[fp, "\t\"RM_FILES = $(RM_FILES)\" \\\n"];
	WriteString[fp, "\t\"LN = $(LN)\" \\\n"];
	WriteString[fp, "\t\"MAKE_DEPEND = $(MAKE_DEPEND)\" \\\n"];
	WriteString[fp, "\t\"DEPEND_PATH = $(DEPEND_PATH)\" \\\n"];
	WriteString[fp, "\t\"DEPEND = $(DEPEND)\" \\\n"];
	WriteString[fp, "\t\"ECHO = $(ECHO)\" \\\n"];
	WriteString[fp, "\t\"CAT = $(CAT)\"\n"];
	
	WriteString[fp, "\n# debugging\n"];
	WriteString[fp, "echo_all:\n"];
	WriteString[fp, "\t@ $(ECHO) \"MAKE = $(MAKE)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"MAKE_OPTS = $(MAKE_OPTS)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"INC_DIR = $(INC_DIR)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"OBJ_DIR = $(OBJ_DIR)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"LIB_DIR = $(LIB_DIR)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"TR = $(TR)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"COMP_C = $(COMP_C)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"COMP_CC = $(COMP_CC)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"COMP_F = $(COMP_F)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"CFLAGS_C = $(CFLAGS_C)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"CFLAGS_CC = $(CFLAGS_CC)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"CFLAGS_F = $(CFLAGS_F)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"LIBRARY = $(LIBRARY)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"AR = $(AR)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"ARFLAGS = $(ARFLAGS)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"RM = $(RM)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"RM_FILES = $(RM_FILES)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"LN = $(LN)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"MAKE_DEPEND = $(MAKE_DEPEND)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"DEPEND_PATH = $(DEPEND_PATH)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"DEPEND = $(DEPEND)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"ECHO = $(ECHO)\"\n"];
	WriteString[fp, "\t@ $(ECHO) \"CAT = $(CAT)\"\n"]
	];
:[font = input; preserveAspect]
WriteClean[fp_OutputStream, dirs_List] := Module[

	{i},
	
	WriteString[fp, "\nclean:\n"];
	WriteString[fp, "\t-@ $(RM) $(RM_FILES)\n"];
	WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = clean\"\n"]
	
	(*
	If[Length[dirs] > 0,
		WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = clean\"\n"],
		WriteString[fp, "#\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = clean\"\n"]
		]
	*)

	(*
	WriteString[fp, "\nclean_all:\n"];
	WriteString[fp, "\t-@ $(CLEAN_RM) $(CLEAN_FILES)\n"];
	If[Length[dirs] > 0,
		WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = clean_all\" \\\n"];
		WriteString[fp, "\t\"CLEAN_RM = $(CLEAN_RM)\" \\\n"];
		WriteString[fp, "\t\"CLEAN_FILES = $(CLEAN_FILES)\"\n"]
		]
	*)
	]; 
:[font = input; preserveAspect]
WriteDepend[fp_OutputStream, src_List, dirs_List] := Module[

	{i}, 
	
	WriteString[fp, "\ndepend_init: all.depend\n"];
	WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = depend_init\"\n"];
	
	(*
	If[Length[dirs] > 0,
		WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = depend_init\"\n"],
		WriteString[fp, "#\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = depend_init\"\n"]
		];
	*)
		
	(* local depend file *)
	WriteString[fp, "\nall.depend: $(DEP)\n"];
	WriteString[fp, "\t@ $(ECHO) \"    $(CAT): $(DEP) > all.depend\"; $(CAT) $(DEP) > all.depend\n"]
	]; 
:[font = input; preserveAspect]
WriteHeaders[fp_OutputStream, hdr_List, dirs_List] := Module[

	{i}, 
	
	WriteString[fp, "\nheaders: $(HDR_LINK)\n"];
	WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = headers\"\n"]
	
	(*
    If[Length[dirs] > 0,
		WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = headers\"\n"],
		WriteString[fp, "#\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = headers\"\n"]
		]
	*)
	]; 
:[font = input; preserveAspect]
WriteObjects[fp_OutputStream, src_List, dirs_List] := Module[

	{}, 
	
    WriteString[fp, "\nobjects: all.depend\n"];
    WriteString[fp, "\t@ $(MAKE) $(MACROS) update_objects\n"];
    WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = objects\"\n"];
    WriteString[fp, "\nupdate_objects: $(OBJ)\n"]
    
    (*
    WriteString[fp, "\nobjects:", If[Length[src] > 0," all.depend $(OBJ)\n","\n"]];
    WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = objects\"\n"]
   	*)
    
    (*
    If[Length[dirs] > 0,
		WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = objects\"\n"],
		WriteString[fp, "#\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = objects\"\n"]
		]
	*)
	]; 
:[font = text; inactive; preserveAspect]
This leads to an order N^2 time to update the library since each call requires a search:
:[font = input; inactive; preserveAspect]
WriteArchive[fp_OutputStream, src_List, dirs_List] := Module[

	{}, 
	
    WriteString[fp, "\nlibrary:\n"];
    WriteString[fp, "\t@ $(AR) $(ARFLAGS) $(LIBRARY) $(OBJ_LINK)\n"];
    WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = library\"\n"]
    
    (*
    WriteString[fp, "library_local:\n"];
    If[Length[src] > 0,
    	WriteString[fp, "\t@ $(AR) $(ARFLAGS) $(LIBRARY) $(OBJ)\n"],
    	WriteString[fp, "#\t@ $(AR) $(ARFLAGS) $(LIBRARY) $(OBJ)\n"]
    	];
    WriteString[fp, "library_subdir:\n"];
    If[Length[dirs] > 0,
		WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = library\"\n"],
		WriteString[fp, "#\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = library\"\n"]
		]
	*)
	]; 
:[font = input; preserveAspect]
WriteArchive[fp_OutputStream, src_List, dirs_List] := Module[

	{}, 
	
    WriteString[fp, "\nlibrary:\n"];
    WriteString[fp, "\t@ $(ECHO) \"\\\\\" >> $(LIBRARY)\n"];
    WriteString[fp, "\t@ $(ECHO) \"\\t$(OBJ_LINK)\\c\" >> $(LIBRARY)\n"];
    WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = library\"\n"]
    
    (*
    WriteString[fp, "library_local:\n"];
    If[Length[src] > 0,
    	WriteString[fp, "\t@ $(AR) $(ARFLAGS) $(LIBRARY) $(OBJ)\n"],
    	WriteString[fp, "#\t@ $(AR) $(ARFLAGS) $(LIBRARY) $(OBJ)\n"]
    	];
    WriteString[fp, "library_subdir:\n"];
    If[Length[dirs] > 0,
		WriteString[fp, "\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = library\"\n"],
		WriteString[fp, "#\t@ $(MAKE) $(MACROS) subdir_driver \"COMMAND = library\"\n"]
		]
	*)
	]; 
:[font = input; preserveAspect]
WriteGroupDrivers[fp_OutputStream, hdr_List, subdir_List] := Module[

	{},

	WriteString[fp, "\n$(SUB_DIR): dummy\n"];
	WriteString[fp, "\t@ $(ECHO) \" cd $@\"; \\\n"];
	WriteString[fp, "\tcd $@; $(MAKE) $(COMMAND) $(MACROS) \"CURR_DIR = $(CURR_DIR)/$@\"\n"]
	
	WriteString[fp, "\ndummy:\n"]
	];
:[font = input; inactive; preserveAspect]
WriteGroupDrivers[fp_OutputStream, hdr_List, subdir_List] := Module[

	{},
	
	If[Length[subdir] > 0,
		WriteString[fp, "\n$(SUB_DIR): force\n"];
		WriteString[fp, "\t@ $(ECHO) \" cd $@\"; \\\n"];
		WriteString[fp, "\tcd $@; $(MAKE) $(COMMAND) $(MACROS) \"CURR_DIR = $(CURR_DIR)/$@\"\n"]
		];
	WriteString[fp, "\nforce:\n"]
	];
:[font = input; preserveAspect; endGroup]
WriteSuffixRules[fp_OutputStream] := Module[

	{}, 
	
	WriteString[fp, "\n# suffix rules\n"];
	
	WriteString[fp, "\n.c.o:\n"];
	WriteString[fp, "\t@ $(ECHO) \"    $(COMP_C): \"$<; \\\n"];
	WriteString[fp, "\t$(COMP_C) $(CFLAGS_C) $<\n"];
	WriteString[fp, "\t@$(LN) $(CURR_DIR)/$@ $(OBJ_DIR)/$@\n"];
	WriteString[fp, "\t@$(ECHO) \"#touch\" > $(LIB_DIR)/lib.junk\n"];
	
	WriteString[fp, "\n.cpp.o:\n"];
	WriteString[fp, "\t@ $(ECHO) \"    $(COMP_CC): \"$<; \\\n"];
	WriteString[fp, "\t$(COMP_CC) $(CFLAGS_CC) $<\n"];
	WriteString[fp, "\t@$(LN) $(CURR_DIR)/$@ $(OBJ_DIR)/$@\n"];
	WriteString[fp, "\t@$(ECHO) \"#touch\" > $(LIB_DIR)/lib.junk\n"];	
	
	WriteString[fp, "\n.C.o:\n"];
	WriteString[fp, "\t@ $(ECHO) \"    $(COMP_CC): \"$<; \\\n"];
	WriteString[fp, "\t$(COMP_CC) $(CFLAGS_CC) $<\n"];
	WriteString[fp, "\t@$(LN) $(CURR_DIR)/$@ $(OBJ_DIR)/$@\n"];
	WriteString[fp, "\t@$(ECHO) \"#touch\" > $(LIB_DIR)/lib.junk\n"];

	WriteString[fp, "\n.f.o:\n"];
	WriteString[fp, "\t@ $(ECHO) \"    $(COMP_F): \"$<; \\\n"];
	WriteString[fp, "\t$(COMP_F) $(CFLAGS_F) $<\n"];
	WriteString[fp, "\t@$(LN) $(CURR_DIR)/$@ $(OBJ_DIR)/$@\n"];
	WriteString[fp, "\t@$(ECHO) \"#touch\" > $(LIB_DIR)/lib.junk\n"];
	
	WriteString[fp, "\n.F.o:\n"];
	WriteString[fp, "\t@ $(ECHO) \"    $(COMP_F): \"$<; \\\n"];
	WriteString[fp, "\t$(COMP_F) $(CFLAGS_F) $<\n"];
	WriteString[fp, "\t@$(LN) $(CURR_DIR)/$@ $(OBJ_DIR)/$@\n"];
	WriteString[fp, "\t@$(ECHO) \"#touch\" > $(LIB_DIR)/lib.junk\n"];
	
	WriteString[fp, "\n.h.h_link:\n"];
	WriteString[fp, "\t@ $(ECHO) \"    $(LN): \"$<; \\\n"];
	WriteString[fp, "\t$(LN) $(CURR_DIR)/$< $(INC_DIR)/$<\n"];
	
	WriteString[fp, "\n.c.d:\n"];
	WriteString[fp, "\t@ $(ECHO) \"    $(MAKE_DEPEND): \"$<; \\\n"];
    WriteString[fp, "\t$(ECHO) \"\" > $@; \\\n"];
    WriteString[fp, "\t$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<\n"];
	
	WriteString[fp, "\n.cpp.d:\n"];
	WriteString[fp, "\t@ $(ECHO) \"    $(MAKE_DEPEND): \"$<; \\\n"];
    WriteString[fp, "\t$(ECHO) \"# DO NOT DELETE\" > $@; \\\n"];
    WriteString[fp, "\t$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<\n"];

	WriteString[fp, "\n.C.d:\n"];
	WriteString[fp, "\t@ $(ECHO) \"    $(MAKE_DEPEND): \"$<; \\\n"];
    WriteString[fp, "\t$(ECHO) \"# DO NOT DELETE\" > $@; \\\n"];
    WriteString[fp, "\t$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<\n"];
    
    WriteString[fp, "\n.f.d:\n"];
	WriteString[fp, "\t@ $(ECHO) \"    $(MAKE_DEPEND): \"$<; \\\n"];
    WriteString[fp, "\t$(ECHO) \"# DO NOT DELETE\" > $@; \\\n"];
    WriteString[fp, "\t$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<\n"];
    
    WriteString[fp, "\n.F.d:\n"];
	WriteString[fp, "\t@ $(ECHO) \"    $(MAKE_DEPEND): \"$<; \\\n"];
    WriteString[fp, "\t$(ECHO) \"# DO NOT DELETE\" > $@; \\\n"];
    WriteString[fp, "\t$(MAKE_DEPEND) $(DEPEND_PATH) -f$@ $<\n"]
	];
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
remove make files:
:[font = input; preserveAspect]
UnMakeMake[] := UnMakeMake[Directory[]]
:[font = input; preserveAspect]
UnMakeMake[homedir_String] := Module[
	
	{files, dirs, makefile, currdir, fullpath},
	
	currdir = Directory[]; 
	fullpath = File /. FileInformation[homedir];
	SetDirectory[fullpath]; 
	Print[fullpath]; 
	files = FileNames[];
	dirs = Select[files, (FileType[#1] == Directory && # != "CVS") & ]; 
    (DeleteIfPresent[files, #1] & ) /@ {"Icon\n", "makefile", "makefile~", "makefile.bak"}; 
    If[Length[Flatten[Position[dirs, "ii_files"]]] > 0,
    	Print[SetDirectory["ii_files"]]; 
    	DeleteFile["*.ii"];
    	SetDirectory[fullpath]; 
    	DeleteDirectory["ii_files"]
    	];
    	
    UnMakeMake /@ dirs; 
    SetDirectory[currdir]
    ]; 
:[font = input; preserveAspect; endGroup; endGroup]
DeleteIfPresent[files_List, killfile_String] := 
   If[Length[Flatten[Position[files, killfile]]] > 0, DeleteFile[killfile]]; 
:[font = subsection; inactive; preserveAspect; startGroup]
projects:
:[font = subsubsection; inactive; preserveAspect; startGroup]
module:
:[font = input; preserveAspect; startGroup]
SetDirectory["uster:users (mac os 9):tahoe:home"]
:[font = output; output; inactive; preserveAspect; endGroup]
"Uster:USERS (Mac OS 9):tahoe:home"
;[o]
Uster:USERS (Mac OS 9):tahoe:home
:[font = input; preserveAspect; startGroup]
TableForm[dirs = FileNames[], TableSpacing -> {0,3}]
:[font = output; output; inactive; preserveAspect; endGroup]
TableForm[{"aztec", "benchmark", "blue", "CBLAS", ".DS_Store", "f2c", "gonzo", "gui\
 
    de.user", "html", "macros", "metis", "spooles", "spoolesMPI", "spoolesMT", "tah\
 
    oe.HEAD", "toolbox"}, TableSpacing -> {0, 3}]
;[o]
aztec
benchmark
blue
CBLAS
.DS_Store
f2c
gonzo
guide.user
html
macros
metis
spooles
spoolesMPI
spoolesMT
tahoe.HEAD
toolbox
:[font = input; preserveAspect]
dirs = Complement[dirs, {"macros"}]
:[font = text; inactive; preserveAspect]
one:
:[font = input; preserveAspect; startGroup]
Map[(SetDirectory[#]; 
	UnMakeMake["src"];
	MakeMake["src"];
	SetDirectory["::"])&, {"gonzo"}];
:[font = print; inactive; preserveAspect; endGroup]
Uster:USERS (Mac OS 9):tahoe:home:gonzo:src
src
:[font = text; inactive; preserveAspect]
process all of the modules:
:[font = input; preserveAspect; endGroup]
Map[(SetDirectory[#]; 
	UnMakeMake["src"];
	MakeMake["src"];
	SetDirectory["::"])&, dirs];
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
blue:
:[font = input; preserveAspect; startGroup]
SetDirectory["uster:users (mac os 9):tahoe:home"]
:[font = output; output; inactive; preserveAspect; endGroup]
"Uster:USERS (Mac OS 9):tahoe:home"
;[o]
Uster:USERS (Mac OS 9):tahoe:home
:[font = input; preserveAspect; startGroup]
TableForm[dirs = FileNames[], TableSpacing -> {0,3}]
:[font = output; output; inactive; preserveAspect; endGroup]
TableForm[{"aztec", "benchmark", "blue", "CBLAS", ".DS_Store", "f2c", "guide.user"\
 
    , "html", "macros", "spooles", "spoolesMPI", "tahoe", "tahoe.stable", "toolbox"}
 
   , TableSpacing -> {0, 3}]
;[o]
aztec
benchmark
blue
CBLAS
.DS_Store
f2c
guide.user
html
macros
spooles
spoolesMPI
tahoe
tahoe.stable
toolbox
:[font = input; preserveAspect; startGroup]
dirs = Complement[dirs, {"macros"}]
:[font = output; output; inactive; preserveAspect; endGroup]
{"aztec", "benchmark", "blue", "CBLAS", ".DS_Store", "f2c", "guide.user", "html"\
 
   , "spooles", "spoolesMPI", "tahoe", "tahoe.stable", "toolbox"}
;[o]
{aztec, benchmark, blue, CBLAS, .DS_Store, f2c, guide.user, html, spooles, 
 
  spoolesMPI, tahoe, tahoe.stable, toolbox}
:[font = text; inactive; preserveAspect]
one:
:[font = input; preserveAspect; startGroup]
Map[(SetDirectory[#]; 
	UnMakeMake["src"];
	MakeMake["src"];
	SetDirectory["::"])&, {"blue"}];
:[font = print; inactive; preserveAspect; endGroup]
Uster:USERS (Mac OS 9):tahoe:home:blue:src
Uster:USERS (Mac OS 9):tahoe:home:blue:src:dynamics
Uster:USERS (Mac OS 9):tahoe:home:blue:src:energy
Uster:USERS (Mac OS 9):tahoe:home:blue:src:lattice
Uster:USERS (Mac OS 9):tahoe:home:blue:src:main
src
dynamics
energy
lattice
main
:[font = text; inactive; preserveAspect]
process all of the modules:
:[font = input; preserveAspect; endGroup; endGroup]
Map[(SetDirectory[#]; 
	UnMakeMake["src"];
	MakeMake["src"];
	SetDirectory["::"])&, dirs];
^*)
