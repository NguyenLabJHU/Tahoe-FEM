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
$Id: MakeClean.ma,v 1.1 2001-09-17 18:20:38 paklein Exp $
:[font = subsection; inactive; Cclosed; preserveAspect; startGroup]
functions:
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
General::spell1: Possible spelling error: new symbol name "Root"
     is similar to existing symbol "Roots".
:[font = input; preserveAspect; endGroup]
Extension[string_String] := StringDrop[string, StringLength[Root[string]]]
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
remove make files:
:[font = input; preserveAspect]
MakeClean[] := MakeClean[Directory[]]
:[font = input; preserveAspect]
MakeClean[homedir_String] := Module[
	
	{files, dirs, makefile, currdir, fullpath},
	
	currdir = Directory[]; 
	fullpath = File /. FileInformation[homedir];
	SetDirectory[fullpath]; 
	Print[fullpath]; 
	files = FileNames[];
	dirs = Select[files, (FileType[#1] == Directory && # != "CVS" && # != "benchmark") & ]; 
    DeleteFile /@ FileNames["*.run"];
    	
    MakeClean /@ dirs; 
    SetDirectory[currdir]
    ]; 
:[font = input; preserveAspect]
DeleteIfPresent[files_List, killfile_String] := 
   If[Length[Flatten[Position[files, killfile]]] > 0, DeleteFile[killfile]]; 
:[font = input; preserveAspect; startGroup]
?*File*
:[font = print; inactive; preserveAspect; endGroup]
ContextToFilename EndOfFile         FileDate          FileNames         RenameFile
CopyFile          File              FileInformation   FileType          SetFileDate
DeleteFile        FileByteCount
:[font = input; preserveAspect; startGroup]
?File
:[font = print; inactive; preserveAspect; endGroup; endGroup; endGroup]
File is a possible value returned by FileType and related functions to indicate the type of
   a file.
:[font = subsection; inactive; preserveAspect; startGroup]
projects:
:[font = subsubsection; inactive; preserveAspect; startGroup]
benchmark:
:[font = input; preserveAspect; startGroup]
SetDirectory["uster:users (mac os 9):tahoe:home:"]
:[font = output; output; inactive; preserveAspect; endGroup]
"Uster:USERS (Mac OS 9):tahoe:home"
;[o]
Uster:USERS (Mac OS 9):tahoe:home
:[font = input; preserveAspect; startGroup]
FileNames[]
:[font = output; output; inactive; preserveAspect; endGroup]
{"admin", "aztec", "benchmark", "blue", "CBLAS", "contrib", "f2c", "gonzo", "guide.user"\
 
   , "html", "macros", "metis", "spooles", "spoolesMPI", "spoolesMT", "tahoe.HEAD", "toolbo\
 
   x"}
;[o]
{admin, aztec, benchmark, blue, CBLAS, contrib, f2c, gonzo, guide.user, html, macros, 
 
  metis, spooles, spoolesMPI, spoolesMT, tahoe.HEAD, toolbox}
:[font = input; preserveAspect; startGroup]
MakeClean["benchmark"]
:[font = print; inactive; preserveAspect]
Uster:USERS (Mac OS 9):tahoe:home:benchmark
Uster:USERS (Mac OS 9):tahoe:home:benchmark:comparator
Uster:USERS (Mac OS 9):tahoe:home:benchmark:comparator:inc
Uster:USERS (Mac OS 9):tahoe:home:benchmark:comparator:lib
Uster:USERS (Mac OS 9):tahoe:home:benchmark:comparator:obj
Uster:USERS (Mac OS 9):tahoe:home:benchmark:comparator:obj:cxx_repository
Uster:USERS (Mac OS 9):tahoe:home:benchmark:comparator:src
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.0
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.0:diffusion
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.0:geometry
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.0:meshfree
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.0:2D.elastodynamic
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.0:2D.elastostatic
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.0:3D.elastodynamic
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.0:3D.elastostatic
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.0:3D.lattice
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:geometry
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.01
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.02
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.03
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.04
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.05
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.06
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.07
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.08
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.09
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.10
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.11
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.16
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.17
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.18
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.30
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.60
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.61
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:2D:material.80
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:geometry
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.01
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.02
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.03
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.04
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.05
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.06
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.07
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.08
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.09
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.10
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.11
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.16
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.17
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.18
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.30
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.45
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.50
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.1:material.solid:3D:material.80
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.2
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.2:contact.simple
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.2:force.controller
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.2:geometry
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.2:K.field
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.2:2DContactPatch
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.2:2DHertzianContact
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.2:3DContactPatch
Uster:USERS (Mac OS 9):tahoe:home:benchmark:level.3
:[font = output; output; inactive; preserveAspect; endGroup; endGroup; endGroup]
"Uster:USERS (Mac OS 9):tahoe:home"
;[o]
Uster:USERS (Mac OS 9):tahoe:home
^*)
