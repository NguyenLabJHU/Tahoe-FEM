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
	fontset = input, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeInput, M42, N23, bold, L-5,  10, "Courier"; 
	fontset = output, output, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23, L-5,  10, "Courier"; 
	fontset = message, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23, R65535, L-5,  10, "Courier"; 
	fontset = print, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23, R65535, B65535, L-5,  10, "Courier"; 
	fontset = info, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23, B65535, L-5,  10, "Courier"; 
	fontset = postscript, PostScript, formatAsPostScript, output, inactive, noPageBreakInGroup, nowordwrap, thinLines, preserveAspect, groupLikeGraphics, M7, l34, w350, h356,  12, "Courier"; 
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
	paletteColors = 128; automaticGrouping; currentKernel; 
]
:[font = subsubsection; inactive; preserveAspect]
 $Id: Z.ma,v 1.1 2001-09-17 18:20:39 paklein Exp $
:[font = subsection; inactive; Cclosed; preserveAspect; startGroup]
Documentation:
:[font = text; inactive; preserveAspect]
How the Z's work:

The strategy is to identify all subexpressions with a Depth less than 2 (ParseCommon).  Some of these subexpressions will have a Length greater than 2, which means there may be some common sub-subexpressions (CommonTermsOne, CommonTermsAll).  z's are chosen from the list of subexpressions with Depth less than 2 and greater than 1 (no atoms).  The expression is then evaluated with the sub-expressions replaced with the z's and the function (ToZ) is recursively applied until the final expression appears the same as its ParseCommon result.
;[s]
21:0,0;73,1;78,0;92,1;103,0;148,1;154,0;228,1;242,0;244,1;258,0;262,1;263,0;314,1;319,0;440,1;441,0;462,1;465,0;541,1;552,0;561,-1;
2:11,12,9,Times,0,10,0,0,0;10,12,9,Courier,1,10,0,0,0;
:[font = text; inactive; preserveAspect; endGroup]
Dealing with If statements:

There are (at least) 2 issues that arise when tryting to apply the z's to If statements.

(1) If has the HoldAll attribute
(2) The conditional expression in the If statement should not be z'd

(1) If must be HoldAll because the there my by improper expressions in the branches of the statement that should not be evaluated if the condition does not apply, ie:

If[Length[a] >= 2, a[[2]], 0]

If a is not a List then part 2 of a will not exist and cannot be evaluated.  The difficulty arises with the z's when trying to eliminate common subexpressions from the branches of the If statements, as well as from the conditional expression.  The sub-expressions returned by ParseCommon will not in general be the same as the unevaluated forms.  The solution is to cause all arguments of the If statements to be evaluated (just once for now).   Evaluation occurs in ToZ which then passes the result onto ToZShell for the sub-expression elimination.  It may be necessary to re-evaluate all the arguments in the If statements during the z-ing process to allow continued elimintation of subexpressions, but this question is still under investigation.  If the routines enter an infinite loop when trying to z an expression containing complicated If statements, the likely solution will be to force evaluation of the arguments again after some sub-expressions have been eliminated.

(2) In general, the conditional expression in an If statements would be z's like any other.  For C/C++ and "if" statement is false if, and only if, the the condition returns integer(0), which may be difficult or unreliable if cast from a double precision z.  Therefore, the conditional statement should not be z'd.  This filtering of conditional statements is achieved by eliminating conditional statements with a Depth of 2 (that is, where each argument is a simple expression) from the sub-expressions returned by ParseCommon, and allowing entire If statements to be treated as a single sub-expression if each of the arguments has a Depth of 2 or less.
;[s]
63:0,0;13,1;15,0;96,1;97,0;103,1;105,0;123,1;125,0;134,1;141,0;190,1;192,0;217,1;218,0;226,1;228,0;237,1;244,0;390,1;421,0;424,1;425,0;435,1;439,0;455,1;456,0;529,1;530,0;605,1;607,0;697,1;708,0;814,1;816,0;888,1;891,0;926,1;934,0;1032,1;1034,0;1057,1;1058,0;1225,1;1226,0;1264,1;1266,0;1449,1;1451,0;1472,1;1473,0;1655,1;1656,0;1710,1;1711,0;1814,1;1819,0;1916,1;1927,0;1949,1;1951,0;2035,1;2040,0;2055,-1;
2:32,12,9,Times,0,10,0,0,0;31,12,9,Courier,1,10,0,0,0;
:[font = subsection; inactive; Cclosed; preserveAspect; startGroup]
To Z:
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
Recursive parsing functions: 
:[font = text; inactive; preserveAspect]
Find common terms with the same Head in the given List:
;[s]
5:0,0;32,1;36,0;50,1;54,0;56,-1;
2:3,12,9,Times,0,10,0,0,0;2,12,9,Courier,1,10,0,0,0;
:[font = input; preserveAspect]
CommonTermsAll[list_,head_] := Module[ 

	{temp1},
		
	If[ SameQ[Head[list],List],
		temp1 = Select[list, (Depth[#] > 1 && Head[#] == head) &];
		temp1 = Flatten[Outer[Intersection,temp1,temp1]];
		Complement[temp1,list],
		{}]
	];
:[font = input; preserveAspect]
CommonTermsOne[list_,head_] := Module[ 

	{temp2, temp1, n, terms, n1, n2},
		
	If[ SameQ[Head[list],List],
		temp1 = Select[list, (Depth[#] > 1 && Head[#] == head) &];
		n = Length[temp1];
		
		terms = {};
		
		For[ n1 = 1, n1 <= n, n1++,
			For[ n2 = n1 + 1, n2 <= n, n2++,
				temp2 = Complement[{Intersection[temp1[[n1]],temp1[[n2]]]},list];
				terms = Join[terms,temp2];
				
				If[ Length[temp2] > 0 && 
					Depth[temp2[[1]]] < 3 && Length[temp2[[1]]] < 3,
					n1 = n2 = n] 
				]
			];
		
		terms ,
		{}]
	];
:[font = text; inactive; preserveAspect]
Remove elements from list after i that contain list[[i]]:
;[s]
7:0,0;21,1;25,0;32,1;33,0;47,1;56,0;58,-1;
2:4,12,9,Times,0,10,0,0,0;3,12,9,Courier,1,10,0,0,0;
:[font = text; inactive; preserveAspect]
Parse driver.  Roughly 10x faster than before.
:[font = input; preserveAspect]
Parse[expr_] := If[Depth[expr] > 2,
	Union[Flatten[Map[ If[Depth[#] > 2, Parse[#], #]&,  Apply[List,expr] ]]],
		{expr}];
:[font = text; inactive; preserveAspect]
special If handling:
;[s]
3:0,0;8,1;10,0;21,-1;
2:2,12,9,Times,0,10,0,0,0;1,12,9,Courier,1,10,0,0,0;
:[font = input; preserveAspect]
Parse[If[x_,y__]] := If[ Max[Map[Depth,{x,y}]] > 2,
								Parse[{x,y}],
								{If[x,y]}];
:[font = input; preserveAspect]
Parse[{}] := {};
:[font = text; inactive; preserveAspect]
Old version that's If-compliant.
;[s]
3:0,0;19,1;21,0;33,-1;
2:2,12,9,Times,0,10,0,0,0;1,12,9,Courier,1,10,0,0,0;
:[font = input; inactive; preserveAspect]
Parse[expr_,res_:{}] := Module[ 

	{n, i, j , res2, restemp, lres, temp1, temp2},

	res2 = res;

	If[ !SameQ[Head[expr], If],

		n = Length[expr];
	
		If[Depth[expr] > 2,
			
			Do[ res2 = Parse[expr[[i]],res2],
				{i,n}],

			res2 = Union[res2,Flatten[{expr}]]
		
			]
	,
	(* else *)

	(* special treatment for "if's" *)
		If[ Depth[ expr[[1]] ] > 2 ||
			Depth[ expr[[2]] ] > 2 ||
			Depth[ expr[[3]] ] > 2 ,

			Do[ res2 = Parse[expr[[i]],res2],
				{i,3}],

			res2 = Union[res2,Flatten[{expr}]]		
		
		]
	];
	
	Select[Sort[res2], (!SkipZQ[Head[#]])&]
	];
:[font = input; Cclosed; preserveAspect; startGroup]
SkipZQ[_] = False;
SkipZFunctions = {Less,Greater,Equal,GreaterEqual,LessEqual};
Unprotect[Evaluate[SkipZFunctions]];
Map[(# /: SkipZQ[#] = True)&, SkipZFunctions]
Protect[Evaluate[SkipZFunctions]]
:[font = output; output; inactive; preserveAspect]
{True, True, True, True, True}
;[o]
{True, True, True, True, True}
:[font = output; output; inactive; preserveAspect; endGroup]
{"Less", "Greater", "Equal", "GreaterEqual", "LessEqua\
 
   l"}
;[o]
{Less, Greater, Equal, GreaterEqual, LessEqual}
:[font = input; preserveAspect]
ParseCommon[expr_,opt_:All] := Module[ {res1},

	res1 = Select[Parse[expr], (!SkipZQ[Head[#]])&];
	
	If[ SameQ[opt, All],
	
		res1 = Union[res1, CommonTermsAll[res1, Times], 
	                       CommonTermsAll[res1, Plus]],
	           
	 	res1 = Union[res1, CommonTermsOne[res1, Times], 
	                       CommonTermsOne[res1, Plus]]
		];
		
	(* Print["parse = ", res1]; *)
	res1	
	
	];
:[font = text; inactive; preserveAspect]
better handling of expressions involving rationals:
:[font = input; preserveAspect]
ParseRational[Times[Rational[n_,d_], x_?AtomQ]] := Times[Rational[n,d],x]
:[font = input; preserveAspect]
ParseRational[Times[Rational[n_,d_], x_]] := Parse[x]
:[font = input; preserveAspect]
ParseRational[x_] := x
:[font = text; inactive; preserveAspect]
Used by ZCount:
:[font = input; preserveAspect; endGroup]
ParseAll[expr_,res_:{}] := Module[ 

	{n, i, j , res2, restemp, lres, temp1, temp2},

	res2 = res;
	n = Length[expr];

	If[Depth[expr] > 1,
		
		Do[ res2 = ParseAll[expr[[i]],res2],
			{i,n}],

		res2 = Union[res2,{expr}]
		];
	
	Sort[Flatten[{res2}]]
	];

:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
ToZ Utilities:
:[font = input; preserveAspect]
CheckName[z_:z,n_] := Module[ {i, expr, zstr, istr, zname, none},

	zstr = ToString[z];
	none = True;	

	Do[
		istr = ToString[i];
		zname = StringJoin[zstr,istr];
		expr = ToExpression[zname];
		If[ !SameQ[zname,ToString[expr]],
			none = False;
			Print[zname," = ",expr]]
		, {i,n}];
		
	If[ none, Print["OK"]]
		];
:[font = input; preserveAspect]
ClearName[z_:z,n_] := Module[ {i, zstr, istr},

	zstr = ToString[z];
	
	Do[
		istr = ToString[i];
		Clear[Evaluate[StringJoin[zstr,istr]]]
		, {i,n}]
		];
:[font = input; preserveAspect]
PrintZ[z_,n_,num_:12] := Module[ {i,str,zstr,prnt},

	If[n > 0,
	zstr = ToString[z];
	str = "{";
	
	Do[
		str = StringJoin[str, zstr, ToString[i],", "];
		If[ Mod[i,num] == 0 && i != n,
			Print[str];
			str = ""
			]
		
		,{i,n}];
		
	str = StringJoin[StringTake[str,
						Max[StringLength[str]-2,0]],"}"];

	Print[str]
	]
];
:[font = text; inactive; preserveAspect]
To skip 1 and 0:
:[font = input; preserveAspect; endGroup]
ZNumberQ[x_] := (Head[x] == Real) && (x != 0) && (x != 1);
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
ToZ Drivers:
:[font = input; preserveAspect]
ZLengthLimit = 10;
:[font = text; inactive; preserveAspect]
rules to clean up expressions with 
:[font = input; preserveAspect]
NegativeQ[n_] := n < 0
:[font = input; preserveAspect]
HoldPower[z_, n_, d_] := HoldForm[Power[Power[z, Rational[n,d]],-1]]
:[font = input; preserveAspect]
HoldProduct[z_, n_, d_] := HoldForm[Times[-1, Times[Rational[n,d], z]]]
:[font = input; preserveAspect]
RationalPowerRule = Power[z_, Rational[n_?NegativeQ,d_]] :> HoldPower[z, -n, d]
:[font = input; preserveAspect]
RationalProductRule = Times[Rational[n_?NegativeQ,d_], z_] :> HoldProduct[z, -n, d]
:[font = text; inactive; preserveAspect]
driver:
:[font = input; preserveAspect]
ToZ[rawexpr_, z_Symbol, nz_Integer, opts___] := Module[ 

	{nz2,parse,copy,rules,n,i,zstr, znstr,temp,newrule,Parseopt,
	 Defaults, expr, nextparse},
	
	Ifattributes = Attributes[If];
	Map[ClearAttributes[If, #]&, Ifattributes];
	expr = rawexpr /. If[_,0,0] -> 0;
	Map[SetAttributes[If, #]&, Ifattributes];

	(* Default options *)
	Defaults = {Parse -> One};
	
	Parseopt = Parse /. {opts} /. Defaults;
	nz2 = nz;
	rules = {};
	zstr = ToString[z];

	parse = ParseCommon[expr,Parseopt];

(*	Print["expr  = ", expr];
	Print["parse = ", parse];
	Print["fl    = ", Union[Flatten[{expr}]]]; *)

	If[ !SameQ[parse,Union[Flatten[{expr}]]] || Max[Map[Depth, parse]] > 1,
(*	If[ !SameQ[parse,RemoveRepeats[expr]],	 *)
(*	If[ !SameQ[parse,Sort[Flatten[{expr}]]], *)
	
		(* parse = Join[Select[parse,ZNumberQ],Select[parse, Depth[#] > 1 &]]; *)
		parse = Select[parse, Depth[#] > 1 &];

		n = Length[parse];	
		While[ n > 0,
		
				(* didn't make much difference *)
				(* parse = Sort[parse, (Length[#1] < Length[#2])&]; *)
			
				parse = Sort[parse];
				
				nz2++;
				znstr = StringJoin[zstr,ToString[nz2]];
				
				nextparse = parse[[1]];
				
				(* limit length of expressions *)
				If[ (Head[nextparse] === Times || Head[nextparse] === Plus) && 
					Length[nextparse] > ZLengthLimit,
					nextparse = Take[nextparse, ZLengthLimit]]; 
				
				newrule = nextparse -> ToExpression[znstr];
				rules = Append[rules, newrule];
				Print[zstr, nz2," = ", InputForm[nextparse],";"];
				(* PrintRule[newrule]; *)
			
				parse = parse /. newrule;
			
				(* find new commmon sub-expressions *)
				(*
				If[SameQ[Parseopt,All],
				parse = Union[parse,ParseCommon[Select[parse,(Length[#] > 2)&],Parseopt]];
					];
				*)
				
				(* parse = Join[Select[parse,ZNumberQ],Select[parse, Depth[#] > 1 &]]; *)
				parse = Select[parse, Depth[#] > 1 &];
								
				n = Length[parse]
				];
		
		
		(* recursion *)
		temp  = ToZ[expr //. rules, z, nz2, opts];
		copy  = temp[[1]];
		rules = Join[rules, temp[[2]]];
		nz2   = temp[[3]],
			
		Print[copy = expr];
		PrintZ[z, nz2]
		];
		
	{copy, rules, nz2}
	
	];
:[font = input; preserveAspect; endGroup]
ToZ[expr_,z_Symbol] := ToZ[expr,z,0];
ToZ[expr_] := ToZ[expr,z];
ToZ[expr_,opts___] := ToZ[expr,z,0,opts];
ToZ[expr_,z_Symbol,opts___] := ToZ[expr,z,0,opts];
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
new testing:
:[font = input; preserveAspect]
ToZTEST[rawexpr_, z_Symbol, nz_Integer, opts___] := Module[ 

	{nz2,parse,copy,rules,n,i,zstr, znstr,temp,newrule,Parseopt,
	 Defaults, expr, nextparse},
	
	Ifattributes = Attributes[If];
	Map[ClearAttributes[If, #]&, Ifattributes];
	expr = rawexpr /. {If[_,0,0] -> 0, RationalPowerRule, RationalProductRule};
	Map[SetAttributes[If, #]&, Ifattributes];

	(* Default options *)
	Defaults = {Parse -> One};
	
	Parseopt = Parse /. {opts} /. Defaults;
	nz2 = nz;
	rules = {};
	zstr = ToString[z];

	parse = ParseCommon[expr,Parseopt];

	If[ !SameQ[parse,Union[Flatten[{expr}]]],
(*	If[ !SameQ[parse,RemoveRepeats[expr]],	 *)
(*	If[ !SameQ[parse,Sort[Flatten[{expr}]]], *)
	
		(* parse = Join[Select[parse,ZNumberQ],Select[parse, Depth[#] > 1 &]]; *)
		parse = Select[parse, Depth[#] > 1 &];

		n = Length[parse];	
		While[ n > 0,
		
				(* didn't make much difference *)
				(* parse = Sort[parse, (Length[#1] < Length[#2])&]; *)
			
				parse = Sort[parse];
				
				nz2++;
				znstr = StringJoin[zstr,ToString[nz2]];
				
				nextparse = parse[[1]];
				
				(* limit length of expressions *)
				If[ (Head[nextparse] === Times || Head[nextparse] === Plus) && 
					Length[nextparse] > ZLengthLimit,
					nextparse = Take[nextparse, ZLengthLimit]]; 
				
				newrule = nextparse -> ToExpression[znstr];
				rules = Append[rules, newrule];
				(* Print[zstr, nz2," = ", InputForm[nextparse],";"]; *)
				PrintRule[newrule];
			
				parse = parse /. newrule;
			
				(* find new commmon sub-expressions *)
				(*
				If[SameQ[Parseopt,All],
				parse = Union[parse,ParseCommon[Select[parse,(Length[#] > 2)&],Parseopt]];
					];
				*)
				
				(* parse = Join[Select[parse,ZNumberQ],Select[parse, Depth[#] > 1 &]]; *)
				parse = Select[parse, Depth[#] > 1 &];
								
				n = Length[parse]
				];

		copy = expr //. rules;
		
		(* recursive execution *)
		If[ Depth[copy] > 1,
					temp = ToZTEST[copy, z, nz2, opts];
					copy  = temp[[1]];
					rules = Join[rules, temp[[2]]];
					nz2   = temp[[3]]
			],
			
		Print[copy = expr];
		PrintZ[z,nz2]
		];
		
	{copy, rules, nz2}
	
	];
:[font = input; preserveAspect; endGroup]
ToZTEST[expr_,z_Symbol] := ToZTEST[expr,z,0];
ToZTEST[expr_] := ToZTEST[expr,z];
ToZTEST[expr_,opts___] := ToZTEST[expr,z,0,opts];
ToZTEST[expr_,z_Symbol,opts___] := ToZTEST[expr,z,0,opts];
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
Simple parsing function: 
:[font = input; preserveAspect]
ParseSimple[expr_,res_:{}] := Module[ 

	{n, i,j , res2, restemp, lres},

	res2 = res;
	n = Length[expr];

	If[Depth[expr] > 2,
		
		Do[ res2 = ParseSimple[expr[[i]],res2],
			{i,n}],

		If[!MemberQ[res2,expr], 
			res2 = Join[res2,{expr}]]
		];
	
	Flatten[{res2}]
	];

:[font = input; preserveAspect; endGroup; endGroup]
ToZSimple[expr_,z_,nz_] := Module[ 

	{nz2,parse,copy,rules,n,i,zstr, znstr,temp},
	
	nz2 = nz;
	rules = {};
	zstr = ToString[z];
	parse = ParseSimple[expr];
	
	If[ !SameQ[parse,Flatten[{expr}]],
		n = Length[parse];
			
		Do[
			If[ Depth[parse[[i]]] > 1,
				nz2++;
				znstr = StringJoin[zstr,ToString[nz2]];
				rules = Join[rules, {parse[[i]] -> ToExpression[znstr]}];
				Print[zstr, nz2," = ", InputForm[parse[[i]]],";"]
				]
				,{i,n}
			];
		
		copy = expr //. rules;
		
		If[ Depth[copy] > 1,
					temp = ToZSimple[copy, z, nz2];
					copy  = temp[[1]];
					rules = Join[rules, temp[[2]]];
					nz2   = temp[[3]]
			],
			
		Print[copy = expr];
		PrintZ[z,nz2]
		];
		
	{copy, rules, nz2}
	
	];
ToZSimple[expr_,z_] := ToZSimple[expr,z,0];
ToZSimple[expr_] := ToZSimple[expr,z];	
:[font = subsection; inactive; Cclosed; preserveAspect; startGroup]
Z Compress:
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
Utilities:
:[font = text; inactive; preserveAspect]
not used!
:[font = input; inactive; preserveAspect]
FindFirstZ[rules_,n_,z_:"z"] := Module[

	{zn, rn, out, rtot},

	rtot = Length[rules];
	out = rtot + 1;
	zn = ToExpression[StringJoin[z,ToString[n]]];
	
	For[ rn = n+1, rn <= rtot, rn++,
		If[!FreeQ[rules[[rn,1]],zn],
			out = rn;
			rn = rtot]
			];
			
	out
	];
:[font = input; inactive; preserveAspect]
FindLastZ[rules_,n_,z_:"z"] := Module[

	{zn, rn, out, rtot},

	rtot = Length[rules];
	out = rtot + 1;
	zn = ToExpression[StringJoin[z,ToString[n]]];
	
	For[ rn = rtot, rn > n, rn--,
		If[!FreeQ[rules[[rn,1]],zn],
			out = rn;
			rn = n]
			];
			
	out
	];
:[font = text; inactive; preserveAspect]
C form without "if" expressions:
:[font = input; inactive; preserveAspect]
PrintRule[rule_] := Print["\t",CForm[rule[[2]]]," = ",
						  CForm[rule[[1]]],";"];
:[font = text; inactive; preserveAspect]
C form with "if" expressions and function tranlsations:
:[font = input; preserveAspect]
MathematicaFunctions = {Power, Sqrt,
                        Log, 
                        Sin, Cos, Tan, 
                        ArcSin, ArcCos, ArcTan,
                        Sinh, Cosh, Tanh};
;[s]
2:0,1;11,0;198,-1;
2:1,11,9,Courier,1,10,0,0,0;1,11,9,Courier,3,10,0,0,0;
:[font = input; preserveAspect]
Off[General::spell1];
CFunctions = {pow, sqrt,
              log,
              sin, cos, tan,
              asin, acos, atan,
              sinh, cosh, tanh}; 
On[General::spell1];   
Unprotect[Evaluate[CFunctions], exp]
Map[Clear, CFunctions];
:[font = input; preserveAspect]
Mathematica2CFunctions = Join[{Power[E,x_] -> exp[x],
				Power[x_,2] -> HoldForm[x x],
				Power[x_,-1] -> HoldForm[1/x],
    Power[x_,-2] -> HoldForm[1/(x x)],
				Power[x_,Rational[1,2]] -> HoldForm[sqrt[x]],
				Power[x_,Rational[-1,2]] -> HoldForm[1/sqrt[x]]},
				MapThread[Rule, {MathematicaFunctions, CFunctions}]]
Protect[Evaluate[CFunctions], exp]            
;[s]
4:0,1;11,0;287,1;298,0;370,-1;
2:2,11,9,Courier,1,10,0,0,0;2,11,9,Courier,3,10,0,0,0;
:[font = input; preserveAspect]
ZForm[x_] := Module[

	{xstr},
	
	xstr = ToString[CForm[x]];
	If[StringLength[xstr] > 4 && StringTake[xstr,4] == "-1.*",
		StringDrop[xstr,{2,4}],
		xstr]
	];
:[font = input; preserveAspect]
PrintRule[rule_] := If[ Head[ rule[[1]] ] === If,

		Print["\t",CForm[rule[[2]]]," = ",
				   "(", CForm[rule[[1,1]]],") ? ",
				   ZForm[rule[[1,2]] /. Mathematica2CFunctions //N], " : ", 
				   ZForm[rule[[1,3]] /. Mathematica2CFunctions //N],";"]
		,
		(* else *)
	
		Print["\t",CForm[rule[[2]]]," = ",
				   ZForm[rule[[1]] /. Mathematica2CFunctions //N],";"]
	];
;[s]
7:0,0;155,1;166,0;220,1;231,0;335,1;346,0;372,-1;
2:4,11,9,Courier,1,10,0,0,0;3,11,9,Courier,3,10,0,0,0;
:[font = input; preserveAspect]
Off[General::spell1];
FortranFunctions = {Power, dsqrt,
              dlog,
              dsin, dcos, dtan,
              dasin, dacos, datan,
              dsinh, dcosh, dtanh};
On[General::spell1];   
Unprotect[Evaluate[FortranFunctions], exp]
Map[Clear, FortranFunctions];
:[font = input; preserveAspect]
Mathematica2FortanFunctions = Join[{Power[E,x_] -> dexp[x],
				Power[x_,2] -> HoldForm[x x],
				Power[x_,-1] -> HoldForm[1/x],
    Power[x_,-2] -> HoldForm[1/(x x)],
				Power[x_,Rational[1,2]] -> HoldForm[dsqrt[x]],
				Power[x_,Rational[-1,2]] -> HoldForm[1/dsqrt[x]]},
				MapThread[Rule, {MathematicaFunctions, FortranFunctions}]]
Protect[Evaluate[FortranFunctions], exp]            
;[s]
4:0,1;11,0;295,1;306,0;390,-1;
2:2,11,9,Courier,1,10,0,0,0;2,11,9,Courier,3,10,0,0,0;
:[font = input; preserveAspect]
ZFortranForm[x_] := Module[

	{xstr},
	
	xstr = ToString[FortranForm[x]];
	If[StringLength[xstr] > 4 && StringTake[xstr,4] == "-1.*",
		StringDrop[xstr,{2,4}],
		xstr]
	];
:[font = input; preserveAspect]
PrintFortranRule[rule_] := 	Print["      ",FortranForm[rule[[2]]],"=",
				   ZFortranForm[rule[[1]] /. Mathematica2FortanFunctions //N]];
;[s]
3:0,0;104,1;115,0;139,-1;
2:2,11,9,Courier,1,10,0,0,0;1,11,9,Courier,3,10,0,0,0;
:[font = input; preserveAspect]
MakeRule[znold_,znnew_,z_:"z"] :=

	ToExpression[StringJoin[z,ToString[znold]]] ->
	ToExpression[StringJoin[z,ToString[znnew]]];
	
:[font = input; preserveAspect; endGroup]
MakeRule[znold_,znnew_,z_:"z",Z_:"Z"] :=

	ToExpression[StringJoin[z,ToString[znold]]] ->
	ToExpression[StringJoin[Z,ToString[znnew]]];
	
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
ZCompress driver:
:[font = input; preserveAspect]
GenZList[n_Integer, z_String:"z"] := Module[

	{outtab, i},
	
	outtab = Table[ToString[i], {i,n}];
	Map[ToExpression[StringJoin[z,#]]&,outtab]
	
	];
:[font = input; preserveAspect]
ReturnAtoms[expr_] := Union[Flatten[{Map[ If[Depth[#] > 1, ReturnAtoms[#], #]&,  Apply[List,expr] ]}]]
:[font = input; preserveAspect]
MakeZTab[zout_List, z_:"z"] := Module[ 

	{zused, zout2, rules, rtot, zlist},
	
	{zout2, rules, rtot} = zout;
	
	zlist = GenZList[rtot,z];
	zused = Append[Map[ReturnAtoms, Map[#[[1]]&,rules]],ReturnAtoms[zout2]];		
	zused = Map[Intersection[zlist,#]&,zused];
	zused = Map[Map[ZNum[#,z]&, zused[[#]]]&, Range[Length[zused]]];
		
	(* suppress error messages *)
	Off[Last::nolast];
	ztab = Map[Last[Position[zused,#]][[1]]&, Range[rtot]];
	On[Last::nolast];

	ztab
	
	];
:[font = text; inactive; preserveAspect]
version before revised flag expired searcher:
:[font = input; inactive; preserveAspect]
ZCompress[zout_,z_:"z"] := Module[ 
	
		{ztab, rules, rtot, rn, zn, repl, zout2, zrt, zcount, rs,
		 shift, zrs, zrs1, empties, znmax},
		
		{zout2, rules, rtot} = zout;
		
		(* generate list of "connectivities" *)
		ztab = MakeZTab[zout, z];
				
		(* remove empty rules *)
		empties = Position[ztab,{}];
		If[ Length[empties] > 0,
			Print["Found ", Length[empties], " unused statements\n"];
			rules = Delete[rules, empties];
			
			(* recalculate ztab *)
			ztab = MakeZTab[{zout2,rules,rtot}, z];
			ztab = Map[ (# /. {} -> 1)& , ztab]

			];		
									
		PrintRule[rules[[1]]];
		(* PrintFortranRule[rules[[1]]]; *)
		
		zcount = 0;
		znmax  = 0;
		For[rn = 2, rn <= Length[rules], rn++,
		
			(* z number on the RHS *)
			zrt = ZNumR[rules[[rn]],z];
			
			(* search ztab up to current label number *)
			For[zn = 1, zn < zrt, zn++,
			
				(* expired last flag *)
				If[ ztab[[zn]] <= rn,
				
					zcount++;
					ztab[[zn]] = ztab[[zrt]];
					
					shift = Join[ {MakeRule[zrt, zn, z]},
							Map[ MakeRule[#, #-1,z]& , Range[zrt,rtot - zcount + 1]]];
											                    (* is this right? *)
					ztab = Join[Take[ztab,zrt-1],
								RotateLeft[Take[ztab,{zrt,rtot}]]];
			
					rules = Join[Take[rules, rn-1],
								 Take[rules, {rn, Length[rules]}] /. shift];
								 
					zout2 = zout2 /. shift;

					 zn = rn  (* why is this rn and not zrt ? *)
					(* zn = zrt *) (* exit *)
					] 
				];
				
			znmax = Max[znmax,ZNumR[rules[[rn]],z]];
			(* PrintFortranRule[rules[[rn]]] *)
			PrintRule[rules[[rn]]]
			];
		
		Print[zout2];	
		PrintZ[z, znmax];
			
		{zout2, rules, znmax}
		];
:[font = input; preserveAspect; endGroup; endGroup]
ZCompress[zout_,z_:"z"] := Module[ 
	
		{ztab, rules, rtot, rn, zn, repl, zout2, zrt, zcount, rs,
		 shift, zrs, zrs1, empties, znmax},
		
		{zout2, rules, rtot} = zout;
		
		(* generate list of "connectivities" *)
		ztab = MakeZTab[zout, z];
				
		(* remove empty rules *)
		empties = Position[ztab,{}];
		If[ Length[empties] > 0,
			Print["Found ", Length[empties], " unused statements\n"];
			rules = Delete[rules, empties];
			
			(* recalculate ztab *)
			ztab = MakeZTab[{zout2,rules,rtot}, z];
			ztab = Map[ (# /. {} -> 1)& , ztab]

			];		
									
		(* PrintFortranRule[rules[[1]]]; *)							
		PrintRule[rules[[1]]];
		
		zcount = 0;
		znmax  = 0;
		For[rn = 2, rn <= Length[rules], rn++,
		
			(* z number on the RHS *)
			zrt = ZNumR[rules[[rn]],z];
			
			(* search ztab up to current label number *)
			For[zn = 1, zn < zrt, zn++,
			
				(* expired last flag *)
				If[ ztab[[zn]] <= rn,
				
					zcount++;
					ztab[[zn]] = ztab[[zrt]];
					ztab[[zrt]] = 0;
					
					shift = MakeRule[zrt, zn, z];
			
					rules = Join[Take[rules, rn-1],
								 Take[rules, {rn, Length[rules]}] /. shift];
								 
					zout2 = zout2 /. shift;

					 zn = zrt  (* exit *)
					] 
				];
				
			znmax = Max[znmax,ZNumR[rules[[rn]],z]];	
			(* PrintFortranRule[rules[[rn]]] *)
			PrintRule[rules[[rn]]]
			];
		
		Print[zout2];	
		PrintZ[z, znmax];
			
		{zout2, rules, znmax}
		];
:[font = subsection; inactive; Cclosed; preserveAspect; startGroup]
Operation Count:
:[font = text; inactive; preserveAspect]
Provides number of operations 
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
Utility:
:[font = input; preserveAspect; endGroup]
PowerCount[expr_] := Module[

	 {ff, ps, count, i, pos, expo},

	ff = FullForm[expr];

	(* list of the positions of terms with powers *)	
	ps = Position[ff,Power];	
	count = 0;
	
	Do[
		(* convert position of power to the index 
		   of the actual exponent *)
		pos = ReplacePart[ps[[i]],2,-1];

		(* grab the exponent by stripping the list
		   of its Head *)
		expo = ff[[Delete[pos,0]]];
		
		If[IntegerQ[expo],
			If[ expo > 0,
				expo -= 2,
				If[ expo < 0,
					expo = Abs[expo] - 1]];
			count += expo
			]
		,{i, Length[ps]}];
		
	count
	];
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
Counting:
:[font = input; preserveAspect]
OpCount[expr_] := Module[ 

	{Functs, ff, pow, str, i, minss, commas},

	Functs = {"Cos","Sin","Tan","Exp","Sec","Csc","Cot"};

	commas = Length[StringPosition[ToString[expr],","]];
	pow = PowerCount[expr];
	str = ToString[FullForm[expr]];
	str = StringJoin[str,Table[",", {pow}]];
	
	Do[ str = StringReplace[str, Functs[[i]] -> ","];
			,{i,Length[Functs]}];
	
	(* correct for all *-1 *)
	minss = Length[StringPosition[str,"Times[-1"]];
	
	Length[StringPosition[str,","]] - minss - commas
	];
:[font = input; preserveAspect; endGroup; endGroup]
OpZCount[zout_,prnt_:False] := Module[ {n, ops, rules, temp},

	ops = 0;
	rules = zout[[2]];
	n = Length[rules];
	
	Do[ 
		temp = OpCount[rules[[i,1]]];
		ops = ops + temp;
		If[ prnt, Print[rules[[i,1]]," = ", temp]]
		, {i,n}];
	
	temp = OpCount[zout[[1]]];
	If[prnt, Print[zout[[1]]," = ", temp]];

	ops = ops + temp
	];
:[font = subsection; inactive; Cclosed; preserveAspect; startGroup]
Plot Z:
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
Utilities:
:[font = input; preserveAspect]
ZNum[zn_,z_:"z"] := Module[

	{nz,t1},
	
	t1 = ToString[zn];
	nz = StringLength[z];
	ToExpression[StringTake[t1,nz-StringLength[t1]]]
	];
:[font = input; preserveAspect]
ZNumR[rule_,z_:"z"] := ZNum[rule[[2]],z];
:[font = input; preserveAspect]
Znums[list_,z_:"z"] := Module[

	{i, nlist, t1, zpos},
	
	nlist = {};
	Do[
		t1   = ToString[list[[i]]];
		zpos = Flatten[StringPosition[t1,z]];
		
		If[ Length[zpos] > 0,
			nlist = Append[nlist,ZNum[list[[i]],z]]
			]
			
		, {i,Length[list]}];
	
	nlist
	];
:[font = input; preserveAspect; endGroup]
ZCount[rules_,z_:"z"] := Module[

	{nrules,outtab, i, t1},
	
	nrules = Length[rules];
	outtab = Table[ {}, {nrules}];
	
	Do[
		t1 = ParseAll[rules[[i,1]]];
		outtab[[ZNumR[rules[[i]],z]]] = 
			Join[outtab[[ZNumR[rules[[i]],z]]],
		 		 Znums[t1,z]]
		, {i,nrules}];
	
	outtab
	];
:[font = subsubsection; inactive; Cclosed; preserveAspect; startGroup]
Functions:
:[font = input; preserveAspect; endGroup; endGroup]
PlotZCount[rules_,z_:"z"] := Module[

	{pts, ztab, i},
	
	pts = {};
	ztab = ZCount[rules,z];
	
	Do[
		pts = Join[pts, Map[Point,
				Flatten[Outer[List,{i},ztab[[i]]],1]]]
		,{i, Length[ztab]}];
			
	Graphics[pts]
	];
^*)
