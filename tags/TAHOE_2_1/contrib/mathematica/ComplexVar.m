BeginPackage["ComplexVar`"]

MakeReal::usage="MakeReal[var sequence] sets the imaginary part of the listed \
variables to zero, so that Re[var]=var."
FMakeReal::usage="FMakeReal definition."
Bar::usage="Bar[expr] takes the complex conjugate of expr, using of MakeReal \
variables where possible."				
ReQ::usage="ReQ is a test to speed simplification with MakeReal variables.  \
ReQ[x] = True if x has been declared Real, if x is a Real combination of \
variables declared Real, if x is an FReQ function with an ReQ argument, or if \
a Number is Real by default."					
FReQ::usage="FReQ is a flag to indicate that a function returns real values \
for real-valued arguments."
Integ::usage="Integ[f, x] is function that circumvents some of special \
integration features of Mathematica dealing with integrals of I. Integ[f, x] \
gives the same result as Integrate[f,x], treating I as a regular constant.  \
Integ[f, {x, xmin, xmax}] gives the definite integral."
ComplexVar::warning = "This package modifies the default definitions for the \
following Protected symbols: {Re, Im, Log, Exp, Abs, Pi, E, Real, Integer, \
[Trig Functions], [Hyperbolic Functions]}."
ComplexVar::info = "Please report any questions or problems with this package \
to paklein@sandia.gov ."
TrigExpand::usage = "TrigExpand is an option for the package ComplexVar.  \
With TrigExpand -> True, terms involving trigonometric and hyperbolic \
functions with arguments involving only MakeReal variables and I, the \
Sqrt[-1], are treated as rational functions of exponentials by Re and Im.  \
TrigExpand -> False causes these functions to be treated as indivisible \
objects."
PolyExpand::usage = "PolyExpand is an option for the package ComplexVar.  \
With PolyExpand -> False, polynomials of the form (x + I y)^n, where x and y \
are in general complex, are left unevaluated by Re or Im if n > 2 or n < -1."
TrigPowerExpand::usage ="TrigPowerExpand is an option for the package \
ComplexVar. With TrigPowerExpand -> False, powers of Trig functions greater \
than 2 are not evaluated by Re and Im.  With TrigPowerExpand -> True, Trig \
functions with arguments whose Re and Im parts are ReQ are expanded as \
exponentials and evaluated.  Note: TrigExpand -> False supersedes any setting \
of TrigPowerExpand."

Begin["`Private`"]

(* Message[ComplexVar::warning]; *)
Message[ComplexVar::info];    

(* Options *)

Options[ComplexVar] = {TrigExpand -> True, PolyExpand ->True, TrigPowerExpand \
-> False}

(* Special Function Groups *)

TrigFunctions = {Sin, Cos, Tan, Csc, Sec, Cot};
InverseTrigFunctions = {ArcSin, ArcCos, ArcTan};
HypFunctions  = {Sinh, Cosh, Tanh, Csch, Sech, Coth};

protected = Unprotect[Re, Im, Log, Exp, Abs, Arg, Pi, E, Real,
																					 Integer, MakeReal, Bar, FReQ, ReQ, Integ,
																					 Evaluate[TrigFunctions], Evaluate[HypFunctions], \
Evaluate[InverseTrigFunctions],
																					 ComplexVar, TrigExpand, PolyExpand, TrigPowerExpand]

(* Function to set the attributes of Real variables *)

NoIm[x_] := Module[ {},
	ReQ[x] ^= True;
	x/: Re[x] = x;
	x/: Im[x] = 0];

MakeReal[{x__}]:= Map[NoIm, {x}];       (* For lists of variables *)
MakeReal[x___] := MakeReal[{x}];        (* For a sequence of not enclosed in \
brackets *)

(* ReQ is a fast flag for marking Real variables *)

ReQ[I]		                = False
ReQ[-I]                 = False
ReQ[ E ]               ^= True
ReQ[ Pi ]        	     ^= True
ReQ[ x_Real ]          ^= True
ReQ[ x_Integer ]       ^= True
ReQ[ x_ + y_ ]          = ReQ[x] && ReQ[y]
ReQ[ x_ y_ ]            = ReQ[x] && ReQ[y]
ReQ[ x_/y_ ]            = ReQ[x] && ReQ[y]
ReQ[ x_^n_ ]	           = ReQ[x] && ReQ[n]
ReQ[ Rational[x_,y_] ]  = ReQ[x] && ReQ[y]
ReQ[ Re[x_] ]           = True
ReQ[ Im[x_] ]           = True
ReQ[ f_?FReQ[x_?ReQ] ]  = True
ReQ[_]                  = False   (* Default False for all others *)

(* Real-Valued Functions given real arguments *)

FNoIm[f_] := f /: FReQ[f] = True;
 
FMakeReal[{x__}] := Map[ FNoIm, {x} ];
FMakeReal[x___]  := FMakeReal[{x}];

FMakeReal[TrigFunctions];
FMakeReal[HypFunctions];
FMakeReal[InverseTrigFunctions];
FReQ[Bar] ^= True
FReQ[Log] ^= True
FReQ[_] = False

Re[ I ] = 0

Re[ Re[x_] ] = Re[x]
Im[ Re[x_] ] = 0 

Re[x_ + y_] = Re[x] + Re[y]
Im[x_ + y_] = Im[x] + Im[y]

Re[x_ y_] = Re[x] Re[y] - Im[x] Im[y]
Im[x_ y_] = Re[x] Im[y] + Re[y] Im[x] 


(* also handles powers of FReQ functions with ReQ arguments *)

Re[ x_?ReQ^n_?ReQ ] = x^n
Im[ x_?ReQ^n_?ReQ ] = 0

(* Suppress evaluations of functions in general *)

Re[x_ x_] := Re[x]^2-Im[x]^2 /; !FReQ[Head[x]]
Im[x_ x_] := 2 Re[x] Im[x]   /; !FReQ[Head[x]]

Re[ Exp[x_] ] = Exp[ Re[x] ] Cos[ Im[x] ]
Im[ Exp[x_] ] = Exp[ Re[x] ] Sin[ Im[x] ]

(* Real-valued functions with real valued arguments *)

Re[ f_?FReQ[x_?ReQ] ] = f[x]
Im[ f_?FReQ[x_?ReQ] ] = 0       


(* Trig functions with non-real valued arguments *)

(* test *)

TrigQ[x_] := MemberQ[ Join[TrigFunctions, HypFunctions], x ]

(* PolyCheck returns True if the argument x is of the form a + I b *)
(* where a and b are ReQ. *)

PolyCheck[x_] := If[ TrigExpand /. Options[ComplexVar],
           ReQ[ExpandDenominator[ Together[I (x - Bar[x])] ]],
           False
                   ] /; !AtomQ[x]

(* Common occurrences of TrigQ functions with PolyCheck arguments *)

Re[Cos[x_?PolyCheck] ]  := Cos[Re[x]] Cosh[Im[x]]
Im[Cos[x_?PolyCheck] ]  :=-Sin[Re[x]] Sinh[Im[x]]

Re[Sin[x_?PolyCheck] ]  := Sin[Re[x]] Cosh[Im[x]]
Im[Sin[x_?PolyCheck] ]  := Cos[Re[x]] Sinh[Im[x]]

Re[Cosh[x_?PolyCheck] ] := Cosh[Re[x]] Cos[Im[x]]
Im[Cosh[x_?PolyCheck] ] := Sinh[Re[x]] Sin[Im[x]]

Re[Sinh[x_?PolyCheck] ] := Sinh[Re[x]] Cos[Im[x]]
Im[Sinh[x_?PolyCheck] ] := Cosh[Re[x]] Sin[Im[x]]

(* all others *)

Re[f_?TrigQ[x_?PolyCheck]] := Module[ {fxtemp},
     fxtemp = ExpandAll[ f[x], Trig->True ];
     Re[ fxtemp ]
     ];
     
Im[f_?TrigQ[x_?PolyCheck]] := Module[ {fxtemp},
     fxtemp = ExpandAll[ f[x], Trig->True ];
     Im[ fxtemp ]
     ];


(* Higher Order Polynomials *)

(* AtomQ test stops evaluation since Expand returns same. *)
(* Algorithm for expanding powers does not work for powers of *)
(* functions with Heads, ie. FReQ, because the expansion produces *)
(* the same argument, therefore causing an infinite loop *)

Re[x_^n_Integer?Positive] := Module[ {xx},
	xx = Expand[x^n];
	Re[ xx ]
	] /; !AtomQ[x] && !FReQ[Head[x]] && (PolyExpand /. Options[ComplexVar])
	
Im[x_^n_Integer?Positive] := Module[ {xx},
	xx = Expand[x^n];
	Im[ xx ]
	] /; !AtomQ[x] && !FReQ[Head[x]] && (PolyExpand /. Options[ComplexVar])

Re[ 1/x_ ] := Module[ {Rex, Imx},
									Rex = Re[x];
									Imx = Im[x];
									Rex/(Rex^2 + Imx^2)
									];
									
Im[ 1/x_ ] := Module[ {Rex, Imx},
									Rex = Re[x];
									Imx = Im[x];
									-Imx/(Rex^2 + Imx^2)
									];

Re[x_^n_Integer?Negative] := Module[ {NumxBar, Numx, Denx},
							  Numx = Expand[x^-n];
	        NumxBar= Bar[Numx];
	        Denx = Expand[NumxBar Numx];
	        Re[ NumxBar ]/Denx
	        ] /; !AtomQ[x] && !FReQ[Head[x]] && (PolyExpand /. \
Options[ComplexVar])
	
Im[x_^n_Integer?Negative] := Module[ {NumxBar, Numx, Denx},
	        Numx = Expand[x^-n];
	        NumxBar= Bar[Numx];
	        Denx = Expand[NumxBar Numx];
	        Im[ NumxBar ]/Denx
	        ] /; !AtomQ[x] && !FReQ[Head[x]] && (PolyExpand /. \
Options[ComplexVar])
		

(* Powers of TrigQ functions with PolyCheck arguments *)

(* Squared *)

Re[ Power[f_?TrigQ[x_?PolyCheck], 2] ] := Re[ f[x] ]^2-Im[ f[x] ]^2 
Im[ Power[f_?TrigQ[x_?PolyCheck], 2] ] := 2 Re[ f[x] ] Im[ f[x] ]   

(* Higher Powers *)

Re[ f_?TrigQ[x_?PolyCheck]^n_Integer ] := Module[ {ftemp},
	    ftemp = ComplexExpand[ f[x]^n ];
	    Re[ftemp]
	    ]  /; (TrigPowerExpand /. Options[ComplexVar])

     
Im[ f_?TrigQ[x_?PolyCheck]^n_Integer ] := Module[ {ftemp},
	    ftemp = ComplexExpand[ f[x]^n ];
	    Im[ftemp]
	    ]  /; (TrigPowerExpand /. Options[ComplexVar])

		
(* Other Rules *)	
		
Abs[ x_ ] = Sqrt[x Bar[x]]
Im[Abs[_]] = 0

Arg/: Im[Arg[_]] = 0


(* Some more rules for simplification *)
(* only split Logs on Re/Im *)

Re[ Log[x_ y_] ] = Re[Log[x]] + Re[Log[y]]
Im[ Log[x_ y_] ] = Im[Log[x]] + Im[Log[y]]

Re[ Log[x_/y_] ] = Re[Log[x]] - Re[Log[y]]
Im[ Log[x_/y_] ] = Im[Log[x]] - Im[Log[y]]

Log/: Re[ Log[x_?ReQ] ] = Log[x]
(* Log/: Im[ Log[x_?ReQ] ] = 0 *) (* not true unless x > 0 *)

Log/: Log[ Exp[x_] ] = x
Log/: Exp[ Log[x_] ] = x

(* ArcTan w/ 2 arguments *)
ArcTan/: Re[ArcTan[x_,y_]] := ArcTan[x,y] /; ReQ[x/y]
ArcTan/: Im[ArcTan[x_,y_]] := 0 										/; ReQ[x/y]

(* Rules for Conjugation *)

Bar/: Bar[ I ]          :=-I
Bar/: Bar[-I ]          := I
Bar/: Bar[ x_?ReQ ]   		:= x
Bar/: Bar[ x_?NumberQ ] := Conjugate[x]
Bar/: Bar[ Bar[x_] ]    := x
Bar/: Bar[ x_ ]         := Map[ Bar, x ] /; !AtomQ[x]


(* Integration Functions *)

Integ[f_, x_] := Module[ {Itemp, ftemp, fInt},
	ftemp = f /. Complex[0,y_] -> y Itemp;
	fInt = Integrate[ ftemp, x ];
	fInt /. Itemp -> I
	];

Integ[f_, {x_, xmin_, xmax_}] := Module [ {fInt, fdef, xtemp},
	fInt = Integ[f, x];
	fdef = (fInt /. x->xmax) - (fInt /. x -> xmin)
	];

Protect[Evaluate[protected]]
Protect[MakeReal, Bar, ReQ, Integ, ComplexVar, PolyExpand, TrigExpand]

End[]

EndPackage[]
