//DEVELOPMENT

#ifndef _FEA_H_
#define _FEA_H_

	#include <stdio.h>
	#include <string.h>
	#include <math.h>
	#include <iostream.h>

	#include "ArrayT.h"
	#include "dArrayT.h"
	#include "dMatrixT.h"

	#include "TensorTransformT.h"

	enum Logic1T 		{ OFF, 		ON 		};
	enum Logic2T 		{ FALSE, 	TRUE 	};
	enum Logic3T 		{ NO, 		YES 	};
	enum AbbrevT 		{ kRow, 	kCol 	};
	enum KineT  		{ kPlaneStrain, kPlaneStress };
	enum MatFillT  	{ kNoZeros, kZeros };
	enum TransposeT { kNonSymetric, kNonSymTranspose, kSymetric };
	enum TimeStepT 	{ kForward_Euler, kBackward_Euler, kCrank_Nicholson, kNewmark, kLax_Friedrich }; 
  enum VersionT   {	kVersionI, kVersionII, kVersionIII, kVersionIV, kVersionV };

	#include "FEA_EquateT.h"
	#include "FEA_dScalarT.h"
	#include "FEA_dMatrixT.h"
	#include "FEA_dVectorT.h"

	#include "FEA_Data_ProcessorT.h"
	#include "FEA_StackT.h"
	#include "FEA_IntegrationT.h"
	#include "FEA_ShapeFunctionT.h"
	#include "FEA_ArraysT.h"
  
#endif

