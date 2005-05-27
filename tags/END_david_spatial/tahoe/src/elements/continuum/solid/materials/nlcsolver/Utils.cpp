/* $Id: Utils.cpp,v 1.6 2004-07-15 08:28:42 paklein Exp $ */
#include "Utils.h"

#include <iostream.h>
#include <ctype.h>

#include "StringT.h"
#include "ifstreamT.h"
#include "LocalArrayT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

namespace Tahoe {

const double sqrt3 = sqrt(3.0);

ifstreamT& OpenExternal(ifstreamT& in, ifstreamT& in2,
			 const char* name)
{
  // check for external file
  char nextchar = in.next_char();
  if (isdigit(nextchar))
      // copy main input file
      return in;
  else
    {
      // open separate file
      StringT filename;
      in >> filename;
      
  		// generate relative path in native format
		filename.ToNativePathName();
		StringT path;
		path.FilePath(in.filename());
		filename.Prepend(path);
      
      OpenExternal(in2, filename, name);
      if (in.skip_comments()) in2.set_marker(in.comment_marker());
      return in2;
    }
}

void OpenExternal(ifstreamT& input, StringT& filename, const char* name)
{
  cout << "\n Utils::OpenExternal: external file for " 
       << name << " : " << filename << endl;
  filename.ToNativePathName();
  input.open(filename);
  if (!input.is_open())
    throwRunTimeError("Utils::OpenExternal: Couldn't open file");
}

void SetStreamPrefs(ostream& stream)
{
  stream.precision(kPrecision);
  stream.setf(ios::showpoint);
  stream.setf(ios::right, ios::adjustfield);
  stream.setf(ios::scientific, ios::floatfield);
}

/* set shape functions and their derivatives at center of element */
/* Na(1 x #nodes), DNa(#sd x #nodes) */
void SetLocalShape_C(dArray2DT& Na, dArray2DT& DNa)
{
  /* dimensions */
  int numsd    = DNa.MajorDim();
  int numnodes = DNa.MinorDim();
  
  /* dimension checks */
  if (numsd != 2 && numsd != 3) 
    throwRunTimeError("Utils::SetLocalShape_C: numsd != 2 or 3");
  if (numnodes != 4 && numnodes != 8)
    throwRunTimeError("Utils::SetLocalShape_C: numnodes != 4 or 8");

  /* initialize */
  Na = 0.0;
  DNa = 0.0;

  if (numsd == 2)
    {
      // nodal coordinates of parent domain 
      double  ra[8] = {-1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 0.0,-1.0};
      double  sa[8] = {-1.0,-1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 0.0};
	
      // coordinates of center of element
      double r = 0.;
      double s = 0.;
  
      // shape funtions and derivatives at center of element
      double* na  = Na(0);
      double* nax = DNa(0);
      double* nay = DNa(1);
  
      for (int lnd = 0; lnd < numnodes; lnd++) 
	{
	  double tempr1 = 1.0 + ra[lnd]*r;
	  double temps1 = 1.0 + sa[lnd]*s;
	  
	  *na++  = 0.25*tempr1*temps1;
	  *nax++ = 0.25*ra[lnd]*temps1;
	  *nay++ = 0.25*tempr1*sa[lnd];
	}
    }
  else       // (nsd == 3)
    {
      // nodal coordinates of parent domain 
      double  ra[] = {-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0};
      double  sa[] = {-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0};
      double  ta[] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0};

      // coordinates of center of element 
      double r = 0.;
      double s = 0.;
      double t = 0.;

      // shape funtion and derivatives at center of element
      double* na  = Na(0);
      double* nax = DNa(0);
      double* nay = DNa(1);
      double* naz = DNa(2);
      
      for (int lnd = 0; lnd < numnodes; lnd++)
	{
	  double tempr1 = 1.0 + ra[lnd]*r;
	  double temps1 = 1.0 + sa[lnd]*s;
	  double tempt1 = 1.0 + ta[lnd]*t;
      
	  *na++  = 0.125*tempr1*temps1*tempt1;
	  *nax++ = 0.125*ra[lnd]*temps1*tempt1;
	  *nay++ = 0.125*tempr1*sa[lnd]*tempt1;
	  *naz++ = 0.125*tempr1*temps1*ta[lnd];
	}
    }
}

// dNa/dE (LDNa) -> dNa/dX (GDNa): coords = Initial_X
// dNa/dE (LDNa) -> dNa/dx (GDNa): coords = Current_X
void ComputeGDNa_C(const LocalArrayT& coords, const dArray2DT& LDNa, dArray2DT& GDNa)
{
  // chain rule Jacobian matrix
  dMatrixT jac(LDNa.MajorDim());
  Jacobian(coords, LDNa, jac); 
  double det = jac.Det();

  // element check
  if (det < kSmall) throwRunTimeError("Utils::ComputeGDNa_C: det(Jac) < kSmall");

  dMatrixT& jac_inv = jac.Inverse();

  // global shape function derivatives
  ChangeOfVariables(jac_inv, LDNa, GDNa);
}

// set shape functions and their derivatives at IPs and the extrapolation matrix
// Na(#IP x #nodes), Na_x[#IP x (#sd x #nodes)], nodal_extrap(#nodes x #IP) 
void SetLocalShape(dArray2DT& Na, ArrayT<dArray2DT>& Na_x, dArray2DT& nodal_extrap)
{
  // fetch some dimensions
  int numint   = Na_x.Length();
  int nsd      = Na_x[0].MajorDim();
  int numnodes = Na_x[0].MinorDim();  // vertex nodes: 4 (Quad), 8 (Hexa)
  
  // check spatial dimensions
  if (nsd != 2 && nsd != 3)
    throwRunTimeError("Utils::SetLocalShape: nsd != 2 or 3");

  //initialize arrays 
  Na = 0.0;
  nodal_extrap = 0.0;
  for (int i = 0; i < Na_x.Length(); i++) Na_x[i] = 0.0;

  // array to hold weights at IPs
  dArrayT weights(numint);
    
  // shape functions and their derivatives in parent domain
  if (nsd == 2) {

    // additional dimension checks
    if (numnodes != 4) 
      throwRunTimeError("Utils::SetLocalShape: numnode != 4 (nsd=2)");
    if (numint != 1 && 
	numint != 4 && 
	numint != 9) 
      throwRunTimeError("Utils::SetLocalShape: numint != 1, 4 or 9 (nsd=2)");

    // integration point coordinates
    double  ra[9] = {-1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 0.0,-1.0, 0.0};
    double  sa[9] = {-1.0,-1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 0.0, 0.0};
    double *xa, *ya;
    double g;
  
    // integration weights
    switch (numint)
      {
      case 1:	
	xa = ra;			  
	ya = sa;			  
	g = 0.0;
	weights[0] = 4.0;
	break;
	
      case 4:
	xa = ra;			  
	ya = sa;			  
	g = 1.0/sqrt3;
	weights[0] = 1.0;
	weights[1] = 1.0;
	weights[2] = 1.0;
	weights[3] = 1.0;
	break;
	
      case 9:
        {
	xa = ra;			  
	ya = sa;			  
	g = sqrt(3.0/5.0);
	double a = 25.0/81.0;
	double b = 40.0/81.0;
	double c = 64.0/81.0;
	weights[0] = a;
	weights[1] = a;
	weights[2] = a;
	weights[3] = a;
	weights[4] = b;
	weights[5] = b;
	weights[6] = b;
	weights[7] = b;
	weights[8] = c;
	break;
        }
	
      default:
	throwRunTimeError("Utils::SetLocalShape: Bad numint (nsd=2)");
      }
  
    // shape functions and derivatives at IPs (using vertex nodes)
    for (int in = 0; in < numint; in++)	
      {	
	double* na  = Na(in);
	double* nax = Na_x[in](0);
	double* nay = Na_x[in](1);
	
	double r = g*xa[in];
	double s = g*ya[in];
	
	for (int lnd = 0; lnd < numnodes; lnd++)
	  {
	    double tempr1 = 1.0 + ra[lnd]*r;
	    double temps1 = 1.0 + sa[lnd]*s;
	    
	    *na++  = 0.25*tempr1*temps1;
	    *nax++ = 0.25*ra[lnd]*temps1;
	    *nay++ = 0.25*tempr1*sa[lnd];
	  }
      }
  }
  else {    // (nsd == 3) 

    /* dimension checks */
    if (numnodes != 8)
      throwRunTimeError("Utils::SetLocalShape: numnodes != 8 (nsd=3)");
    if (numint != 1 && 
	numint != 8 &&
	numint != 27)
      throwRunTimeError("Utils::SetLocalShape: numint != 1, 8 or 27 (nsd=3)");
  
    /* integration point coordinates */
    double  ra[] = {-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0};
    double  sa[] = {-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0};
    double  ta[] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0};
    double  xa_27[27], ya_27[27], za_27[27];
    double  *xa, *ya, *za;
    double  g;
    
    /* integration weights */
    switch (numint)
      {
      case 1:	
	g = 0.0;
	xa = ra;
	ya = sa;
	za = ta;
	weights[0] = 8.0;
	break;
	
      case 8:
	g = 1.0/sqrt3;
	xa = ra;
	ya = sa;
	za = ta;
	weights = 1.0;
	break;
	
      case 27:
        {
	/* coordinates */
	double b1 = sqrt(3.0/5.0);
	double b_1D[3] = {-b1, 0.0, b1}; 
	
	/* weights */
	double w1 = 5.0/9.0;
	double w2 = 8.0/9.0;
	double w_1D[3] = {w1, w2, w1};
	int x_i = 0;
	int y_i = 0;
	int z_i = 0;
	for (int i = 0; i < 27; i++)
	  {
	    xa_27[i]   = b_1D[x_i];
	    ya_27[i]   = b_1D[y_i];
	    za_27[i]   = b_1D[z_i];
	    weights[i] = w_1D[x_i]*w_1D[y_i]*w_1D[z_i];
	    
	    if (++x_i == 3)
	      {
		x_i = 0;
		if (++y_i == 3)
		  {
		    y_i = 0;
		    z_i++;
		  }
	      }
	  }						
	
	xa = xa_27;
	ya = ya_27;
	za = za_27;
	g  = 1.0;		
	break;
        }
	
      default:
	throwRunTimeError("Utils::SetLocalShape: Bad numint (nsd=3)");
      }
    
    // shape functions and derivatives at IPs (using vertex nodes)
    for (int in = 0; in < numint; in++)	
      {
	double* na  = Na(in);
	double* nax = Na_x[in](0);
	double* nay = Na_x[in](1);
	double* naz = Na_x[in](2);
	
	double r = g*xa[in];
	double s = g*ya[in];
	double t = g*za[in];
	
	for (int lnd = 0; lnd < numnodes; lnd++)
	  {
	    double tempr1 = 1.0 + ra[lnd]*r;
	    double temps1 = 1.0 + sa[lnd]*s;
	    double tempt1 = 1.0 + ta[lnd]*t;
	    
	    *na++  = 0.125*tempr1*temps1*tempt1;
	    *nax++ = 0.125*ra[lnd]*temps1*tempt1;
	    *nay++ = 0.125*tempr1*sa[lnd]*tempt1;
	    *naz++ = 0.125*tempr1*temps1*ta[lnd];
	  }
      }
  }

  // nodal extrapolation matrix (least square smoothing)

  dMatrixT NaxNaT     (numnodes);
  dMatrixT NaxNaT_inv (numnodes);
  dArrayT  v1 (numnodes), v2 (numnodes);

  NaxNaT = 0.0;
  for (int i = 0; i < numint; i++)
    {
      Na.RowCopy(i, v1);
      NaxNaT_inv.Outer(v1, v1);   // NaxNaT_inv is dummy here
      NaxNaT.AddScaled(weights[i], NaxNaT_inv);
    }

  NaxNaT_inv = MatrixInversion(NaxNaT);

  for (int i = 0; i < numint; i++)
    {
      Na.RowCopy(i, v1);
      NaxNaT_inv.Multx(v1, v2);
      v2 *= weights[i];
      nodal_extrap.SetColumn(i, v2);
    }

  
  /*NaxNaT = 0.;
  for (int i = 0; i < numnodes; i++)
    for (int j = 0; j < numnodes; j++)
      for (int k = 0; k < numint; k++)
	NaxNaT(i, j) += Na(k, i) * Na(k, j) * weights[k];

  NaxNaT_inv = MatrixInversion(NaxNaT);

  nodal_extrap = 0.;
  for (int i = 0; i < numnodes; i++)
    for (int j = 0; j < numint; j++)
      for (int k = 0; k < numnodes; k++)
	nodal_extrap(i, j) += NaxNaT_inv(i, k) * Na(j, k) * weights[j];
  */ 
  /* 
                double data[64] = {
        (1. + 3.*sqrt(3.))/8.  ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.,
        (1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. - 3.*sqrt(3.))/8.  ,(1. - sqrt(3.))/8.,
        (1. + sqrt(3.))/8.     ,(1. + 3.*sqrt(3.))/8.  ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.,
        (1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. - 3.*sqrt(3.))/8.,
        (1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. + 3.*sqrt(3.))/8.  ,(1. + sqrt(3.))/8.,
        (1. - 3.*sqrt(3.))/8.  ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.,
        (1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. + 3.*sqrt(3.))/8.,
        (1. - sqrt(3.))/8.     ,(1. - 3.*sqrt(3.))/8.  ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.,
        (1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. - 3.*sqrt(3.))/8.  ,(1. - sqrt(3.))/8.,
        (1. + 3.*sqrt(3.))/8.  ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.,
        (1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. - 3.*sqrt(3.))/8.,
        (1. + sqrt(3.))/8.     ,(1. + 3.*sqrt(3.))/8.  ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.,
        (1. - 3.*sqrt(3.))/8.  ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.,
        (1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. + 3.*sqrt(3.))/8.  ,(1. + sqrt(3.))/8.,
        (1. - sqrt(3.))/8.     ,(1. - 3.*sqrt(3.))/8.  ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.,
        (1. + sqrt(3.))/8.     ,(1. - sqrt(3.))/8.     ,(1. + sqrt(3.))/8.     ,(1. + 3.*sqrt(3.))/8.};

                        dMatrixT smooth(8,8,data);
      //                  nodal_extrap = smooth;
        for (int i = 0; i < numnodes; i++)
	   for (int j = 0; j < numint; j++)
               nodal_extrap(i, j) = smooth(i, j);
 */ 
}

// dNa/dE (LDNa) -> dNa/dX (GDNa): coords = Initial_X
// dNa/dE (LDNa) -> dNa/dx (GDNa): coords = Current_X
void ComputeGDNa(const LocalArrayT& coords, const ArrayT<dArray2DT>& LDNa, 
		 ArrayT<dArray2DT>& GDNa)
{
  // fetch some dimensions
  int numint = LDNa.Length();
  int nsd    = LDNa[0].MajorDim();

  // allocate space for Jacobian
  dMatrixT jac(nsd);

  // loop over all integration points
  for (int i = 0; i < numint; i++)
    {
      // compute Jacobian matrix
      Jacobian(coords, LDNa[i], jac); 
      double det = jac.Det();
      
      // element check
      if (det < kSmall) throwRunTimeError("Utils::ComputeGDNa: det(Jac) < ksmall");

      dMatrixT& jac_inv = jac.Inverse();

      // global shape function derivatives
      ChangeOfVariables(jac_inv, LDNa[i], GDNa[i]);
    }
}

void Jacobian(const LocalArrayT& nodal, const dArray2DT& DNa, dMatrixT& jac)
{
  // fetch dimensions
  int nsd      = DNa.MajorDim();
  int numnodes = DNa.MinorDim();

  // work spaces
  dArrayT va(nsd);
  dArrayT vb(nsd);
  dMatrixT dummy(nsd);

  // jacobian
  jac = 0.;
  for (int i = 0; i < numnodes; i++)
    {
      DNa.ColumnCopy(i, vb);
      for (int j = 0; j < nsd; j++) va[j] = nodal(i,j);
      dummy.Outer(va, vb);
      jac += dummy;
    }
}

/* dNa/dE -> dNa/dX : tensor = Jac^(-1) */
/* dNa/dX -> dNa/dx : tensor = F^(-1)   */
void ChangeOfVariables(const dMatrixT& tensor, const dArray2DT& DNa_known,
		      dArray2DT& DNa_unknown)
{
  // fetch dimensions
  int nsd      = DNa_known.MajorDim();
  int numnodes = DNa_known.MinorDim();

  // work space arrays
  dArrayT va(nsd);
  dArrayT vb(nsd);

  // change of variables in derivatives
  for (int i = 0; i < numnodes; i++)
    {
      DNa_known.ColumnCopy(i, va);
      tensor.MultTx(va, vb);
      DNa_unknown.SetColumn(i, vb);
    }
}

// nodal(#nodes x (#sd x #sd)), interp(#IP x (#sd x #sd)), Na(#IP x #nodes)
void Interpolate(const ArrayT<dMatrixT>& nodal, ArrayT<dMatrixT>& interp,
		 const dArray2DT& Na)
{
  // fetch dimensions
  int numint   = Na.MajorDim();
  int numnodes = Na.MinorDim();

  // interpolate to all integration points
  for (int i = 0; i < numint; i++)
    {
      interp[i] = 0.0;
      for (int j = 0; j < numnodes; j++)
	interp[i].AddScaled(Na(i, j), nodal[j]);
    }
}

// ipvalues(#IP x (#sd x #sd)), extrap(#nodes x (#sd x #sd)), Ma(#nodes x #IP)
void Extrapolate(const ArrayT<dMatrixT>& ipvalues, ArrayT<dMatrixT>& extrap,
		 const dArray2DT& Ma)
{
  // fetch dimensions
  int numnodes = Ma.MajorDim();
  int numint   = Ma.MinorDim();

  // extrapolate to all nodal points
  for (int i = 0; i < numnodes; i++)
    {
      extrap[i] = 0.0;
      for (int j = 0; j < numint; j++)
	extrap[i].AddScaled(Ma(i, j), ipvalues[j]);
    }
}

void LUDecomposition(dMatrixT& a, const int& n, iArrayT& indx, double& d)
{
   const double TINY = 1.e-20;

   int i, imax, j, k;
   double big, dum, sum, temp;
   dArrayT vv(n);

   d = 1.;
   for (i=0; i<n; i++) {
      big = 0.0;
      for (j=0; j<n; j++) 
         if ((temp=fabs(a(i,j))) > big) big = temp; 
      if (big == 0.0)
	throwRunTimeError("Utils::LUDecomposition: Singular matrix");
      vv[i] = 1.0/big;
   }
   for (j=0; j<n; j++) {
      for (i=0; i<j; i++) {
         sum = a(i,j);
         for (k=0; k<i; k++) sum -= a(i,k)*a(k,j);
         a(i,j) = sum;
      }
      big = 0.0;
      for (i=j; i<n; i++) {
         sum = a(i,j);
         for (k=0; k<j; k++)
            sum -= a(i,k)*a(k,j);
         a(i,j) = sum;
         if ( (dum=vv[i]*fabs(sum)) >= big ) {
            big = dum;
            imax = i;
         }    
      }
      if (j != imax) {
         for (k=0; k<n; k++) {
            dum = a(imax,k);
            a(imax,k) =a(j,k);
            a(j,k) = dum;
         }
         d = -(d);
         vv[imax] = vv[j];
      }
      indx[j] = imax;
      if (a(j,j) == 0.0) a(j,j) = TINY;
      if (j!= (n-1)) {
         dum = 1.0/(a(j,j));
         for (i=j+1; i<n; i++) a(i,j) *= dum;
      }
   }
}
       
void LUBackSubstitution(dMatrixT& a, const int& n, const iArrayT& indx, dArrayT& b)
{
   int i, ii, ip, j;
   double sum;

   ii = -1;
   for (i=0; i<n; i++) {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if (ii > -1)
         for (j=ii; j<=i-1; j++) sum -= a(i,j)*b[j];
      else if (sum) ii = i;
      b[i] = sum;
   }
   for (i=n-1; i>=0;i--) {
      sum = b[i];
      for (j=i+1; j<n; j++) sum -= a(i,j)*b[j];
      b[i] = sum/a(i,i);
   }
}
         
dMatrixT MatrixInversion(dMatrixT& a)
{
   int i, j, n;
   double d;
   
   if ( (n=a.Rows()) != a.Cols() )
     throwRunTimeError("Utils::MatrixInversion(a): a.Rows != a.Cols");

   iArrayT indx(n);
   dArrayT col(n);
   dMatrixT ainv(n,n);

   LUDecomposition(a, n, indx, d);
   for (j=0; j<n; j++){
      for (i=0; i<n; i++) col[i] = 0.0;
      col[j] = 1.0;
      LUBackSubstitution(a, n, indx, col);
      for(i=0; i<n; i++) ainv(i,j) = col[i];
   }
   return ainv;
}

void writeMessage(const char* msg)
{
  cout << "\n tahoe::material message";

  if (msg)
     cout << ":\n" << " " << msg << '\n';
  else
     cout << ". (no message)\n";
}

void writeWarning(const char* msg)
{
  cout << "\n tahoe::material warning";

  if (msg)
     cout << ":\n" << " " << msg << '\n';
  else
     cout << ". (no warning)\n";
}

void throwRunTimeError(const char* msg)
{
  cout << "\n tahoe::material run time error";

  if (msg)
     cout << ":\n" << " " << msg << '\n';
  else
     cout << ".\n";

  cout << " exiting tahoe with general fail status\n";

  throw ExceptionT::kGeneralFail;
}

void throwMemoryError(const char* msg)
{
  cout << "\n tahoe::material memory error";

  if (msg)
     cout << ":\n" << " " << msg << '\n';
  else
     cout << ".\n";

  cout << " exiting tahoe with memory fail status\n";

  throw ExceptionT::kOutOfMemory;
}

} // namespace Tahoe
